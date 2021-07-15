//
//  particle.c
//  softBidispersityMetricJamming
//
//  Created by Zhaoyu Xie on 4/8/20.
//  Copyright © 2020 Zhaoyu Xie. All rights reserved.
//

#include "particle.h"
#include "surface.h"
#include "energy.h"
#include "hessian.h"
#include <math.h>
#include <stdlib.h>

void partitionOneParticle(particle *p, int nPart[], double partBoundX[], double partBoundY[], double partBoundZ[]) {
    double x = (p->rr)*sin(p->th)*cos(p->ph);
    double y = (p->rr)*sin(p->th)*sin(p->ph);
    double z = (p->rr)*cos(p->th);
    int i=0, j=0, k=0;
    while(i<nPart[0]&&x>=partBoundX[i])
        i++;
    p->coord[0] = i-1;
    while(j<nPart[1]&&y>=partBoundY[j])
        j++;
    p->coord[1] = j-1;
    while(k<nPart[2]&&z>=partBoundZ[k])
        k++;
    p->coord[2] = k-1;
}

void coordinateTransform(double x, double y, double z, double a, double b, double *th, double *ph) {
    *th = acos(z/a);
    if(*th<MACHINE_EPSILON || PI-(*th)<MACHINE_EPSILON)
        *ph = 0;
    else
        if(y>=0)
            *ph = acos(x/b/sin(*th));
        else
            *ph = 2*PI-acos(x/b/sin(*th));
}

void addForces(int count, double nm[][3], particle p[], int np, double epsilon) {
    for(int i=0; i<np; i++) {
        p[i].forceTh = 0;
        p[i].forcePh = 0;
    }
    for(int i=0; i<np; i++)
        for(int j=0; j<np; j++) {
            if(j==i||abs(p[i].coord[0]-p[j].coord[0])>1||abs(p[i].coord[1]-p[j].coord[1])>1||abs(p[i].coord[2]-p[j].coord[2])>1)
                continue;
            p[i].forceTh -= Vth(count, nm, &(p[i]), &(p[j]), epsilon);
            p[i].forcePh -= Vph(count, nm, &(p[i]), &(p[j]), epsilon);
        }
}

void moveSystem(int count, double nm[][3], particle p[], int np, double cForce[], double dt, int nPart[], double partBoundX[], double partBoundY[], double partBoundZ[]) {
    double volume = volumeFromCoefficients(count, nm);
    for(int i=0; i<count; i++)
        nm[i][2] += cForce[i]*dt;
    projectOntoConstraint(count, nm, volume);
    for(int i=0; i<np; i++) {
        p[i].th += p[i].forceTh*dt;
        p[i].ph += p[i].forcePh*dt;
        p[i].rr = radiusFromCoefficients(count, nm, p[i].th, p[i].ph);
        partitionOneParticle(&(p[i]), nPart, partBoundX, partBoundY, partBoundZ);
    }
}

void deformSurface(int count, double nm[][3], particle p[], int np, double cForce[], double dt, double volume, int nPart[], double partBoundX[], double partBoundY[], double partBoundZ[]) {
    for(int i=0; i<count; i++)
        nm[i][2] += cForce[i]*dt;
/*
    for(int c=0; c<count; c++)
        printf("%d %d %.15lf\n", (int)(nm[c][0]), (int)(nm[c][1]), nm[c][2]);
*/
    projectOntoConstraint(count, nm, volume);
/*
    for(int c=0; c<count; c++)
        printf("%d %d %.15lf\n", (int)(nm[c][0]), (int)(nm[c][1]), nm[c][2]);
*/    
    for(int i=0; i<np; i++) {
        p[i].rr = radiusFromCoefficients(count, nm, p[i].th, p[i].ph);
        partitionOneParticle(&(p[i]), nPart, partBoundX, partBoundY, partBoundZ);
    }
}

double energyForTimestep(int count, double nm[][3], particle p[], int np, double epsilon, double cForce[], double dt, int nPart[], double partBoundX[], double partBoundY[], double partBoundZ[]) {
    double coeff[count];
    for(int i=0; i<count; i++)
        coeff[i] = nm[i][2];
    for(int i=0; i<np; i++) {
        p[i].rrOld = p[i].rr;
        p[i].thOld = p[i].th;
        p[i].phOld = p[i].ph;
        for(int j=0; j<3; j++)
            p[i].oldCoord[j] = p[i].coord[j];
    }
    moveSystem(count, nm, p, np, cForce, dt, nPart, partBoundX, partBoundY, partBoundZ);
    double energy = totalEnergy(count, nm, p, np, epsilon);
    for(int i=0; i<count; i++)
        nm[i][2] = coeff[i];
    for(int i=0; i<np; i++) {
        p[i].rr = p[i].rrOld;
        p[i].th = p[i].thOld;
        p[i].ph = p[i].phOld;
        for(int j=0; j<3; j++)
            p[i].coord[j] = p[i].oldCoord[j];
    }
    return energy;
}

double goldenSearch(int count, double nm[][3], particle p[], int np, double epsilon, double cForce[], int nPart[], double partBoundX[], double partBoundY[], double partBoundZ[]){
    double phi = 1.618033988749895;//goldenRatio+1
    double w = 0.3819660112501051;//1-goldenRatio
    double tol = MACHINE_EPSILON;
    double currentEnergy = totalEnergy(count, nm, p, np, epsilon);
    
    double maxBracket = fabs(nm[0][2]/cForce[0]);
    for(int i=1; i<count; i++)
        if(fabs(nm[i][2]/cForce[i])<maxBracket)
            maxBracket = fabs(nm[i][2]/cForce[i]);
    double initialBracketGuess = 0.001*maxBracket;
    double bracketUpperlimit = 0.1*maxBracket;
    
    double dtA, dtB, enA, enB;
    dtA = initialBracketGuess;
    enA = energyForTimestep(count, nm, p, np, epsilon, cForce, dtA, nPart, partBoundX,partBoundY, partBoundZ);
    dtB = dtA*phi;
    enB = energyForTimestep(count, nm, p, np, epsilon, cForce, dtB, nPart, partBoundX,partBoundY, partBoundZ);
    do {
        if(enA >= currentEnergy) {
            if(dtA>tol||(enA-currentEnergy)/currentEnergy>tol){
                if(enA == enB)
                    return 0;
                dtB = dtA;
                enB = enA;
                dtA = dtA/phi;
                enA = energyForTimestep(count, nm, p, np, epsilon, cForce, dtA, nPart, partBoundX,partBoundY, partBoundZ);
            } else
                return 0;
        } else if(enB <= enA) {
            if(dtB>bracketUpperlimit)
                return bracketUpperlimit;
            dtA = dtB;
            enA = enB;
            dtB = dtB*phi;
            enB = energyForTimestep(count, nm, p, np, epsilon, cForce, dtB, nPart, partBoundX,partBoundY, partBoundZ);
        }
    } while(enA>=currentEnergy || enB<=enA);
    //Now the initial bracket is (0, dtA, dtB). Find the minimum next.
    double dtx, dty;
    double enx, eny;
    dtx = dtA;
    enx = enA;
    dtA = 0;
    enA = currentEnergy;
    //(dtA, dtx, dty, dtB)
    do {
        if(dtB-dtx > dtx-dtA) {
            dty = dtx + w*(dtB-dtx);
            eny = energyForTimestep(count, nm, p, np, epsilon, cForce, dty, nPart, partBoundX,partBoundY, partBoundZ);
        }
        else {
            dty = dtx;
            eny = enx;
            dtx = dtx - w*(dtx-dtA);
            enx = energyForTimestep(count, nm, p, np, epsilon, cForce, dtx, nPart, partBoundX,partBoundY, partBoundZ);
        }
        if(enx<eny) {
            dtB = dty;
            enB = eny;
        } else {
            dtA = dtx;
            enA = enx;
            dtx = dty;
            enx = eny;
        }
    } while(dtB-dtA > tol*(dtx+dty));
    //printf("dt: %lf\n", dtx);
    return dtx;
}

void calculateConjugateForces(int count, double nm[][3], particle p[], int np, double epsilon) {
    for(int i=0; i<np; i++) {
        p[i].forceTh = 0;
        p[i].forcePh = 0;
    }
    addForces(count, nm, p, np, epsilon);
    double gammaNum = 0;
    double gammaDenom = 0;
    for(int i=0; i<np; i++) {
        gammaNum += p[i].forceTh*p[i].forceTh + p[i].forcePh*p[i].forcePh;
        gammaDenom += p[i].forceThOld*p[i].forceThOld + p[i].forcePhOld*p[i].forcePhOld;
    }
    for(int i=0; i<np; i++) {
        p[i].conjGradTh = p[i].forceTh + gammaNum/gammaDenom*p[i].conjGradTh;
        p[i].conjGradPh = p[i].forcePh + gammaNum/gammaDenom*p[i].conjGradPh;
        p[i].forceThOld = p[i].forceTh;
        p[i].forcePhOld = p[i].forcePh;
    }
}

double particleEnergyForTimestep(int count, double nm[][3], particle p[], int np, double epsilon, double dt, int nPart[], double partBoundX[], double partBoundY[], double partBoundZ[]) {
    for(int i=0; i<np; i++) {
        p[i].rrOld = p[i].rr;
        p[i].thOld = p[i].th;
        p[i].phOld = p[i].ph;
        for(int j=0; j<3; j++)
            p[i].oldCoord[j] = p[i].coord[j];
    }
    for(int i=0; i<np; i++) {
        p[i].th += p[i].conjGradTh*dt;
        p[i].ph += p[i].conjGradPh*dt;
        p[i].rr = radiusFromCoefficients(count, nm, p[i].th, p[i].ph);
        partitionOneParticle(&(p[i]), nPart, partBoundX, partBoundY, partBoundZ);
    }
    double energy = totalParticleEnergy(p, np, epsilon);
    for(int i=0; i<np; i++) {
        p[i].rr = p[i].rrOld;
        p[i].th = p[i].thOld;
        p[i].ph = p[i].phOld;
        for(int j=0; j<3; j++)
            p[i].coord[j] = p[i].oldCoord[j];
    }
    return energy;
}

double goldenSearchParticles(int count, double nm[][3], particle p[], int np, double epsilon, int nPart[], double partBoundX[], double partBoundY[], double partBoundZ[]){
    double phi = 1.618033988749895;//goldenRatio+1
    double w = 0.3819660112501051;//1-goldenRatio
    double tol = MACHINE_EPSILON;
    double currentEnergy = totalParticleEnergy(p, np, epsilon);
    
    double maxForce = 0;
    for(int i=1; i<np; i++) {
        if(fabs(p[i].conjGradTh)>maxForce)
            maxForce = fabs(p[i].conjGradTh);
        if(fabs(p[i].conjGradPh)>maxForce)
            maxForce = fabs(p[i].conjGradPh);
    }
    if(maxForce == 0)
        return 1; //any finite number
    double initialBracketGuess = 0.01/maxForce;
    
    double dtA, dtB, enA, enB;
    dtA = initialBracketGuess;
    enA = particleEnergyForTimestep(count, nm, p, np, epsilon, dtA, nPart, partBoundX,partBoundY, partBoundZ);
    dtB = dtA*phi;
    enB = particleEnergyForTimestep(count, nm, p, np, epsilon, dtB, nPart, partBoundX,partBoundY, partBoundZ);
    if(enA == 0)
        return dtA;
    do {
        if(enA >= currentEnergy) {
            if(dtA >= tol){
                dtB = dtA;
                enB = enA;
                dtA = dtA/phi;
                enA = particleEnergyForTimestep(count, nm, p, np, epsilon, dtA, nPart, partBoundX,partBoundY, partBoundZ);
                if(enA == 0)
                    return dtA;
            } else
                return 0;
        } else if(enB <= enA) {
            dtA = dtB;
            enA = enB;
            dtB = dtB*phi;
            enB = particleEnergyForTimestep(count, nm, p, np, epsilon, dtB, nPart, partBoundX,partBoundY, partBoundZ);
            if(enA == 0)
                return dtA;
        }
    } while(enA>=currentEnergy || enB<=enA);
    //Now the initial bracket is (0, dtA, dtB). Find the minimum next.
    double dtx, dty;
    double enx, eny;
    dtx = dtA;
    enx = enA;
    dtA = 0;
    enA = currentEnergy;
    //(dtA, dtx, dty, dtB)
    do {
        if(dtB-dtx > dtx-dtA) {
            dty = dtx + w*(dtB-dtx);
            eny = particleEnergyForTimestep(count, nm, p, np, epsilon, dty, nPart, partBoundX,partBoundY, partBoundZ);
        }
        else {
            dty = dtx;
            eny = enx;
            dtx = dtx - w*(dtx-dtA);
            enx = particleEnergyForTimestep(count, nm, p, np, epsilon, dtx, nPart, partBoundX,partBoundY, partBoundZ);
        }
        if(enx<eny) {
            dtB = dty;
            enB = eny;
        } else {
            dtA = dtx;
            enA = enx;
            dtx = dty;
            enx = eny;
        }
    } while(dtB-dtA > tol*(dtx+dty));
    //printf("dt: %lf\n", dtx);
    return dtx;
}

void conjGradDescentParticles(int count, double nm[][3], particle p[], int np, double epsilon, int nPart[], double partBoundX[], double partBoundY[], double partBoundZ[]) {
    double energy = totalParticleEnergy(p, np, epsilon);
    double energyOld = energy;
    for(int i=0; i<np; i++) {
        p[i].forceTh = 0;
        p[i].forcePh = 0;
    }
    addForces(count, nm, p, np, epsilon);
    for(int i=0; i<np; i++) {
        p[i].conjGradTh = p[i].forceTh;
        p[i].conjGradPh = p[i].forcePh;
        p[i].forceThOld = p[i].forceTh;
        p[i].forcePhOld = p[i].forcePh;
    }
    
    int conjGradCount = 0;
    double dt = 0;
    
    while(1) {
        dt = goldenSearchParticles(count, nm, p, np, epsilon, nPart, partBoundX, partBoundY, partBoundZ);
        if(dt == 0) {
            if(conjGradCount==0)
                break;
            else {
                conjGradCount = 0;
                for(int i=0; i<np; i++) {
                    p[i].conjGradTh = p[i].forceTh;
                    p[i].conjGradPh = p[i].forcePh;
                    p[i].forceThOld = p[i].forceTh;
                    p[i].forcePhOld = p[i].forcePh;
                }
            }
        }
        else {
            for(int i=0; i<np; i++) {
                p[i].th += p[i].conjGradTh*dt;
                p[i].ph += p[i].conjGradPh*dt;
                p[i].rr = radiusFromCoefficients(count, nm, p[i].th, p[i].ph);
                partitionOneParticle(&(p[i]), nPart, partBoundX, partBoundY, partBoundZ);
            }
            conjGradCount++;
            if(conjGradCount>=1000) {
/*                do{
                int info;
                solveHessian(count, nm, np, p, epsilon, &info);
                if(info==0) {
                    for(int i=0; i<np; i++) {
                        p[i].conjGradTh = p[i].forceTh;
                        p[i].conjGradPh = p[i].forcePh;
                    }
                    dt = goldenSearchParticles(count, nm, p, np, epsilon, nPart, partBoundX, partBoundY, partBoundZ);
                    for(int i=0; i<np; i++) {
                        p[i].th += p[i].conjGradTh*dt;
                        p[i].ph += p[i].conjGradPh*dt;
                        p[i].rr = radiusFromCoefficients(count, nm, p[i].th, p[i].ph);
                        partitionOneParticle(&(p[i]), nPart, partBoundX, partBoundY, partBoundZ);
                    }
                }
                else {
                    printf("Newton method failed!\n");
                    break;
                }
                energy = totalParticleEnergy(p, np, epsilon);
                printf("%d: %.25lf\n", conjGradCount, energy);
                } while(dt!=0);
*/
                conjGradCount = 0;
                energyOld = totalParticleEnergy(p, np, epsilon);
                for(int i=0; i<np; i++) {
                    p[i].forceTh = 0;
                    p[i].forcePh = 0;
                }
                addForces(count, nm, p, np, epsilon);
                for(int i=0; i<np; i++) {
                    p[i].conjGradTh = p[i].forceTh;
                    p[i].conjGradPh = p[i].forcePh;
                    p[i].forceThOld = p[i].forceTh;
                    p[i].forcePhOld = p[i].forcePh;
                }
                continue;
            }
            energy = totalParticleEnergy(p, np, epsilon);
            //printf("%d: %.25lf\n", conjGradCount, energy);
            if(energy==0 || fabs(energy-energyOld)/energyOld<1e-16)
                break;
            energyOld = energy;
            calculateConjugateForces(count, nm, p, np, epsilon);
        }
    }
}

void conjGradDescentParticlesSteps(int count, double nm[][3], particle p[], int np, double epsilon, int nPart[], double partBoundX[], double partBoundY[], double partBoundZ[], int  steps) {
    double energy = totalParticleEnergy(p, np, epsilon);
    double energyOld = energy;
    for(int i=0; i<np; i++) {
        p[i].forceTh = 0;
        p[i].forcePh = 0;
    }
    addForces(count, nm, p, np, epsilon);
    for(int i=0; i<np; i++) {
        p[i].conjGradTh = p[i].forceTh;
        p[i].conjGradPh = p[i].forcePh;
        p[i].forceThOld = p[i].forceTh;
        p[i].forcePhOld = p[i].forcePh;
    }
    
    int conjGradCount = 0;
    double dt = 0;
    int totalCount = 0;
    
    while(1) {
        dt = goldenSearchParticles(count, nm, p, np, epsilon, nPart, partBoundX, partBoundY, partBoundZ);
        if(dt == 0) {
            if(conjGradCount==0)
                break;
            else {
                conjGradCount = 0;
                for(int i=0; i<np; i++) {
                    p[i].conjGradTh = p[i].forceTh;
                    p[i].conjGradPh = p[i].forcePh;
                    p[i].forceThOld = p[i].forceTh;
                    p[i].forcePhOld = p[i].forcePh;
                }
            }
        }
        else {
            for(int i=0; i<np; i++) {
                p[i].th += p[i].conjGradTh*dt;
                p[i].ph += p[i].conjGradPh*dt;
                p[i].rr = radiusFromCoefficients(count, nm, p[i].th, p[i].ph);
                partitionOneParticle(&(p[i]), nPart, partBoundX, partBoundY, partBoundZ);
            }
            conjGradCount++;
            totalCount++;
            energy = totalParticleEnergy(p, np, epsilon);
            //printf("%d: %.25lf\n", conjGradCount, energy);
            if(energy==0 || (conjGradCount==1&&fabs(energy-energyOld)/energyOld<1e-10) || totalCount>steps)
                break;
            if(conjGradCount>=1000 || fabs(energy-energyOld)/energyOld<1e-10) {
                conjGradCount = 0;
                //energyOld = totalParticleEnergy(p, np, epsilon);
                for(int i=0; i<np; i++) {
                    p[i].forceTh = 0;
                    p[i].forcePh = 0;
                }
                addForces(count, nm, p, np, epsilon);
                for(int i=0; i<np; i++) {
                    p[i].conjGradTh = p[i].forceTh;
                    p[i].conjGradPh = p[i].forcePh;
                    p[i].forceThOld = p[i].forceTh;
                    p[i].forcePhOld = p[i].forcePh;
                }
                //continue;
            }
            else {
                calculateConjugateForces(count, nm, p, np, epsilon);
            }
            energyOld = energy;
        }
    }
}


