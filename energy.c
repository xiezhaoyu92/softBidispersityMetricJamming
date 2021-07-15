//
//  energy.c
//  softBidispersityMetricJamming
//
//  Created by Zhaoyu Xie on 4/8/20.
//  Copyright © 2020 Zhaoyu Xie. All rights reserved.
//

#include "energy.h"
#include "surface.h"
#include "particle.h"
#include <math.h>
#include <stdlib.h>

double V(particle *p1, particle *p2, double epsilon) {
    double r1 = p1->rr;
    double r2 = p2->rr;
    double s1 = sin(p1->th);
    double s2 = sin(p2->th);
    double c1 = cos(p1->th);
    double c2 = cos(p2->th);
    double c12 = cos(p1->ph-p2->ph);
    double s12 = sin(p1->ph-p2->ph);
    double R = sqrt(r1*r1+r2*r2-2*r1*r2*c1*c2-2*r1*r2*s1*s2*c12);
    
    double sigma = p1->rad+p2->rad;
    if(R<sigma)
        return epsilon/2*(1-R/sigma)*(1-R/sigma);
    else
        return 0;
}

double Vth(int count, double nm[][3], particle *p1, particle *p2, double epsilon) {
    double r1 = p1->rr;
    double r2 = p2->rr;
    double s1 = sin(p1->th);
    double s2 = sin(p2->th);
    double c1 = cos(p1->th);
    double c2 = cos(p2->th);
    double c12 = cos(p1->ph-p2->ph);
    double s12 = sin(p1->ph-p2->ph);
    double R = sqrt(r1*r1+r2*r2-2*r1*r2*c1*c2-2*r1*r2*s1*s2*c12);
    
    double drth1 = 0;
    for(int i=0; i<count; i++)
        drth1 += nm[i][2]*tesseralHarmonicsZDth((int)(nm[i][0]), (int)(nm[i][1]), p1->th, p1->ph);
    double Rth = (r1*drth1-(drth1*s1+r1*c1)*r2*s2*c12-(drth1*c1-r1*s1)*r2*c2)/R;
    
    double sigma = p1->rad+p2->rad;
    if(R<sigma)
        return -epsilon/sigma*(1-R/sigma)*Rth;
    else
        return 0;
}
double Vph(int count, double nm[][3], particle *p1, particle *p2, double epsilon) {
    double r1 = p1->rr;
    double r2 = p2->rr;
    double s1 = sin(p1->th);
    double s2 = sin(p2->th);
    double c1 = cos(p1->th);
    double c2 = cos(p2->th);
    double c12 = cos(p1->ph-p2->ph);
    double s12 = sin(p1->ph-p2->ph);
    double R = sqrt(r1*r1+r2*r2-2*r1*r2*c1*c2-2*r1*r2*s1*s2*c12);
    
    double drph1 = 0;
    for(int i=0; i<count; i++)
        drph1 += nm[i][2]*tesseralHarmonicsZDph((int)(nm[i][0]), (int)(nm[i][1]), p1->th, p1->ph);
    double Rph = (r1*drph1-(drph1*c12-r1*s12)*r2*s1*s2-drph1*r2*c1*c2)/R;
    
    double sigma = p1->rad+p2->rad;
    if(R<sigma)
        return -epsilon/sigma*(1-R/sigma)*Rph;
    else
        return 0;
}

double Vc(int count, double nm[][3], int c, particle *p1, particle *p2, double epsilon) {
    double Y1 = tesseralHarmonicsZ((int)(nm[c][0]), (int)(nm[c][1]), p1->th, p1->ph);
    double Y2 = tesseralHarmonicsZ((int)(nm[c][0]), (int)(nm[c][1]), p2->th, p2->ph);
    double r1 = p1->rr;
    double r2 = p2->rr;
    double s1 = sin(p1->th);
    double s2 = sin(p2->th);
    double c1 = cos(p1->th);
    double c2 = cos(p2->th);
    double c12 = cos(p1->ph-p2->ph);
    double s12 = sin(p1->ph-p2->ph);
    double R = sqrt(r1*r1+r2*r2-2*r1*r2*c1*c2-2*r1*r2*s1*s2*c12);
    
    double Rc = (r1*Y1+r2*Y2-(Y1*r2+r1*Y2)*(s1*s2*c12+c1*c2))/R;
    
    double sigma = p1->rad+p2->rad;
    if(R<sigma)
        return -epsilon/sigma*(1-R/sigma)*Rc;
    else
        return 0;
}

double Ec(int count, double nm[][3], int c, particle p[], int np, double epsilon) {
    double Vcp = 0;
    for(int i=0; i<np; i++)
        for(int j=i+1; j<np; j++) {
            if(abs(p[i].coord[0]-p[j].coord[0])>1||abs(p[i].coord[1]-p[j].coord[1])>1||abs(p[i].coord[2]-p[j].coord[2])>1)
                continue;
            Vcp += Vc(count, nm, c, &(p[i]), &(p[j]), epsilon);
        }
    return Vcp+areaFromCoefficientsD(count, nm, c);
}
/*
void EcAll(int count, double nm[][3], double cForce[], particle p[], int np, double epsilon) {
    double areaD[count];
    areaFromCoefficientsDAll(count, nm, areaD);
    for(int c=0; c<count; c++){
        double Vcp = 0;
        for(int i=0; i<np; i++)
            for(int j=i+1; j<np; j++) {
                if(abs(p[i].coord[0]-p[j].coord[0])>1||abs(p[i].coord[1]-p[j].coord[1])>1||abs(p[i].coord[2]-p[j].coord[2])>1)
                    continue;
                Vcp += Vc(count, nm, c, &(p[i]), &(p[j]), epsilon);
        }
        cForce[c] = -(Vcp+areaD[c]);
    }
}
*/
double totalEnergy(int count, double nm[][3], particle p[], int np, double epsilon) {
    double Ep = 0;
    for(int i=0; i<np; i++)
        for(int j=i+1; j<np; j++) {
            if(abs(p[i].coord[0]-p[j].coord[0])>1||abs(p[i].coord[1]-p[j].coord[1])>1||abs(p[i].coord[2]-p[j].coord[2])>1)
                continue;
            Ep += V(&(p[i]), &(p[j]), epsilon);
        }
    double Es = areaFromCoefficients(count, nm);
    return Ep+Es;
}

double totalParticleEnergy(particle p[], int np, double epsilon) {
    double Ep = 0;
    for(int i=0; i<np; i++)
        for(int j=i+1; j<np; j++) {
            if(abs(p[i].coord[0]-p[j].coord[0])>1||abs(p[i].coord[1]-p[j].coord[1])>1||abs(p[i].coord[2]-p[j].coord[2])>1)
                continue;
            Ep += V(&(p[i]), &(p[j]), epsilon);
        }
    return Ep;
}

void coeffForceConstraint(int count, double nm[][3], particle p[], int np, double epsilon, double cForce[]){
    for(int c=0; c<count; c++) {
        cForce[c] = -Ec(count, nm, c, p, np, epsilon);
        //printf("%.15lf\n", cForce[c]);
    }
    double nv[count];
    volumeNormal(count, nm, nv);
    double dotProduct = 0;
    for(int c=0; c<count; c++)
        dotProduct += cForce[c]*nv[c];
    for(int c=0; c<count; c++)
        cForce[c] = cForce[c]-dotProduct*nv[c];
}

void projectOntoConstraint(int count, double nm[][3], double volume) {
    double gradNew[count];
    double sigma;
    double lambda;
    double sigmaOld;
    
    sigma = volumeFromCoefficients(count, nm)-volume;
    do {
        sigmaOld = sigma;
        volumeNormal(count, nm, gradNew);
        lambda = sigma/volumeNormalNorm(count, nm);
        for(int c=0; c<count; c++)
            nm[c][2] = nm[c][2]-lambda*gradNew[c];
        sigma = volumeFromCoefficients(count, nm)-volume;
    } while(fabs(sigma)>1e-6&&fabs(sigma)!=fabs(sigmaOld));
}

double volumeNormalProjection(int count, double nm[][3], particle p[], int np, double epsilon) {
    double cForce[count];
    for(int c=0; c<count; c++) {
        cForce[c] = -Ec(count, nm, c, p, np, epsilon);
    }
    double cForceNorm = 0;
    for(int c=0; c<count; c++)
        cForceNorm += cForce[c]*cForce[c];
    cForceNorm = sqrt(cForceNorm);
    double nv[count];
    volumeNormal(count, nm, nv);
    double dotProduct = 0;
    for(int c=0; c<count; c++)
        dotProduct += cForce[c]*nv[c]/cForceNorm;
    return dotProduct;
}


