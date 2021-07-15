//
//  main.c
//  softBidispersityMetricJamming
//
//  Created by Zhaoyu Xie on 4/8/20.
//  Copyright © 2020 Zhaoyu Xie. All rights reserved.
//

#ifndef TRUE
#define TRUE -1
#endif

#ifndef FALSE
#define FALSE 0
#endif

#ifndef PI
#define PI 3.141592653589793
#endif

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <getopt.h>
#include "particle.h"
#include "energy.h"
#include "surface.h"

int main(int argc, char * argv[]) {
    int debugOutQ = FALSE;
    int debugInQ = FALSE;
    
    double epsilon = 1;
    int readUnfinishedQ = FALSE;
    int animationQ = FALSE;
    
    static struct option longOptions[]={
        {"magnitude",          required_argument,  NULL, 'e'},
        {"readUnfinished",     no_argument,        NULL, 'u'},
        {0,                    0,                  0,     0 }
    };
    int optIndex = 0;
    int opt;
    while ((opt = getopt_long(argc, argv,"e:u",
                              longOptions, &optIndex )) != -1) {
        switch (opt) {
            case 0:
                break;
            case 'e' :
                epsilon = atof(optarg);
                break;
            case 'u' :
                readUnfinishedQ = TRUE;
                break;
            default:
                exit(1);
        }
    }
    
    FILE *eps;
    if(readUnfinishedQ) {
        eps = fopen("epsilon.dat", "r");
        fscanf(eps, "%lf", &epsilon);
    }
    else {
        eps = fopen("epsilon.dat", "w");
        fprintf(eps, "%.15lf\n", epsilon);
    }
    fclose(eps);
    
    double rad;
    FILE *radFile=NULL;
    radFile = fopen("rad.dat","r");
    if (radFile) {
        fscanf(radFile, "%lf", &rad);
        fclose(radFile);
    } else {
        printf("radFile pointer is null\n");
        exit(1);
    }
    
    int np;
    FILE *npFile=NULL;
    npFile = fopen("npts.dat","r");
    if (npFile) {
        fscanf(npFile, "%i", &np);
        fclose(npFile);
    } else {
        printf("npFile pointer is null\n");
        exit(1);
    }
    
    int count = 0;
    FILE *coeffFile=NULL;
    if(readUnfinishedQ)
        coeffFile = fopen("coeffOutTemp.dat","r");
    else
        coeffFile = fopen("coeff.dat","r");
    if (coeffFile) {
        for(char c=getc(coeffFile); c!=EOF; c=getc(coeffFile))
            if(c=='\n')
                count++;
        count++;
        fclose(coeffFile);
    } else {
        printf("coeffFile pointer is null\n");
        exit(1);
    }
    
    double nm[count][3];
    if(readUnfinishedQ)
        coeffFile = fopen("coeffOutTemp.dat","r");
    else
        coeffFile = fopen("coeff.dat","r");
    if (coeffFile) {
        for(int i=0; i<count; i++)
            fscanf(coeffFile, "%lf %lf %lf", nm[i], nm[i]+1, nm[i]+2);
        fclose(coeffFile);
    } else {
        printf("coeffFile pointer is null\n");
        exit(1);
    }
        
    for(int i=0; i<count; i++)
        printf("%d %d %lf\n", (int)(nm[i][0]), (int)(nm[i][1]), nm[i][2]);
    
    particle p[np];
    FILE *configuration=NULL;
    FILE *radiusFile=NULL;
    if(readUnfinishedQ)
        configuration = fopen("configOutTemp.asc","r");
    else
        configuration = fopen("configParam.asc","r");
    radiusFile = fopen("radii.dat", "r");
    if (configuration&&radiusFile) {
        for (int j = 0; j < np; j++) {
            fscanf(configuration, "%lf %lf", &(p[j].th), &(p[j].ph));
            p[j].rr = radiusFromCoefficients(count, nm, p[j].th, p[j].ph);
            p[j].x = p[j].rr*sin(p[j].th)*cos(p[j].ph);
            p[j].y = p[j].rr*sin(p[j].th)*sin(p[j].ph);
            p[j].z = p[j].rr*cos(p[j].th);
            fscanf(radiusFile, "%lf", &(p[j].rad));
        }
        fclose(configuration);
        fclose(radiusFile);
    }
    else {
        printf("configuration or radii pointer is null\n");
        exit(1);
    }

    double L=0;//L represent one side boundary
    for(int i=0; i<np; i++) {
        if(fabs(p[i].x)>L)
            L = fabs(p[i].x);
        if(fabs(p[i].y)>L)
            L = fabs(p[i].y);
        if(fabs(p[i].z)>L)
            L = fabs(p[i].z);
    }
    L = 2*L;//enlarge the box to make sure the deformation is inside.
    
    int nPart[3];
    nPart[0] = (int)(L/(rad))+1;
    nPart[1] = (int)(L/(rad))+1;
    nPart[2] = (int)(L/(rad))+1;
    double partBoundX[nPart[0]];
    double partBoundY[nPart[1]];
    double partBoundZ[nPart[2]];
    for(int i=0; i<nPart[0]; i++)
        partBoundX[i] = -L+i*2*rad;
    for(int i=0; i<nPart[1]; i++)
        partBoundY[i] = -L+i*2*rad;
    for(int i=0; i<nPart[2]; i++)
        partBoundZ[i] = -L+i*2*rad;
    
    for (int j = 0; j < np; j++) {
        partitionOneParticle(&(p[j]), nPart,partBoundX, partBoundY, partBoundZ);
    }
    
    double cForce[count];
    
    conjGradDescentParticles(count, nm, p, np, epsilon, nPart, partBoundX, partBoundY, partBoundZ);
    
    double particleEnergy = totalParticleEnergy(p, np, epsilon);
    double energy = totalEnergy(count, nm, p, np, epsilon);
    double volumeOrig = volumeFromCoefficients(count, nm);
    printf("The original total energy is %.15lf\n", energy);
    printf("The particle interaction energy is %.15lf\n", particleEnergy);
    printf("The volume is %.15lf\n", volumeOrig);
    
    //addForces(count, nm, p, np, epsilon);
    //for(int i=0; i<np; i++)
    //    printf("%d: %.15lf, %.15lf\n", i, p[i].forceTh, p[i].forcePh);
    
    double coeffLast[count];
    for(int i=0; i<count; i++)
        coeffLast[i] = nm[i][2];
    for(int i=0; i<np; i++) {
        p[i].thLast = p[i].th;
        p[i].phLast = p[i].ph;
        p[i].rrLast = p[i].rr;
        for(int j=0; j<3; j++)
            p[i].coordLast[j] = p[i].coord[j];
    }
    double energyOld = energy;
    double t = 0;
    FILE *animation = fopen("animation.dat", "a");
    FILE *energyFile = fopen("energy.dat", "a");
    fprintf(energyFile, "%.15lf %.15lf\n", t, energy);
    if(animationQ){
        for(int i=0; i<count; i++)
            fprintf(animation, "%d %d %.15lf\n", (int)(nm[i][0]), (int)(nm[i][1]), nm[i][2]);
        for(int i=0; i<np; i++)
            fprintf(animation, "%.15lf %.15lf\n", p[i].th, p[i].ph);
    }
    fclose(animation);
    fclose(energyFile);
    
    FILE *configurationOutputTemp = fopen("configOutTemp.asc", "w");
    if(configurationOutputTemp) {
        for(int i=0; i<np; i++)
            fprintf(configurationOutputTemp, "%.15lf %.15lf\n", (p[i].th), (p[i].ph));
        fclose(configurationOutputTemp);
    }
    else {
        printf("configuration output temp file failed to open!\n");
        exit(1);
    }
    FILE *coefficientTempFile = fopen("coeffOutTemp.dat", "w");
    if(coefficientTempFile) {
        for(int i=0; i<count-1; i++)
            fprintf(coefficientTempFile, "%d %d %.15lf\n", (int)(nm[i][0]), (int)(nm[i][1]), nm[i][2]);
        fprintf(coefficientTempFile, "%d %d %.15lf", (int)(nm[count-1][0]), (int)(nm[count-1][1]), nm[count-1][2]);
        fclose(coefficientTempFile);
    }
    else {
        printf("coefficient file temp failed to open!\n");
        exit(1);
    }
    
    int rollBackQ = FALSE;
    double dt = 0;
    double dtOld = 0;
    if(readUnfinishedQ) {
        FILE *dtTempFile = fopen("dtTemp.dat","r");
        fscanf(dtTempFile, "%lf", &dt);
        fclose(dtTempFile);
    }
    while(1) {
        if(!rollBackQ) {
            coeffForceConstraint(count, nm, p, np, epsilon, cForce);

            if(debugOutQ) {
                FILE *cForceFile = fopen("cForceRecord.dat","w");
                for(int c=0; c<count; c++)
                    //printf("%d %d %.15lf\n", (int)(nm[c][0]), (int)(nm[c][1]), cForce[c]);
                    fprintf(cForceFile, "%.15lf\n", cForce[c]);
                fclose(cForceFile);
            }

            FILE *projection = fopen("volumeNormalProjection.dat","a");
            fprintf(projection,"%.15lf %.15lf\n", t, volumeNormalProjection(count, nm, p, np, epsilon));
            fclose(projection);
        
            double dtMax = fabs(nm[0][2]/cForce[0]);
            for(int i=1; i<count; i++){
                if(cForce[i]>0) {
                    if(dtMax>(nm[0][2]-nm[i][2])/cForce[i])
                        dtMax = (nm[0][2]-nm[i][2])/cForce[i];
                    }
                else {
                    if(dtMax>(nm[0][2]+nm[i][2])/fabs(cForce[i]))
                        dtMax = (nm[0][2]+nm[i][2])/fabs(cForce[i]);
                }
            }
            if(dt==0||dt>0.01*dtMax)
                dt = 0.01*dtMax;
        }
        
        if(debugOutQ) {
            FILE *dtDebug = fopen("dtDebug.dat", "w");
            fprintf(dtDebug, "%.15lf\n", dt);
            fclose(dtDebug);
        }
        FILE *dtTempFile = fopen("dtTemp.dat","w");
        fprintf(dtTempFile, "%.15lf\n", dt);
        fclose(dtTempFile);
        deformSurface(count, nm, p, np, cForce, dt, volumeOrig, nPart, partBoundX, partBoundY, partBoundZ);
        if(debugOutQ) {
            FILE* coeffAfterProj = fopen("coeffAfterProj.asc","w");
            for(int i=0; i<count-1; i++)
                fprintf(coeffAfterProj, "%d %d %.15lf\n", (int)(nm[i][0]), (int)(nm[i][1]), nm[i][2]);
            fprintf(coeffAfterProj, "%d %d %.15lf", (int)(nm[count-1][0]), (int)(nm[count-1][1]), nm[count-1][2]);
            fclose(coeffAfterProj);
            FILE* configAfterProj = fopen("configAfterProj.asc","w");
            for(int i=0; i<np; i++)
                fprintf(configAfterProj, "%.15lf %.15lf\n", (p[i].th), (p[i].ph));
            fclose(configAfterProj);
        }
        //conjGradDescentParticles(count, nm, p, np, epsilon, nPart, partBoundX, partBoundY, partBoundZ);
        conjGradDescentParticlesSteps(count, nm, p, np, epsilon, nPart, partBoundX, partBoundY, partBoundZ, 1e6);
        energy = totalEnergy(count, nm, p, np, epsilon);
        double volume = volumeFromCoefficients(count, nm);
        for(int i=0; i<count; i++)
            printf("%d %d %.15lf\n", (int)(nm[i][0]), (int)(nm[i][1]), nm[i][2]);
        printf("energy: %.20lf\n", energy);
        if(debugOutQ) {
            FILE* configAfterRelax = fopen("configAfterRelax.asc","w");
            for(int i=0; i<np; i++)
                fprintf(configAfterRelax, "%.15lf %.15lf\n", (p[i].th), (p[i].ph));
            fclose(configAfterRelax);
            FILE *energyRecord = fopen("energyRecord.dat", "a");
            fprintf(energyRecord, "%.15lf %.15lf\n", dt, energy);
            fclose(energyRecord);
        }
        if(energy<energyOld) {
            if(rollBackQ&&fabs(energy-energyOld)/energyOld<1e-10)
                break;

            for(int i=0; i<count; i++)
                coeffLast[i] = nm[i][2];
            for(int i=0; i<np; i++) {
                p[i].thLast = p[i].th;
                p[i].phLast = p[i].ph;
                p[i].rrLast = p[i].rr;
                for(int j=0; j<3; j++)
                    p[i].coordLast[j] = p[i].coord[j];
            }
            t = t+dt;
            energyOld = energy;
            if(rollBackQ) {
                rollBackQ = FALSE;
                dtOld = dt;
            }
            else {
                if(dt==dtOld)
                    dt = dt*2;
                else
                    dtOld = dt;
            }
            animation = fopen("animation.dat", "a");
            energyFile = fopen("energy.dat", "a");
            fprintf(energyFile, "%.15lf %.15lf\n", t, energy);
            if(animationQ){
                for(int i=0; i<count; i++)
                    fprintf(animation, "%d %d %.15lf\n", (int)(nm[i][0]), (int)(nm[i][1]), nm[i][2]);
                for(int i=0; i<np; i++)
                    fprintf(animation, "%.15lf %.15lf\n", p[i].th, p[i].ph);
            }
            fclose(animation);
            fclose(energyFile);
            
            configurationOutputTemp = fopen("configOutTemp.asc", "w");
            if(configurationOutputTemp) {
                for(int i=0; i<np; i++)
                    fprintf(configurationOutputTemp, "%.15lf %.15lf\n", (p[i].th), (p[i].ph));
                fclose(configurationOutputTemp);
            }
            else {
                printf("configuration output temp file failed to open!\n");
                exit(1);
            }
            coefficientTempFile = fopen("coeffOutTemp.dat", "w");
            if(coefficientTempFile) {
                for(int i=0; i<count-1; i++)
                    fprintf(coefficientTempFile, "%d %d %.15lf\n", (int)(nm[i][0]), (int)(nm[i][1]), nm[i][2]);
                fprintf(coefficientTempFile, "%d %d %.15lf", (int)(nm[count-1][0]), (int)(nm[count-1][1]), nm[count-1][2]);
                fclose(coefficientTempFile);
            }
            else {
                printf("coefficient file temp failed to open!\n");
                exit(1);
            }
        }
        else {
            rollBackQ = TRUE;
            dt = dt/2;
            for(int i=0; i<count; i++)
                nm[i][2] = coeffLast[i];
            for(int i=0; i<np; i++) {
                p[i].th = p[i].thLast;
                p[i].ph = p[i].phLast;
                p[i].rr = p[i].rrLast;
                for(int j=0; j<3; j++)
                    p[i].coord[j] = p[i].coordLast[j];
            }
            if(dt < MACHINE_EPSILON)
                break;
        }
    }
    
    FILE *forceSurface = fopen("forcesSurface.txt", "w");
    for(int i=0; i<count; i++) {
        fprintf(forceSurface, "%.15lf\n", cForce[i]);
    }
    fclose(forceSurface);
    FILE *forceParticle = fopen("forcesParticles.txt", "w");
    for(int i=0; i<np; i++)
        fprintf(forceParticle, "%.15lf %.15lf\n", p[i].forceTh, p[i].forcePh);
    fclose(forceParticle);
    
    
    FILE *configurationOutput = NULL;
    configurationOutput = fopen("configOut.asc", "w");
    if(configurationOutput) {
        for(int i=0; i<np; i++)
            fprintf(configurationOutput, "%.15lf %.15lf\n", (p[i].th), (p[i].ph));
        fclose(configurationOutput);
    }
    else {
        printf("configuration output file failed to open!\n");
        exit(1);
    }
    FILE *coefficientFile = NULL;
    coefficientFile = fopen("coeffOut.dat", "w");
    if(coefficientFile) {
        for(int i=0; i<count-1; i++)
            fprintf(coefficientFile, "%d %d %.15lf\n", (int)(nm[i][0]), (int)(nm[i][1]), nm[i][2]);
        fprintf(coefficientFile, "%d %d %.15lf", (int)(nm[count-1][0]), (int)(nm[count-1][1]), nm[count-1][2]);
        fclose(coefficientFile);
    }
    else {
        printf("coefficient file failed to open!\n");
        exit(1);
    }

    printf("Hello, World!\n");
    return 0;
}
