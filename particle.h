//
//  particle.h
//  softBidispersityMetricJamming
//
//  Created by Zhaoyu Xie on 4/8/20.
//  Copyright © 2020 Zhaoyu Xie. All rights reserved.
//

#ifndef particle_h
#define particle_h

#include <stdio.h>

#ifndef PI
#define PI 3.141592653589793
#endif

#ifndef MACHINE_EPSILON
#define MACHINE_EPSILON 1e-15
#endif

#ifndef TRUE
#define TRUE -1
#endif

#ifndef FALSE
#define FALSE 0
#endif

typedef struct {
    double rad;
    double r;
    double x;
    double y;
    double z;
    
    double rr;
    double rrLast;
    double rrOld;
    
    double th;
    double ph;
    double thLast;
    double phLast;
    double thOld;
    double phOld;
    
    double forceTh;
    double forcePh;
    double forceThOld;
    double forcePhOld;
    double conjGradTh;
    double conjGradPh;
    
    int oldCoord[3];
    int coord[3];
    int coordLast[3];
} particle;

void coordinateTransform(double x, double y, double z, double a, double b, double *th, double *ph);
void partitionOneParticle(particle *p, int nPart[], double partBoundX[], double partBoundY[], double partBoundZ[]);
void addForces(int count, double nm[][3], particle p[], int np, double epsilon);
void moveSystem(int count, double nm[][3], particle p[], int np, double cForce[], double dt, int nPart[], double partBoundX[], double partBoundY[], double partBoundZ[]);
void deformSurface(int count, double nm[][3], particle p[], int np, double cForce[], double dt, double volume, int nPart[], double partBoundX[], double partBoundY[], double partBoundZ[]);
double goldenSearch(int count, double nm[][3], particle p[], int np, double epsilon, double cForce[], int nPart[], double partBoundX[], double partBoundY[], double partBoundZ[]);

void conjGradDescentParticles(int count, double nm[][3], particle p[], int np, double epsilon, int nPart[], double partBoundX[], double partBoundY[], double partBoundZ[]);
void conjGradDescentParticlesSteps(int count, double nm[][3], particle p[], int np, double epsilon, int nPart[], double partBoundX[], double partBoundY[], double partBoundZ[], int  steps);

#endif /* particle_h */
