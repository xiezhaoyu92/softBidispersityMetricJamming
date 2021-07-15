//
//  energy.h
//  softBidispersityMetricJamming
//
//  Created by Zhaoyu Xie on 4/8/20.
//  Copyright © 2020 Zhaoyu Xie. All rights reserved.
//

#ifndef energy_h
#define energy_h

#include <stdio.h>
#include "particle.h"

double V(particle *p1, particle *p2, double epsilon);
double Vth(int count, double nm[][3], particle *p1, particle *p2, double epsilon);
double Vph(int count, double nm[][3], particle *p1, particle *p2, double epsilon);
double Vc(int count, double nm[][3], int c, particle *p1, particle *p2, double epsilon);
double Ec(int count, double nm[][3], int c, particle p[], int np, double epsilon);
//void EcAll(int count, double nm[][3], double cForce[], particle p[], int np, double epsilon);
double totalEnergy(int count, double nm[][3], particle p[], int np, double epsilon);
double totalParticleEnergy(particle p[], int np, double epsilon);
void coeffForceConstraint(int count, double nm[][3], particle p[], int np, double epsilon, double cForce[]);
void projectOntoConstraint(int count, double nm[][3], double volume);
double volumeNormalProjection(int count, double nm[][3], particle p[], int np, double epsilon);
#endif /* energy_h */
