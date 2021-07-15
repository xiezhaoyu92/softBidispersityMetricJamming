//
//  surface.h
//  softBidispersityMetricJamming
//
//  Created by Zhaoyu Xie on 4/8/20.
//  Copyright © 2020 Zhaoyu Xie. All rights reserved.
//

#ifndef surface_h
#define surface_h

#include <stdio.h>

double tesseralHarmonicsZ(int n, int m, double th, double ph);
double tesseralHarmonicsZDth(int n, int m, double th, double ph);
double tesseralHarmonicsZDph(int n, int m, double th, double ph);
double radiusFromCoefficients(int count, double nm[][3], double th, double ph);
double areaFromCoefficients(int count, double nm[][3]);
double areaFromCoefficientsD(int count, double nm[][3], int c);
//void areaFromCoefficientsDAll(int count, double nm[][3], double areaD[]);
double volumeFromCoefficients(int count, double nm[][3]);
void volumeNormal(int count, double nm[][3], double nv[]);
double volumeNormalNorm(int count, double nm[][3]);
#endif /* surface_h */
