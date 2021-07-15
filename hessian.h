//
//  hessian.h
//  softBidispersityMetricJamming
//
//  Created by Zhaoyu Xie on 6/6/20.
//  Copyright © 2020 Zhaoyu Xie. All rights reserved.
//

#ifndef hessian_h
#define hessian_h

#include <stdio.h>
#include "particle.h"

void solveHessian(int count, double nm[][3], int np, particle p[], double epsilon, int *info);

#endif /* hessian_h */
