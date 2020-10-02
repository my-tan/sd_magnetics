//
//  repulsive.h
//  asd_sphere_mag_20200512
//
//  Created by Mingyang Tan on 5/12/20.
//  Copyright © 2020 Mingyang Tan. All rights reserved.
//

#ifndef repulsive_h
#define repulsive_h

#include <stdio.h>
#include "lubrication.h"

void repulsive_hard_potential(struct lub * lub_sys,
                              const int np,
                              const double f0,
                              const double tau,
                              const double deltam, const double delta0,
                              double * f);

#endif /* repulsive_h */
