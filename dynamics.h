//
//  dynamics.h
//  asd_sphere_20200805
//
//  Created by Mingyang Tan on 8/6/20.
//  Copyright © 2020 Mingyang Tan. All rights reserved.
//

#ifndef dynamics_h
#define dynamics_h

#include <stdio.h>
#include "sd.h"
#include "ewald.h"
#include "lubrication.h"

void update_config(struct sd * sys,
                   const double delt,
                   const int ts,
                   const double t,
                   double * d,
                   double * fr, double * fm,
                   double * frx, double * fmx,
                   double * fhn, double * fhf,
                   double * f,
                   double * ulub, double * uoe,
                   double * u, double * s);

void calc_uos(struct sd * sys,
              struct ewald * ewald_sys,
              struct lub * lub_sys,
              double * fr, double * fm,
              double * fhn, double * fhf,
              double * f,
              double * ulub, double * uoe,
              double * u, double * s);

#endif /* dynamics_h */
