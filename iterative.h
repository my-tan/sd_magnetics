//
//  iterative.h
//  asd_sphere_20200508
//
//  Created by Mingyang Tan on 5/11/20.
//  Copyright © 2020 Mingyang Tan. All rights reserved.
//

#ifndef iterative_h
#define iterative_h

#include <stdio.h>
#include "sd.h"
#include "lubrication.h"
#include "ewald.h"

void cg(struct lub * lub_sys,
        const int np,
        const int itmax, const double tol,
        const double * b,
        double * x);

void pcg_iccf(struct lub * lub_sys,
              const int np,
              const int itmax, const double tol,
              const double * b,
              double * x);

void pgmres_mob(struct sd * sys,
                struct ewald * ewald_sys,
                struct lub * lub_sys,
                const double * b,
                double * x);

void pgmres_mag(struct sd * sys,
                struct ewald * ewald_sys,
                const double * b,
                double * x);

#endif /* iterative_h */
