//
//  magnetics.c
//  asd_sphere_mag_20200512
//
//  Created by Mingyang Tan on 5/12/20.
//  Copyright © 2020 Mingyang Tan. All rights reserved.
//

#include "magnetics.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "sd.h"
#include "ewald.h"
#include "iterative.h"

/*
 Calculate magnetic force
 */
void magnetics(struct sd * sys,
               struct ewald * ewald_sys,
               const double t,
               double * d,
               double * fm)
{
    int n3;
    int i, j;
    
    double b[3];
    
    n3 = sys->np * 3;
    
    /* Obtain magnetic field */
    if (sys->omega == 0.0)         {    // static field in x-direction
        if (sys->direct == 0) {
            b[0] = 1.0;
            b[1] = 0.0;
            b[2] = 0.0;
        }
        else if (sys->direct == 1) {    // static field in y-direction
            b[0] = 0.0;
            b[1] = 1.0;
            b[2] = 0.0;
        }
        else if (sys->direct == 2) {    // static field in z-direction
            b[0] = 0.0;
            b[1] = 0.0;
            b[2] = 1.0;
        }
        
    } else {
        if (sys->direct == 0)      {    // rotating field in xy-plane
            b[0] = cos(sys->omega * t);
            b[1] = sin(sys->omega * t);
            b[2] = 0.0;
        }
        else if (sys->direct == 1) {    // rotating field in xz-plane
            b[0] = cos(sys->omega * t);
            b[1] = 0.0;
            b[2] = sin(sys->omega * t);
        }
        else if (sys->direct == 2) {    // rotating field in yz-plane
            b[0] = 0.0;
            b[1] = cos(sys->omega * t);
            b[2] = sin(sys->omega * t);
        }
    }
    
    /*
     Obtain magnetic dipole moments of particles
     d = M^{-1} . B
     using GMRES
     */
    pgmres_mag(sys,
               ewald_sys,
               b,
               d);
    
    /*
     Calculate magnetic force
     */
    force_wrap(ewald_sys,
               sys->np,
               sys->xi, sys->xi2,
               sys->p, sys->p3,
               sys->nkx, sys->nky, sys->nkz,
               sys->nk,
               sys->gaussCo, sys->gaussExp,
               sys->vk, sys->eta,
               d,
               fm);
    
    for (i = 0; i < n3; i++) {
        fm[i] /= sys->Ma;
    }
    
}
