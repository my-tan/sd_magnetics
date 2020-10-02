//
//  repulsive.c
//  asd_sphere_mag_20200512
//
//  Created by Mingyang Tan on 5/12/20.
//  Copyright © 2020 Mingyang Tan. All rights reserved.
//

#include "repulsive.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "lubrication.h"

void repulsive_hard_potential(struct lub * lub_sys,
                              const int np,
                              const double f0,
                              const double tau,
                              const double deltam, const double delta0,
                              double * f)
{
    int ipair;
    int i, j;
    
    double r;
    double ex, ey, ez;
    double delta;
    double exptau;
    double fac;
    
    for (i = 0; i < np; i++) {
        
        for (ipair = 0; ipair < lub_sys->pair_sys[i].ah_count; ipair++) {
            j = lub_sys->pair_sys[i].ah_list[ipair];
            r = lub_sys->pair_sys[i].r[ipair];
            ex = lub_sys->pair_sys[i].e[ipair * 3    ];
            ey = lub_sys->pair_sys[i].e[ipair * 3 + 1];
            ez = lub_sys->pair_sys[i].e[ipair * 3 + 2];
            
            delta = r - 2.0;
            if (delta > deltam) {
                continue;
            }
            if (delta < delta0) {
                delta = delta0;
            }
            
            exptau = exp(- tau * delta);
            
            fac = f0 * tau * exptau / (1.0 + exptau);
            
            f[i * 3 + 0] -= fac * ex;
            f[i * 3 + 1] -= fac * ey;
            f[i * 3 + 2] -= fac * ez;
            
            f[j * 3 + 0] += fac * ex;
            f[j * 3 + 1] += fac * ey;
            f[j * 3 + 2] += fac * ez;
            
        } // for ipair
        
    } // for i
    
}
