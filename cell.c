//
//  cell.c
//  asd_sphere_20200805
//
//  Created by Mingyang Tan on 8/5/20.
//  Copyright © 2020 Mingyang Tan. All rights reserved.
//

#include "cell.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

/*
 Calculates the fractional position of particles in a unit cell
 
 Inputs:
    lx              dimension of simulation box in x direction
    ly              dimension of simulation box in y direction
    lz              dimension of simulation box in z direction
    xlh             - lx / 2.0
    ylh             - ly / 2.0
    zlh             - lz / 2.0
    strain          deformation of simulation box, suppose in xy-plane
    pos             position of particles
 
 Outputs:
    fpos            fractional position of particles
 */
void pos_fraction(const double lx, const double ly, const double lz,
                  const double xlh, const double ylh, const double zlh,
                  const double strain,
                  const double * pos,
                  double * fpos)
{
    double dx[3];
    dx[0] = pos[0] - xlh;
    dx[1] = pos[1] - ylh;
    dx[2] = pos[2] - zlh;
    
    dx[0] -= strain * pos[1];
    
    fpos[0] = dx[0] / lx;
    fpos[1] = dx[1] / ly;
    fpos[2] = dx[2] / lx;
}

/*
 Return particles back to simulation box
 
 Inputs:
 lx              dimension of simulation box in x direction
 ly              dimension of simulation box in y direction
 lz              dimension of simulation box in z direction
 xlh             - lx / 2.0
 xrh             + lx / 2.0
 ylh             - ly / 2.0
 yrh             + ly / 2.0
 zlh             - lz / 2.0
 zrh             - lz / 2.0
 strain          deformation of simulation box, suppose in xy-plane
 
 Outputs:
 pos             position of particles
 */
void return_to_cell(const double lx, const double ly, const double lz,
                    const double xlh, const double xrh,
                    const double ylh, const double yrh,
                    const double zlh, const double zrh,
                    const double strain,
                    double * pos)
{
    double dx;
    dx = strain * ly;
    
    if (pos[2] < zlh) {
        pos[2] += lz;
    }
    if (pos[2] >= zrh) {
        pos[2] -= lz;
    }
    
    if (pos[1] < ylh) {
        pos[1] += ly;
        pos[0] += dx;
    }
    if (pos[1] >= yrh) {
        pos[1] -= ly;
        pos[0] -= dx;
    }
    
    if (pos[0] < xlh) {
        pos[0] += lx;
    }
    if (pos[0] >= xrh) {
        pos[0] -= lx;
    }
    
}

/*
 Calculates the minimum-image distance between two points
 
 Inputs:
 lx              dimension of simulation box in x direction
 ly              dimension of simulation box in y direction
 lz              dimension of simulation box in z direction
 xlh             - lx / 2.0
 xrh             + lx / 2.0
 ylh             - ly / 2.0
 yrh             + ly / 2.0
 zlh             - lz / 2.0
 zrh             - lz / 2.0
 strain          deformation of simulation box, suppose in xy-plane
 
 Outputs:
 dist            distance between two points
 */
void min_image(const double lx, const double ly, const double lz,
               const double xlh, const double xrh,
               const double ylh, const double yrh,
               const double zlh, const double zrh,
               const double strain,
               double * dist)
{
    if (dist[2] >= lzhigh) {
        dist[2] -= lz;
    }
    if (dist[2] < lzlow) {
        dist[2] += lz;
    }
    
    if (dist[1] >= lyhigh) {
        dist[1] -= ly;
        dist[0] -= ly * strain;
    }
    if (dist[1] < lylow) {
        dist[1] += ly;
        dist[0] += ly * strain;
    }
    
    if (dist[0] >= lxhigh) {
        dist[0] -= lx;
    }
    if (dist[0] < lxlow) {
        dist[0] += lx;
    }
    
}
