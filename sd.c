//
//  sd.c
//  asd_sphere_20200805
//
//  Created by Mingyang Tan on 8/5/20.
//  Copyright © 2020 Mingyang Tan. All rights reserved.
//

#include "sd.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

struct sd * sd_initialize(void)
{
    struct sd * sys = NULL;
    sys = (struct sd * ) malloc(sizeof(struct sd));
    
    sys->np = 0;
    
    sys->phi = 0.0;
    sys->lx = 0.0;
    sys->ly = 0.0;
    sys->lz = 0.0;
    sys->xlh = 0.0;
    sys->xrh = 0.0;
    sys->ylh = 0.0;
    sys->yrh = 0.0;
    sys->zlh = 0.0;
    sys->zrh = 0.0;
    sys->volume = 0.0;
    sys->pivol = 0.0;
    
    sys->pos = NULL;
    
    sys->lcut = 0.0;
    sys->lcut2 = 0.0;
    
    sys->rcut = 0.0;
    sys->rcut2 = 0.0;
    
    sys->xi = 0.0;
    sys->xi2 = 0.0;
    
    sys->self_a = 0.0;
    sys->self_c = 0.0;
    sys->self_m = 0.0;
    sys->self_d = 0.0;
    
    sys->nkx = 0;
    sys->nky = 0;
    sys->nkz = 0;
    sys->nk = 0;
    sys->dkx = 0.0;
    sys->dky = 0.0;
    sys->dkz = 0.0;
    sys->vk = 0.0;
    
    sys->eta = 0.0;
    sys->p = 0;
    sys->p3 = 0;
    
    sys->gaussCo = 0.0;
    sys->gaussExp = 0.0;
    
    sys->lambda = 0.0;
    sys->itmax = 0;
    sys->tol = 0.0;
    
    sys->deltam = 0.0;
    sys->delta0 = 0.0;
    sys->f0 = 0.0;
    sys->tau = 0.0;
    
    sys->Ma = 0.0;
    sys->direct = 0;
    sys->omega = 0.0;
    
    sys->shear_rate = 0.0;
    
    sys->vginf[0] = 0.0;
    sys->vginf[1] = 0.0;
    sys->vginf[2] = 0.0;
    sys->vginf[3] = 0.0;
    sys->vginf[4] = 0.0;
    sys->vginf[5] = 0.0;
    sys->vginf[6] = 0.0;
    sys->vginf[7] = 0.0;
    sys->vginf[8] = 0.0;
    
    sys->einf[0] = 0.0;
    sys->einf[1] = 0.0;
    sys->einf[2] = 0.0;
    sys->einf[3] = 0.0;
    sys->einf[4] = 0.0;
    
    sys->uinf[0] = 0.0;
    sys->uinf[1] = 0.0;
    sys->uinf[2] = 0.0;
    sys->uinf[3] = 0.0;
    sys->uinf[4] = 0.0;
    sys->uinf[5] = 0.0;
    
    sys->strain = 0.0;
    sys->shift = 0.0;
    
    sys->viscosity = 0.0;
    
    return sys;
};

void sd_np(struct sd * sys,
           int np)
{
    sys->np = np;
}

void sd_cell(struct sd * sys,
             double phi,
             double lx, double ly, double lz)
{
    sys->phi = phi;
    sys->lx = lx;
    sys->ly = ly;
    sys->lz = lz;
    sys->xlh = - lx / 2.0;
    sys->xrh = + lx / 2.0;
    sys->ylh = - ly / 2.0;
    sys->yrh = + ly / 2.0;
    sys->zlh = - lz / 2.0;
    sys->zrh = + lz / 2.0;
    sys->volume = lx * ly * lz;
    sys->pivol = M_PI / sys->volume;
}

void sd_pos(struct sd * sys,
            double * pos)
{
    int n3;
    int i;
    
    n3 = sys->np * 3;
    sys->pos = malloc(n3 * sizeof(double));
    
    for (i = 0; i < n3; i++) {
        sys->pos[i] = pos[i];
    }
    
}

void sd_rcut(struct sd * sys,
             double rcut)
{
    sys->rcut = rcut;
    sys->rcut2 = rcut * rcut;
}

void sd_lcut(struct sd * sys,
             double lcut)
{
    sys->lcut = lcut;
    sys->lcut2 = lcut * lcut;
}

void sd_ewald(struct sd * sys,
              double xi)
{
    double pisqr;
    pisqr = sqrt(M_PI);
    
    sys->xi = xi;
    sys->xi2 = xi * xi;
    
    sys->self_a = 1.0 - xi / pisqr * (6.0 - 40.0 / 3.0 * sys->xi2);
    sys->self_c = 0.75 - xi / pisqr * sys->xi2 * 10.0;
    sys->self_m = 0.9 - xi / pisqr * sys->xi2 * (12.0 - 30.24 * sys->xi2);
    sys->self_d = 1.0 - 4.0 / 3.0 * sys->xi2 / pisqr * sys->xi;
}

void sd_fft(struct sd * sys,
            int nkx, int nky, int nkz)
{
    sys->nkx = nkx;
    sys->nky = nky;
    sys->nkz = nkz;
    sys->nk = nkx * nky * nkz;
    
    sys->dkx = sys->lx / (double) nkx;
    sys->dky = sys->ly / (double) nky;
    sys->dkz = sys->lz / (double) nkz;
    sys->vk = (sys->dkx) * (sys->dky) * (sys->dkz);
}

void sd_gauss(struct sd * sys,
              double eta, int p)
{
    sys->p = p;
    sys->p3 = p * p * p;
    sys->eta = eta;
    
    sys->gaussCo = 2.0 * (sys->xi2) / (M_PI * eta) * sqrt(2.0 * sys->xi2 / (M_PI * eta) );
    sys->gaussExp = 2.0 * sys->xi2 / eta;
}

void sd_repulsive(struct sd * sys,
                  double deltam, double delta0,
                  double f0, double tau)
{
    sys->deltam = deltam;
    sys->delta0 = delta0;
    sys->f0 = f0;
    sys->tau = tau;
}

void sd_magnetic(struct sd * sys,
                 double Ma,
                 int direct,
                 double omega)
{
    sys->Ma = Ma;
    sys->direct = direct;
    sys->omega = omega;
}

void sd_uinf(struct sd * sys,
             double shear_rate,
             double ux, double uy, double uz,
             double ox, double oy, double oz)
{
    sys->shear_rate = shear_rate;
    
    sys->vginf[1] = shear_rate;
    
    sys->einf[0] = 2.0 * sys->vginf[0] + sys->vginf[4];
    sys->einf[1] = sys->vginf[1];
    sys->einf[2] = sys->vginf[2];
    sys->einf[3] = sys->vginf[5];
    sys->einf[4] = 2.0 * sys->vginf[4] + sys->vginf[0];
    
    sys->uinf[0] = ux;
    sys->uinf[1] = uy;
    sys->uinf[2] = uz;
    sys->uinf[3] = 0.5 * (sys->vginf[7] - sys->vginf[5]) + ox;
    sys->uinf[4] = 0.5 * (sys->vginf[2] - sys->vginf[6]) + oy;
    sys->uinf[5] = 0.5 * (sys->vginf[3] - sys->vginf[1]) + oz;
}

void sd_iterative(struct sd * sys,
                  double lambda,
                  int itmax, double tol)
{
    sys->lambda = lambda;
    sys->itmax = itmax;
    sys->tol = tol;
}

void sd_strain(struct sd * sys,
               double t)
{
    int i;
    double shift;
    double strain;
    
    strain = (sys->shear_rate) * t;
    
    i = (int) (strain + 0.5);
    strain = strain - (double) i;
    shift = (sys->ly) * strain;
    
    sys->strain = strain;
    sys->shift = shift;
}

void sd_set_xi(double error,
               double lx, double ly, double lz,
               double * xi,
               double * rcut, int * kcut)
{
    double xi_default;
    double rmax_temp;
    double lmin;
    
    lmin = lx;
    if (lmin > ly) {
        lmin = ly;
    }
    if (lmin > lz) {
        lmin = lz;
    }
    
    xi_default = 0.5;
    rmax_temp = sqrt( - log( error ) ) / xi_default;
    
    if (rmax_temp < lmin / 3.0) {
        * rcut = rmax_temp;
        * xi = xi_default;
    } else {
        * rcut = lmin / 3.0;
        * xi = sqrt( - log( error ) ) / (* rcut);
    }
    
    (* kcut) = (int) ( 2.0 * sqrt( - log( error ) ) * (*xi) ) + 1;
}

void sd_set_fft_point(int kcut,
                      double lx, double ly, double lz,
                      int * nkx, int * nky, int * nkz)
{
    int n;
    int base2;
    
    n = (int) ( (double) kcut * lx / M_PI) + 1;
    base2 = (int) ( log( (double) n ) / log(2.0) + 1.0 );
    * nkx = (int) pow(2.0, (double) base2);
    
    n = (int) ( (double) kcut * ly / M_PI) + 1;
    base2 = (int) ( log( (double) n ) / log(2.0) + 1.0 );
    * nky = (int) pow(2.0, (double) base2);
    
    n = (int) ( (double) kcut * lz / M_PI) + 1;
    base2 = (int) ( log( (double) n ) / log(2.0) + 1.0 );
    * nkz = (int) pow(2.0, (double) base2);
}

void sd_set_eta(struct sd * sys,
                double error,
                double * eta, int * p)
{
    double m;
    double w;
    
    m = 1.0;
    while (erfc( m / sqrt(2.0) ) > error ) {
        m += 0.01;
    }
    (* p) = (int) (m * m / M_PI);
    if ((* p) % 2 == 0) {
        (* p) += 1;
    }
    w = (double) (* p) * (sys->dkx) / 2.0;
    
    (* eta) = 4.0 * w * w * sys->xi2 / (m * m);
    
}
