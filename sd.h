//
//  sd.h
//  asd_sphere_20200805
//
//  Created by Mingyang Tan on 8/5/20.
//  Copyright © 2020 Mingyang Tan. All rights reserved.
//

#ifndef sd_h
#define sd_h

#include <stdio.h>

struct sd {
    int np;
    
    double phi;
    double volume;
    double lx, ly, lz;
    double xlh, xrh;
    double ylh, yrh;
    double zlh, zrh;
    double pivol;
    
    double * pos;
    
    double lcut, lcut2;
    double rcut, rcut2;
    
    double xi, xi2;
    double self_a, self_c, self_m, self_d;
    
    int nkx, nky, nkz;
    int nk;
    double dkx, dky, dkz;
    double vk;
    
    double eta;
    int p, p3;
    
    double gaussCo;
    double gaussExp;
    
    double lambda;
    int itmax;
    double tol;
    
    double deltam;
    double delta0;
    double f0;
    double tau;
    
    double Ma;
    int direct;
    double omega;
    
    double shear_rate;
    double vginf[9];
    double einf[5];
    double uinf[6];
    
    double strain;
    double shift;
    
    double viscosity;
};

struct sd * sd_initialize(void);

void sd_np(struct sd * sys,
           int np);

void sd_cell(struct sd * sys,
             double phi,
             double lx, double ly, double lz);

void sd_pos(struct sd * sys,
            double * pos);

void sd_rcut(struct sd * sys,
             double rcut);

void sd_lcut(struct sd * sys,
             double lcut);

void sd_ewald(struct sd * sys,
              double xi);

void sd_fft(struct sd * sys,
            int nkx, int nky, int nkz);

void sd_gauss(struct sd * sys,
              double eta, int p);

void sd_repulsive(struct sd * sys,
                  double deltam, double delta0,
                  double f0, double tau);

void sd_magnetic(struct sd * sys,
                 double Ma,
                 int direct,
                 double omega);

void sd_uinf(struct sd * sys,
             double shear_rate,
             double ux, double uy, double uz,
             double ox, double oy, double oz);

void sd_iterative(struct sd * sys,
                  double lambda,
                  int itmax, double tol);

void sd_strain(struct sd * sys,
               double t);

void sd_set_xi(double error,
               double lx, double ly, double lz,
               double * xi,
               double * rcut, int * kcut);

void sd_set_fft_point(int kcut,
                      double lx, double ly, double lz,
                      int * nkx, int * nky, int * nkz);

void sd_set_eta(struct sd * sys,
                double error,
                double * eta, int * p);

#endif /* sd_h */
