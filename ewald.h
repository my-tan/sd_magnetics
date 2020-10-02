//
//  ewald.h
//  asd_sphere_20200805
//
//  Created by Mingyang Tan on 8/5/20.
//  Copyright © 2020 Mingyang Tan. All rights reserved.
//

#ifndef ewald_h
#define ewald_h

#include <stdio.h>
#include <fftw3.h>
#include "sd.h"

struct mcell {
    int tot_boxes;
    
    int * inter_box;
};

struct gridK {
    double kx;
    double ky;
    double kz;
    double k;
    
    double ya;
    double yb;
    double yc;
    double yg;
    double yh;
    double ym;
    
    double w;
};

struct gridP {
    int * ids;
    double * r2;
    double * dx;
    double * dy;
    double * dz;
};

struct mpair {
    int bh_count;
    int * bh_list;
    
    int ah_count;
    int * ah_list;
    
    double * r;
    double * e;
};

struct ewald {
    int ndx, ndy, ndz;
    int nd;
    double rx, ry, rz;
    
    int tot_pairs;
    int * row_id;
    int * col_id;
    int * index;
    
    int * head;
    int * list;
    
    double * muf;
    double * mus;
    double * mef;
    double * mes;
    
    double * mbd;
    
    struct mcell * cell_sys;
    struct gridP * gridp_sys;
    struct gridK * gridk_sys;
    struct mpair * pair_sys;
};

struct ewald * ewald_initialize(struct sd * sys);

void gridk_value(struct ewald * ewald_sys,
                 const double lx, const double ly, const double lz,
                 const int nkx, const int nky, const int nkz,
                 const int nk,
                 const double xi2, const double eta,
                 const double strain);

void id_part_grid(struct ewald * ewald_sys,
                  const int np,
                  const int p, const int p3,
                  const double xlh, const double xrh,
                  const double ylh, const double yrh,
                  const double zlh, const double zrh,
                  const double lx, const double ly, const double lz,
                  const int nkx, const int nky, const int nkz,
                  const double dkx, const double dky, const double dkz,
                  const double strain,
                  const double * pos);

void near_rcell_neighbor_shear(struct ewald * ewald_sys,
                               const int np,
                               const double strain);

void rcell_list(struct ewald * ewald_sys,
                const int np,
                const double * pos);

void mob_pairs_list(struct ewald * ewald_sys,
                    const int np,
                    const double lxlow, const double lxhigh,
                    const double lylow, const double lyhigh,
                    const double lzlow, const double lzhigh,
                    const double lx, const double ly, const double lz,
                    const double rcut2,
                    const double strain,
                    double * pos);

void ewald_real_scalars(double xi, double r,
                        double * xa, double * ya,
                        double * yb,
                        double * xc, double * yc,
                        double * xg, double * yg,
                        double * yh,
                        double * xm, double * ym, double * zm);

void mob_abc(double xa, double ya,
             double yb,
             double xc, double yc,
             double * e,
             double * mabc);

void mob_gh(double xg, double yg,
            double yh,
            double * e,
            double * mgh);

void mob_zm(double xm, double ym, double zm,
            double * e,
            double * mzm);

void ewald_real_mat(struct ewald * ewald_sys,
                    const int np,
                    const double xi);

void ewald_recip_spread(struct ewald * ewald_sys,
                        const int np, const int p3,
                        const double gaussCo, const double gaussExp,
                        double * f, double * s,
                        fftw_complex * gridFx, fftw_complex * gridFy, fftw_complex * gridFz,
                        fftw_complex * gridTx, fftw_complex * gridTy, fftw_complex * gridTz,
                        fftw_complex * gridSxx, fftw_complex * gridSxy, fftw_complex * gridSxz, fftw_complex * gridSyz, fftw_complex * gridSyy);

void ewald_recip_green(struct ewald * ewald_sys,
                       const int nk,
                       fftw_complex * gridFx, fftw_complex * gridFy, fftw_complex * gridFz,
                       fftw_complex * gridTx, fftw_complex * gridTy, fftw_complex * gridTz,
                       fftw_complex * gridSxx, fftw_complex * gridSxy, fftw_complex * gridSxz, fftw_complex * gridSyz, fftw_complex * gridSyy);

void ewald_recip_interp(struct ewald * ewald_sys,
                        const int np, const int p3,
                        const double gaussCo, const double gaussExp,
                        const double vk,
                        fftw_complex * gridFx, fftw_complex * gridFy, fftw_complex * gridFz,
                        fftw_complex * gridTx, fftw_complex * gridTy, fftw_complex * gridTz,
                        fftw_complex * gridSxx, fftw_complex * gridSxy, fftw_complex * gridSxz, fftw_complex * gridSyz, fftw_complex * gridSyy,
                        double * u, double * e);

void ewald_recip_wrap(struct ewald * ewald_sys,
                      const int np, const int p3,
                      const int nkx, const int nky, const int nkz,
                      const int nk,
                      const double gaussCo, const double gaussExp,
                      const double vk,
                      double * f, double * s,
                      double * u, double * e);

void ewald_mob_prob_wrap(struct ewald * ewald_sys,
                         const int np, const int p3,
                         const int nkx, const int nky, const int nkz,
                         const int nk,
                         const double gaussCo, const double gaussExp,
                         const double vk,
                         const double self_a, const double self_c, const double self_m,
                         double * f, double * s,
                         double * u, double * e);

void dipole_real_scalars(double xi ,double r,
                         double * xd, double * yd);

void dip_a(double xd, double yd,
           double * e,
           double * mb);

void dipole_real_mat(struct ewald * ewald_sys,
                     const int np,
                     const double xi);

void dipole_recip_spread(struct ewald * ewald_sys,
                         const int np, const int p3,
                         const double gaussCo, const double gaussExp,
                         double * d,
                         fftw_complex * gridDx, fftw_complex * gridDy, fftw_complex * gridDz);

void dipole_recip_green(struct ewald * ewald_sys,
                        const int nk,
                        fftw_complex * gridDx, fftw_complex * gridDy, fftw_complex * gridDz);

void field_recip_interp(struct ewald * ewald_sys,
                        const int np, const int p3,
                        const double gaussCo, const double gaussExp,
                        const double vk,
                        fftw_complex * gridDx, fftw_complex * gridDy, fftw_complex * gridDz,
                        double * b);

void force_recip_interp(struct ewald * ewald_sys,
                        const int np, const int p3,
                        const double xi2, const double eta,
                        const double gaussCo, const double gaussExp,
                        const double vk,
                        fftw_complex * gridDx, fftw_complex * gridDy, fftw_complex * gridDz,
                        const double * d,
                        double * f);

void field_recip_wrap(struct ewald * ewald_sys,
                      const int np, const int p3,
                      const int nkx, const int nky, const int nkz,
                      const int nk,
                      const double gaussCo, const double gaussExp,
                      const double vk,
                      double * d,
                      double * b);

void field_wrap(struct ewald * ewald_sys,
                const int np, const int p3,
                const int nkx, const int nky, const int nkz,
                const int nk,
                const double gaussCo, const double gaussExp,
                const double vk,
                const double self_d,
                double * d,
                double * b);

void force_recip_wrap(struct ewald * ewald_sys,
                      const int np, const int p3,
                      const int nkx, const int nky, const int nkz,
                      const int nk,
                      const double xi2, const double eta,
                      const double gaussCo, const double gaussExp,
                      const double vk,
                      double * d,
                      double * f);

void force_wrap(struct ewald * ewald_sys,
                const int np,
                const double xi, const double xi2,
                const int p, const int p3,
                const int nkx, const int nky, const int nkz,
                const int nk,
                const double gaussCo, const double gaussExp,
                const double vk, const double eta,
                double * d,
                double * f);

#endif /* ewald_h */
