//
//  lubrication.h
//  asd_sphere_20200508
//
//  Created by Mingyang Tan on 5/9/20.
//  Copyright © 2020 Mingyang Tan. All rights reserved.
//

#ifndef lubrication_h
#define lubrication_h

#include <stdio.h>
#include "sd.h"

struct lcell {
    int tot_boxes;
    int * inter_box;
};

struct lpair {
    int bh_count;
    int * bh_list;
    
    int ah_count;
    int * ah_list;
    
    double * r;
    double * e;
};

struct lub {
    int ndx, ndy, ndz;
    int nd;
    double rx, ry, rz;
    
    int tot_pairs;
    int * row_id;
    int * col_id;
    int * index;
    
    int * head;
    int * list;
    
    double * lfu;
    double * lfe;
    double * lsu;
    double * lse;
    
    int precond;
    double * iccl;
    
    double * lfu_l;
    double * lfu_u;
    
    struct lcell * cell_sys;
    struct lpair * pair_sys;
    
};

struct lub * lub_initialize(struct sd * sys);

void near_lcell_neighbor_shear(struct lub * lub_sys,
                               const int np,
                               const double strain);

void lcell_list(struct lub * lub_sys,
                const int np,
                const double * pos);

void lub_pairs_list(struct lub * lub_sys,
                    const int np,
                    const double xlh, const double xrh,
                    const double ylh, const double yrh,
                    const double zlh, const double zrh,
                    const double lx, const double ly, const double lz,
                    const double lcut2,
                    const double strain,
                    double * pos);

void lub_scalars(double r,
                 double * x11a, double * y11a,
                 double * y11b,
                 double * x11c, double * y11c,
                 double * x11g, double * y11g,
                 double * y11h,
                 double * x12a, double * y12a,
                 double * y12b,
                 double * x12c, double * y12c,
                 double * x12g, double * y12g,
                 double * y12h,
                 double * xm, double * ym, double * zm);

void lub_abc(double x11a, double y11a,
             double y11b,
             double x11c, double y11c,
             double x12a, double y12a,
             double y12b,
             double x12c, double y12c,
             double * e,
             double * labc);

void lub_gh(double x11g, double y11g,
            double y11h,
            double x12g, double y12g,
            double y12h,
            double * e,
            double * lgh);

void lub_zm(double xm, double ym, double zm,
            double * e,
            double * lzm);

void lub_mat(struct lub * lub_sys,
             const int np,
             const double lambda);

void lub_precond(struct lub * lub_sys,
                 const int np, const double lambda);

void lub_precond_2(struct lub * lub_sys,
                   const int np, const double lambda);

#endif /* lubrication_h */
