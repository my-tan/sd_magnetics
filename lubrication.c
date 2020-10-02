//
//  lubrication.c
//  asd_sphere_20200508
//
//  Created by Mingyang Tan on 5/9/20.
//  Copyright © 2020 Mingyang Tan. All rights reserved.
//

#include "lubrication.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "sd.h"
#include "mynla.h"
#include "cell.h"

extern double *rabc;
extern double *x11as;
extern double *x12as;
extern double *y11as;
extern double *y12as;
extern double *y11bs;
extern double *y12bs;
extern double *x11cs;
extern double *x12cs;
extern double *y11cs;
extern double *y12cs;

extern double *rghm;
extern double *x11gs;
extern double *x12gs;
extern double *y11gs;
extern double *y12gs;
extern double *y11hs;
extern double *y12hs;
extern double *xms;
extern double *yms;
extern double *zms;

extern double * kron;
extern double * levi;

struct lub * lub_initialize(struct sd * sys)
{
    int i;
    
    struct lub * lub_sys = NULL;
    lub_sys = (struct lub * ) malloc(sizeof(struct lub));
    
    lub_sys->ndx = floor(sys->lx / sys->lcut);
    lub_sys->ndy = floor(sys->ly / sys->lcut);
    lub_sys->ndz = floor(sys->lz / sys->lcut);
    lub_sys->nd = lub_sys->ndx * lub_sys->ndy * lub_sys->ndz;
    lub_sys->rx = sys->lx / (double) lub_sys->ndx;
    lub_sys->ry = sys->ly / (double) lub_sys->ndy;
    lub_sys->rz = sys->lz / (double) lub_sys->ndz;
    
    lub_sys->tot_pairs = 0;
    lub_sys->row_id = NULL;
    lub_sys->col_id = NULL;
    lub_sys->index = NULL;
    lub_sys->head = NULL;
    lub_sys->list = NULL;
    
    lub_sys->lfu = NULL;
    lub_sys->lfe = NULL;
    lub_sys->lsu = NULL;
    lub_sys->lse = NULL;
    
    lub_sys->iccl = NULL;
    
    lub_sys->lfu_l = NULL;
    lub_sys->lfu_u = NULL;
    
    lub_sys->cell_sys = (struct lcell * ) malloc(lub_sys->nd * sizeof(struct lcell));
    for (i = 0; i < lub_sys->nd; i++) {
        lub_sys->cell_sys[i].tot_boxes = 0;
        lub_sys->cell_sys[i].inter_box = NULL;
    }
    
    lub_sys->pair_sys = (struct lpair * ) malloc(sys->np * sizeof(struct lpair));
    for (i = 0; i < sys->np; i++) {
        lub_sys->pair_sys[i].bh_count = 0;
        lub_sys->pair_sys[i].bh_list = NULL;
        
        lub_sys->pair_sys[i].ah_count = 0;
        lub_sys->pair_sys[i].ah_list = NULL;
        
        lub_sys->pair_sys[i].r = NULL;
        lub_sys->pair_sys[i].e = NULL;
    }
    
    return lub_sys;
};

static int cell_id(int ix, int iy, int iz,
                   int ndx, int ndy, int ndz)
{
    int ndx3, ndy3, ndz3;
    int id;
    
    ndx3 = ndx * 3;
    ndy3 = ndy * 3;
    ndz3 = ndz * 3;
    
    id = (ix + ndx3) % ndx + ((iy + ndy3) % ndy) * ndx + ((iz + ndz3) % ndz) * (ndx * ndy);
    
    return id;
}

void near_lcell_neighbor_shear(struct lub * lub_sys,
                               const int np,
                               const double strain)
{
    int ndx, ndy, ndz;
    int px;
    int ix, iy, iz;
    int id;
    
    double sx;
    
    ndx = lub_sys->ndx;
    ndy = lub_sys->ndy;
    ndz = lub_sys->ndz;
    
    for (iz = 0; iz < ndz; iz++) {
        
        for (iy = 0; iy < ndy - 1; iy++) {
            
            for (ix = 0; ix < ndx; ix++) {
                id = cell_id(ix, iy, iz,
                             ndx, ndy, ndz);
                
                lub_sys->cell_sys[id].tot_boxes = 16;
                lub_sys->cell_sys[id].inter_box = malloc(16 * sizeof(int));
                
                lub_sys->cell_sys[id].inter_box[0] = cell_id(ix + 1, iy    , iz,
                                                             ndx, ndy, ndz);
                lub_sys->cell_sys[id].inter_box[1] = cell_id(ix + 1, iy + 1, iz,
                                                             ndx, ndy, ndz);
                lub_sys->cell_sys[id].inter_box[2] = cell_id(ix    , iy + 1, iz,
                                                             ndx, ndy, ndz);
                lub_sys->cell_sys[id].inter_box[3] = cell_id(ix - 1, iy + 1, iz,
                                                             ndx, ndy, ndz);
                lub_sys->cell_sys[id].inter_box[4] = cell_id(ix + 1, iy    , iz - 1,
                                                             ndx, ndy, ndz);
                lub_sys->cell_sys[id].inter_box[5] = cell_id(ix + 1, iy + 1, iz - 1,
                                                             ndx, ndy, ndz);
                lub_sys->cell_sys[id].inter_box[6] = cell_id(ix    , iy + 1, iz - 1,
                                                             ndx, ndy, ndz);
                lub_sys->cell_sys[id].inter_box[7] = cell_id(ix - 1, iy + 1, iz - 1,
                                                             ndx, ndy, ndz);
                lub_sys->cell_sys[id].inter_box[8] = cell_id(ix + 1, iy    , iz + 1,
                                                             ndx, ndy, ndz);
                lub_sys->cell_sys[id].inter_box[9] = cell_id(ix + 1, iy + 1, iz + 1,
                                                             ndx, ndy, ndz);
                lub_sys->cell_sys[id].inter_box[10] = cell_id(ix    , iy + 1, iz + 1,
                                                              ndx, ndy, ndz);
                lub_sys->cell_sys[id].inter_box[11] = cell_id(ix - 1, iy + 1, iz + 1,
                                                              ndx, ndy, ndz);
                lub_sys->cell_sys[id].inter_box[12] = cell_id(ix    , iy    , iz + 1,
                                                              ndx, ndy, ndz);
                lub_sys->cell_sys[id].inter_box[13] = - 1;
                lub_sys->cell_sys[id].inter_box[14] = - 1;
                lub_sys->cell_sys[id].inter_box[15] = - 1;
                
            } // for ix
            
        } // for iy
        
    } // for iz
    
    // Top layers
    iy = ndy - 1;
    
    sx = strain - round(strain);
    px = (int) ((1.0 + sx) * (double) ndx);
    
    for (iz = 0; iz < ndz; iz++) {
        
        for (ix = 0; ix < ndx; ix++) {
            
            id = cell_id(ix, iy, iz,
                         ndx, ndy, ndz);
            
            lub_sys->cell_sys[id].tot_boxes = 16;
            lub_sys->cell_sys[id].inter_box = malloc(16 * sizeof(int));
            
            lub_sys->cell_sys[id].inter_box[0] = cell_id(ix + 1,      iy    , iz,
                                                         ndx, ndy, ndz);
            lub_sys->cell_sys[id].inter_box[1] = cell_id(ix + 1 - px, iy + 1, iz,
                                                         ndx, ndy, ndz);
            lub_sys->cell_sys[id].inter_box[2] = cell_id(ix     - px, iy + 1, iz,
                                                         ndx, ndy, ndz);
            lub_sys->cell_sys[id].inter_box[3] = cell_id(ix - 1 - px, iy + 1, iz,
                                                         ndx, ndy, ndz);
            lub_sys->cell_sys[id].inter_box[4] = cell_id(ix + 1,      iy    , iz - 1,
                                                         ndx, ndy, ndz);
            lub_sys->cell_sys[id].inter_box[5] = cell_id(ix + 1 - px, iy + 1, iz - 1,
                                                         ndx, ndy, ndz);
            lub_sys->cell_sys[id].inter_box[6] = cell_id(ix     - px, iy + 1, iz - 1,
                                                         ndx, ndy, ndz);
            lub_sys->cell_sys[id].inter_box[7] = cell_id(ix - 1 - px, iy + 1, iz - 1,
                                                         ndx, ndy, ndz);
            lub_sys->cell_sys[id].inter_box[8] = cell_id(ix + 1,      iy    , iz + 1,
                                                         ndx, ndy, ndz);
            lub_sys->cell_sys[id].inter_box[9] = cell_id(ix + 1 - px, iy + 1, iz + 1,
                                                         ndx, ndy, ndz);
            lub_sys->cell_sys[id].inter_box[10] = cell_id(ix     - px, iy + 1, iz + 1,
                                                          ndx, ndy, ndz);
            lub_sys->cell_sys[id].inter_box[11] = cell_id(ix - 1 - px, iy + 1, iz + 1,
                                                          ndx, ndy, ndz);
            lub_sys->cell_sys[id].inter_box[12] = cell_id(ix         , iy    , iz + 1,
                                                          ndx, ndy, ndz);
            
            lub_sys->cell_sys[id].inter_box[13] = cell_id(ix - 2 - px, iy + 1, iz,
                                                          ndx, ndy, ndz);
            lub_sys->cell_sys[id].inter_box[14] = cell_id(ix - 2 - px, iy + 1, iz - 1,
                                                          ndx, ndy, ndz);
            lub_sys->cell_sys[id].inter_box[15] = cell_id(ix - 2 - px, iy + 1, iz + 1,
                                                          ndx, ndy, ndz);
            
            if (lub_sys->cell_sys[id].inter_box[13] == lub_sys->cell_sys[id].inter_box[1]) {
                lub_sys->cell_sys[id].inter_box[13] = -1;
            }
            
            if (lub_sys->cell_sys[id].inter_box[14] == lub_sys->cell_sys[id].inter_box[5]) {
                lub_sys->cell_sys[id].inter_box[14] = -1;
            }
            
            if (lub_sys->cell_sys[id].inter_box[15] == lub_sys->cell_sys[id].inter_box[9]) {
                lub_sys->cell_sys[id].inter_box[15] = -1;
            }
        } // for ix
        
    } // for iz
    
}

void lcell_list(struct lub * lub_sys,
                const int np,
                const double * pos)
{
    int ndx, ndy, ndz;
    int ix, iy, iz;
    int id;
    int i;
    
    ndx = lub_sys->ndx;
    ndy = lub_sys->ndy;
    ndz = lub_sys->ndz;
    
    lub_sys->head = calloc(lub_sys->nd, sizeof(int));
    lub_sys->list = calloc(np, sizeof(int));
    
    for (i = 0; i < lub_sys->nd; i++) {
        lub_sys->head[i] = -1;
    }
    
    for (i = 0; i < np; i++) {
        ix = (int) (pos[i * 3 + 0] / lub_sys->rx);
        iy = (int) (pos[i * 3 + 1] / lub_sys->ry);
        iz = (int) (pos[i * 3 + 2] / lub_sys->rz);
        
        id = cell_id(ix, iy, iz,
                     ndx, ndy, ndz);
        lub_sys->list[i] = lub_sys->head[id];
        lub_sys->head[id] = i;
    }
    
}

/*
 Merge sort array a[]
 */
static void merge_sort_pair(const int lb, const int ub,
                            int a[],
                            int aux[])
{
    int mid;
    int i;
    int pl, pr;
    
    if (ub <= lb) {
        return;
    }
    
    mid = (int) (lb + ub) / 2;
    
    merge_sort_pair(lb, mid,
                    a,
                    aux);
    merge_sort_pair(mid + 1, ub,
                    a,
                    aux);
    
    pl = lb;
    pr = mid + 1;
    
    for (i = lb; i <= ub; i++) {
        if (pl == mid + 1) {
            aux[i] = a[pr];
            pr++;
        }
        else if (pr == ub + 1) {
            aux[i] = a[pl];
            pl++;
        }
        else if (a[pl] < a[pr]) {
            aux[i] = a[pl];
            pl++;
        }
        else {
            aux[i] = a[pr];
            pr++;
        } // if ... else ...
        
    } // for i
    
    for (i = lb; i <= ub; i++) {
        a[i] = aux[i];
    }
    
}

/*
 Merge sort array a[] and sort r[] and e[] corresponding to the order of a[]
 */
static void merge_sort_pair_all(const int lb, const int ub,
                                int a[], double r[], double e[],
                                int aux[], double aux_r[], double aux_e[])
{
    int mid;
    int i;
    int pl, pr;
    
    if (ub <= lb) {
        return;
    }
    
    mid = (int) (lb + ub) / 2;
    
    merge_sort_pair_all(lb, mid,
                        a, r, e,
                        aux, aux_r, aux_e);
    merge_sort_pair_all(mid + 1, ub,
                        a, r, e,
                        aux, aux_r, aux_e);
    
    pl = lb;
    pr = mid + 1;
    
    for (i = lb; i <= ub; i++) {
        
        if (pl == mid + 1) {
            aux[i] = a[pr];
            aux_r[i] = r[pr];
            aux_e[i * 3 + 0] = e[pr * 3 + 0];
            aux_e[i * 3 + 1] = e[pr * 3 + 1];
            aux_e[i * 3 + 2] = e[pr * 3 + 2];
            pr++;
        }
        else if (pr == ub + 1) {
            aux[i] = a[pl];
            aux_r[i] = r[pl];
            aux_e[i * 3 + 0] = e[pl * 3 + 0];
            aux_e[i * 3 + 1] = e[pl * 3 + 1];
            aux_e[i * 3 + 2] = e[pl * 3 + 2];
            pl++;
        }
        else if (a[pl] < a[pr]) {
            aux[i] = a[pl];
            aux_r[i] = r[pl];
            aux_e[i * 3 + 0] = e[pl * 3 + 0];
            aux_e[i * 3 + 1] = e[pl * 3 + 1];
            aux_e[i * 3 + 2] = e[pl * 3 + 2];
            pl++;
        }
        else {
            aux[i] = a[pr];
            aux_r[i] = r[pr];
            aux_e[i * 3 + 0] = e[pr * 3 + 0];
            aux_e[i * 3 + 1] = e[pr * 3 + 1];
            aux_e[i * 3 + 2] = e[pr * 3 + 2];
            pr++;
        } // if ... else ...
        
    } // for i
    
    for (i = lb; i <= ub; i++) {
        a[i] = aux[i];
        r[i] = aux_r[i];
        e[i * 3 + 0] = aux_e[i * 3 + 0];
        e[i * 3 + 1] = aux_e[i * 3 + 1];
        e[i * 3 + 2] = aux_e[i * 3 + 2];
    } // for i
    
}

void lub_pairs_list(struct lub * lub_sys,
                    const int np,
                    const double xlh, const double xrh,
                    const double ylh, const double yrh,
                    const double zlh, const double zrh,
                    const double lx, const double ly, const double lz,
                    const double lcut2,
                    const double strain,
                    double * pos)
{
    int i, j, k;
    int ib, jb;
    int ipair;
    int npair_ah, pair_acc;
    int indent;
    
    double dx, dy, dz;
    double s2, s;
    double dist[3];
    
    int * aux;
    double * aux_r;
    double * aux_e;
    
    aux = NULL;
    aux_r = NULL;
    aux_e = NULL;
    
    lub_sys->tot_pairs = 0;
    lub_sys->index = calloc(np * (np + 1) / 2, sizeof(int));
    
    for (i = 0; i < np; i++) {
        lub_sys->pair_sys[i].bh_count = 0;
        lub_sys->pair_sys[i].bh_list = malloc(lub_sys->pair_sys[i].bh_count * sizeof(int));
        
        lub_sys->pair_sys[i].ah_count = 0;
        lub_sys->pair_sys[i].ah_list = malloc(lub_sys->pair_sys[i].ah_count * sizeof(int));
        lub_sys->pair_sys[i].r = malloc(lub_sys->pair_sys[i].ah_count * sizeof(double));
        lub_sys->pair_sys[i].e = malloc(lub_sys->pair_sys[i].ah_count * 3 * sizeof(double));
    } // for i
    
    for (ib = 0; ib < lub_sys->nd; ib++) {
        i = lub_sys->head[ib];
        // loop over all particles in current box
        while (i >= 0) {
            
            // loop over all particles below i in the current box
            j = lub_sys->list[i];
            while (j >= 0) {
                dist[0] = pos[j * 3 + 0] - pos[i * 3 + 0];
                dist[1] = pos[j * 3 + 1] - pos[i * 3 + 1];
                dist[2] = pos[j * 3 + 2] - pos[i * 3 + 2];
                
                min_image(lx, ly, lz,
                          xlh, xrh,
                          ylh, yrh,
                          zlh, zrh,
                          strain,
                          dist);
                
                s2 = dist[0] * dist[0] + dist[1] * dist[1] + dist[2] * dist[2];
                
                if (s2 < lcut2) {
                    lub_sys->tot_pairs += 1;
                    s = sqrt(s2);
                    if (s <= 2.0001) {
                        s = 2.0001;
                    }
                    
                    if (i > j) {
                        // pair afterhand
                        lub_sys->pair_sys[j].ah_count += 1;
                        lub_sys->pair_sys[j].ah_list = realloc(lub_sys->pair_sys[j].ah_list, lub_sys->pair_sys[j].ah_count * sizeof(int));
                        lub_sys->pair_sys[j].r = realloc(lub_sys->pair_sys[j].r, lub_sys->pair_sys[j].ah_count * sizeof(double));
                        lub_sys->pair_sys[j].e = realloc(lub_sys->pair_sys[j].e, lub_sys->pair_sys[j].ah_count * 3 * sizeof(double));
                        
                        lub_sys->pair_sys[j].ah_list[lub_sys->pair_sys[j].ah_count - 1] = i;
                        lub_sys->pair_sys[j].r[lub_sys->pair_sys[j].ah_count - 1] = s;
                        lub_sys->pair_sys[j].e[(lub_sys->pair_sys[j].ah_count - 1) * 3 + 0] = - dist[0] / s;
                        lub_sys->pair_sys[j].e[(lub_sys->pair_sys[j].ah_count - 1) * 3 + 1] = - dist[1] / s;
                        lub_sys->pair_sys[j].e[(lub_sys->pair_sys[j].ah_count - 1) * 3 + 2] = - dist[2] / s;
                        
                        // pair beforehand
                        lub_sys->pair_sys[i].bh_count += 1;
                        lub_sys->pair_sys[i].bh_list = realloc(lub_sys->pair_sys[i].bh_list, lub_sys->pair_sys[i].bh_count * sizeof(int));
                        
                        lub_sys->pair_sys[i].bh_list[lub_sys->pair_sys[i].bh_count - 1] = j;
                        
                    } else {
                        // pair afterhand
                        lub_sys->pair_sys[i].ah_count += 1;
                        lub_sys->pair_sys[i].ah_list = realloc(lub_sys->pair_sys[i].ah_list, lub_sys->pair_sys[i].ah_count * sizeof(int));
                        lub_sys->pair_sys[i].r = realloc(lub_sys->pair_sys[i].r, lub_sys->pair_sys[i].ah_count * sizeof(double));
                        lub_sys->pair_sys[i].e = realloc(lub_sys->pair_sys[i].e, lub_sys->pair_sys[i].ah_count * 3 * sizeof(double));
                        
                        lub_sys->pair_sys[i].ah_list[lub_sys->pair_sys[i].ah_count - 1] = j;
                        lub_sys->pair_sys[i].r[lub_sys->pair_sys[i].ah_count - 1] = s;
                        lub_sys->pair_sys[i].e[(lub_sys->pair_sys[i].ah_count - 1) * 3 + 0] = dist[0] / s;
                        lub_sys->pair_sys[i].e[(lub_sys->pair_sys[i].ah_count - 1) * 3 + 1] = dist[1] / s;
                        lub_sys->pair_sys[i].e[(lub_sys->pair_sys[i].ah_count - 1) * 3 + 2] = dist[2] / s;
                        
                        // pair beforehand
                        lub_sys->pair_sys[j].bh_count += 1;
                        lub_sys->pair_sys[j].bh_list = realloc(lub_sys->pair_sys[j].bh_list, lub_sys->pair_sys[j].bh_count * sizeof(int));
                        
                        lub_sys->pair_sys[j].bh_list[lub_sys->pair_sys[j].bh_count - 1] = i;
                    } // if i > j ...
                    
                } // if s2 < lcut2
                
                j = lub_sys->list[j];
            } // while j >= 0
            
            // loop over neighbor boxes
            for (k = 0; k < lub_sys->cell_sys[ib].tot_boxes; k++) {
                
                jb = lub_sys->cell_sys[ib].inter_box[k];
                if (jb >= 0) {
                    // loop over all particles in neighbor boxes
                    j = lub_sys->head[jb];
                    
                    while (j >= 0) {
                        dist[0] = pos[j * 3 + 0] - pos[i * 3 + 0];
                        dist[1] = pos[j * 3 + 1] - pos[i * 3 + 1];
                        dist[2] = pos[j * 3 + 2] - pos[i * 3 + 2];
                        
                        min_image(lx, ly, lz,
                                  xlh, xrh,
                                  ylh, yrh,
                                  zlh, zrh,
                                  strain,
                                  dist);
                        s2 = dist[0] * dist[0] + dist[1] * dist[1] + dist[2] * dist[2];
                        
                        if (s2 < lcut2) {
                            lub_sys->tot_pairs += 1;
                            s = sqrt(s2);
                            if (s <= 2.0001) {
                                s = 2.0001;
                            }
                            
                            if (i > j) {
                                // pair afterhand
                                lub_sys->pair_sys[j].ah_count += 1;
                                lub_sys->pair_sys[j].ah_list = realloc(lub_sys->pair_sys[j].ah_list, lub_sys->pair_sys[j].ah_count * sizeof(int));
                                lub_sys->pair_sys[j].r = realloc(lub_sys->pair_sys[j].r, lub_sys->pair_sys[j].ah_count * sizeof(double));
                                lub_sys->pair_sys[j].e = realloc(lub_sys->pair_sys[j].e, lub_sys->pair_sys[j].ah_count * 3 * sizeof(double));
                                
                                lub_sys->pair_sys[j].ah_list[lub_sys->pair_sys[j].ah_count - 1] = i;
                                lub_sys->pair_sys[j].r[lub_sys->pair_sys[j].ah_count - 1] = s;
                                lub_sys->pair_sys[j].e[(lub_sys->pair_sys[j].ah_count - 1) * 3 + 0] = - dist[0] / s;
                                lub_sys->pair_sys[j].e[(lub_sys->pair_sys[j].ah_count - 1) * 3 + 1] = - dist[1] / s;
                                lub_sys->pair_sys[j].e[(lub_sys->pair_sys[j].ah_count - 1) * 3 + 2] = - dist[2] / s;
                                
                                // pair beforehand
                                lub_sys->pair_sys[i].bh_count += 1;
                                lub_sys->pair_sys[i].bh_list = realloc(lub_sys->pair_sys[i].bh_list, lub_sys->pair_sys[i].bh_count * sizeof(int));
                                
                                lub_sys->pair_sys[i].bh_list[lub_sys->pair_sys[i].bh_count - 1] = j;
                                
                            } else {
                                // pair afterhand
                                lub_sys->pair_sys[i].ah_count += 1;
                                lub_sys->pair_sys[i].ah_list = realloc(lub_sys->pair_sys[i].ah_list, lub_sys->pair_sys[i].ah_count * sizeof(int));
                                lub_sys->pair_sys[i].r = realloc(lub_sys->pair_sys[i].r, lub_sys->pair_sys[i].ah_count * sizeof(double));
                                lub_sys->pair_sys[i].e = realloc(lub_sys->pair_sys[i].e, lub_sys->pair_sys[i].ah_count * 3 * sizeof(double));
                                
                                lub_sys->pair_sys[i].ah_list[lub_sys->pair_sys[i].ah_count - 1] = j;
                                lub_sys->pair_sys[i].r[lub_sys->pair_sys[i].ah_count - 1] = s;
                                lub_sys->pair_sys[i].e[(lub_sys->pair_sys[i].ah_count - 1) * 3 + 0] = dist[0] / s;
                                lub_sys->pair_sys[i].e[(lub_sys->pair_sys[i].ah_count - 1) * 3 + 1] = dist[1] / s;
                                lub_sys->pair_sys[i].e[(lub_sys->pair_sys[i].ah_count - 1) * 3 + 2] = dist[2] / s;
                                
                                // pair beforehand
                                lub_sys->pair_sys[j].bh_count += 1;
                                lub_sys->pair_sys[j].bh_list = realloc(lub_sys->pair_sys[j].bh_list, lub_sys->pair_sys[j].bh_count * sizeof(int));
                                
                                lub_sys->pair_sys[j].bh_list[lub_sys->pair_sys[j].bh_count - 1] = i;
                            } // if i > j ...
                            
                        } // if s2 < lcut2
                        
                        j = lub_sys->list[j];
                    } // while j >= 0
                    
                } // if jb > 0
                
            } // for k
            
            i = lub_sys->list[i];
        } // while i >= 0
        
    } // for ib
    
    /* Sort pairs */
    for (i = 0; i < np; i++) {
        if (lub_sys->pair_sys[i].ah_count < 2) {
            continue;
        }
        
        aux = calloc(lub_sys->pair_sys[i].ah_count, sizeof(int));
        aux_r = calloc(lub_sys->pair_sys[i].ah_count, sizeof(double));
        aux_e = calloc(lub_sys->pair_sys[i].ah_count * 3, sizeof(double));
        
        merge_sort_pair_all(0, lub_sys->pair_sys[i].ah_count - 1,
                            lub_sys->pair_sys[i].ah_list, lub_sys->pair_sys[i].r, lub_sys->pair_sys[i].e,
                            aux, aux_r, aux_e);
        
        free(aux);
        free(aux_r);
        free(aux_e);
    } // for i
    
    for (i = 0; i < np; i++) {
        if (lub_sys->pair_sys[i].bh_count < 2) {
            continue;
        }
        
        aux = calloc(lub_sys->pair_sys[i].bh_count, sizeof(int));
        
        merge_sort_pair(0, lub_sys->pair_sys[i].bh_count - 1,
                        lub_sys->pair_sys[i].bh_list,
                        aux);
        
        free(aux);
    } // for i
    
    /* Make list of interacting pairs */
    lub_sys->row_id = malloc(lub_sys->tot_pairs * sizeof(int));
    lub_sys->col_id = malloc(lub_sys->tot_pairs * sizeof(int));
    
    pair_acc = 0;
    for (i = 0; i < np; i++) {
        npair_ah = lub_sys->pair_sys[i].ah_count;
        indent = i * (i + 1) / 2;
        
        lub_sys->index[i * np + i - indent] = i;
        
        for (ipair = 0; ipair < npair_ah; ipair++) {
            j = lub_sys->pair_sys[i].ah_list[ipair];
            
            lub_sys->row_id[ipair + pair_acc] = i;
            lub_sys->col_id[ipair + pair_acc] = j;
            
            lub_sys->index[i * np + j - indent] = np + pair_acc + ipair;
        } // for ipair
        
        pair_acc += npair_ah;
    } // for i
    
}

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
                 double * xm, double * ym, double * zm)
{
    double xi;
    double xi1;
    double dlx;
    double xdlx;
    double dlx1;
    
    int ida;
    int ib;
    int ia;
    
    double c1, c2;
    
    xi = r - 2.0;
    
    if (r <= 2.1) {
        xi1 = 1.0 / xi;
        dlx = log(xi1);
        xdlx = xi * dlx;
        dlx1 = dlx + xdlx;
        
        (* x11a) = 0.25 * xi1 + 0.225 * dlx - 1.23041 + 3.0 / 112.0 * xdlx + 1.8918 * xi;
        (* x12a) = -(* x11a) + 0.00312 - 0.0011 * xi;
        (* y11a) = 1.0 / 6.0 * dlx - 0.39394 + 0.95665 * xi;
        (* y12a) = -(* y11a) + 0.00463606 - 0.007049 * xi;
        (* y11b) = -1.0 / 6.0 * dlx + 0.408286 - 1.0 / 12.0 * xdlx - 0.84055 * xi;
        (* y12b) = -(* y11b) + 0.00230818 - 0.007508 * xi;
        (* x11c) = 0.0479 - 1.0 / 6.0 * xdlx + 0.12494 * xi;
        (* x12c) = -0.031031 + 1.0 / 6.0 * xdlx - 0.174476 * xi;
        (* y11c) = 4.0 * 1.0 / 15.0 * dlx - 0.605434 + 94.0 / 375.0 * xdlx + 0.939139 * xi;
        (* y12c) = 1.0 / 15.0 * dlx - 0.212032 + 31.0 / 375.0 * xdlx + 0.452843 * xi;
        (* x11g) = 0.25 * xi1 + 0.225 * dlx + 39.0 / 280.0 * xdlx - 1.16897 + 1.47882 * xi;
        (* x12g) = -(0.25 * xi1 + 0.225 * dlx + 39.0 / 280.0 * xdlx) + 1.178967 - 1.480493 * xi;
        (* y11g) = 1.0 / 12.0 * dlx + 1.0 / 24.0 * xdlx - 0.2041 + 0.442226 * xi;
        (* y12g) = -(1.0 / 12.0 * dlx + 1.0 / 24.0 * xdlx) + 0.216365 - 0.46983 * xi;
        (* y11h) = 0.5 * (1.0 / 15.0 * dlx) - 0.143777 + 137.0 / 1500.0 * xdlx + 0.264207 * xi;
        (* y12h) = 2.0 * (1.0 / 15.0 * dlx) - 0.298166 + 113.0 / 1500.0 * xdlx + 0.534123 * xi;
        (* xm) = 1.0 / 3.0 * xi1 + 0.3 * dlx - 1.48163 + 0.335714 * xdlx + 1.413604 * xi;
        (* ym) = 1.0 / 6.0 * dlx1 - 0.423489 + 0.827286 * xi;
        (* zm) = 0.0129151 - 0.042284 * xi;
    } else {
        ida = floor(20.0 * xi);
        ib = ida - 2;
        ia = ib + 1;
        c1 = (r - rabc[ib]) / (rabc[ia] - rabc[ib]);
        (* x11a) = (x11as[ia] - x11as[ib]) * c1 + x11as[ib];
        (* x12a) = (x12as[ia] - x12as[ib]) * c1 + x12as[ib];
        (* y11a) = (y11as[ia] - y11as[ib]) * c1 + y11as[ib];
        (* y12a) = (y12as[ia] - y12as[ib]) * c1 + y12as[ib];
        (* y11b) = (y11bs[ia] - y11bs[ib]) * c1 + y11bs[ib];
        (* y12b) = (y12bs[ia] - y12bs[ib]) * c1 + y12bs[ib];
        (* x11c) = (x11cs[ia] - x11cs[ib]) * c1 + x11cs[ib];
        (* x12c) = (x12cs[ia] - x12cs[ib]) * c1 + x12cs[ib];
        (* y11c) = (y11cs[ia] - y11cs[ib]) * c1 + y11cs[ib];
        (* y12c) = (y12cs[ia] - y12cs[ib]) * c1 + y12cs[ib];
        
        if (r < 2.2) {
            ib = -10 + floor(100.0 * xi);
        } else {
            ib = ida + 6;
        }
        ia = ib + 1;
        
        c2 = (r - rghm[ib]) / (rghm[ia] - rghm[ib]);
        (* x11g) = (x11gs[ia] - x11gs[ib]) * c2 + x11gs[ib];
        (* x12g) = (x12gs[ia] - x12gs[ib]) * c2 + x12gs[ib];
        (* y11g) = (y11gs[ia] - y11gs[ib]) * c2 + y11gs[ib];
        (* y12g) = (y12gs[ia] - y12gs[ib]) * c2 + y12gs[ib];
        (* y11h) = (y11hs[ia] - y11hs[ib]) * c2 + y11hs[ib];
        (* y12h) = (y12hs[ia] - y12hs[ib]) * c2 + y12hs[ib];
        (* xm) = (xms[ia] - xms[ib]) * c2 + xms[ib];
        (* ym) = (yms[ia] - yms[ib]) * c2 + yms[ib];
        (* zm) = (zms[ia] - zms[ib]) * c2 + zms[ib];
    }
    
}

void lub_abc(double x11a, double y11a,
             double y11b,
             double x11c, double y11c,
             double x12a, double y12a,
             double y12b,
             double x12c, double y12c,
             double * e,
             double * labc)
{
    int i, j, k;
    double bdot;
    
    for (i = 0; i < 3; i++) {
        for (j = 0; j < 3; j++) {
            labc[(i + 0) * 12 + (j + 0)] = x11a * e[i] * e[j] + y11a * (kron[i * 3 + j] - e[i] * e[j]);     // a11
            bdot = 0.0;
            for (k = 0; k < 3; k++) {
                bdot += y11b * levi[i * 9 + j * 3 + k] * e[k];
            }
            labc[(i + 3) * 12 + (j + 0)] = +bdot;       // b11
            labc[(i + 0) * 12 + (j + 3)] = -bdot;       // bt11
            labc[(i + 3) * 12 + (j + 3)] = x11c * e[i] * e[j] + y11c * (kron[i * 3 + j] - e[i] * e[j]);     // c11
            
            labc[(i + 0) * 12 + (j + 6)] = x12a * e[i] * e[j] + y12a * (kron[i * 3 + j] - e[i] * e[j]);     // a12
            bdot = 0.0;
            for (k = 0; k < 3; k++) {
                bdot += y12b * levi[i * 9 + j * 3 + k] * e[k];
            }
            labc[(i + 3) * 12 + (j + 6)] = +bdot;       // b12
            labc[(i + 0) * 12 + (j + 9)] = +bdot;       // bt12
            labc[(i + 3) * 12 + (j + 9)] = x12c * e[i] * e[j] + y12c * (kron[i * 3 + j] - e[i] * e[j]);     // c12
            
            // a21 = a12, b21 = -b12, bt21 = -bt12, c21 = c12
            labc[(i + 6) * 12 + (j + 0)] = +labc[(i + 0) * 12 + (j + 6)];
            labc[(i + 9) * 12 + (j + 0)] = -labc[(i + 3) * 12 + (j + 6)];
            labc[(i + 6) * 12 + (j + 3)] = -labc[(i + 0) * 12 + (j + 9)];
            labc[(i + 9) * 12 + (j + 3)] = +labc[(i + 3) * 12 + (j + 9)];
            
            // a22 = a11, b22 = -b11, bt22 = -bt11, c22 = c11
            labc[(i + 6) * 12 + (j + 6)] = +labc[(i + 0) * 12 + (j + 0)];
            labc[(i + 9) * 12 + (j + 6)] = -labc[(i + 3) * 12 + (j + 0)];
            labc[(i + 6) * 12 + (j + 9)] = -labc[(i + 0) * 12 + (j + 3)];
            labc[(i + 9) * 12 + (j + 9)] = +labc[(i + 3) * 12 + (j + 3)];
        }
    }
    
}

void lub_gh(double x11g, double y11g,
            double y11h,
            double x12g, double y12g,
            double y12h,
            double * e,
            double * lgh)
{
    int i, j, k, l;
    double hdot1, hdot2;
    double *g = calloc(27, sizeof(double));
    double *h = calloc(27, sizeof(double));
    
    // g11 and h11
    for (i = 0; i < 3; i++) {
        for (j = 0; j < 3; j++) {
            for (k = 0; k < 3; k++) {
                g[i * 9 + j * 3 + k] = x11g * (e[i] * e[j] - 1.0 / 3.0 * kron[i * 3 + j]) * e[k] + y11g * (e[i] * kron[j * 3 + k] + e[j] * kron[i * 3 + k] - 2.0 * e[i] * e[j] * e[k]);
                hdot1 = 0.0;    hdot2 = 0.0;
                for (l = 0; l < 3; l++) {
                    hdot1 += y11h * levi[j * 9 + k * 3 + l] * e[l];
                    hdot2 += y11h * levi[i * 9 + k * 3 + l] * e[l];
                }
                h[i * 9 + j * 3 + k] = hdot1 * e[i] + hdot2 * e[j];
            } // for k
        } // for j
    } // for i
    
    for (i = 0; i < 3; i++) {
        lgh[0 * 12 + i] = g[9 * 0 + 3 * 0 + i];
        lgh[1 * 12 + i] = g[9 * 0 + 3 * 1 + i];
        lgh[2 * 12 + i] = g[9 * 0 + 3 * 2 + i];
        lgh[3 * 12 + i] = g[9 * 1 + 3 * 2 + i];
        lgh[4 * 12 + i] = g[9 * 1 + 3 * 1 + i];
        
        lgh[0 * 12 + i + 3] = h[9 * 0 + 3 * 0 + i];
        lgh[1 * 12 + i + 3] = h[9 * 0 + 3 * 1 + i];
        lgh[2 * 12 + i + 3] = h[9 * 0 + 3 * 2 + i];
        lgh[3 * 12 + i + 3] = h[9 * 1 + 3 * 2 + i];
        lgh[4 * 12 + i + 3] = h[9 * 1 + 3 * 1 + i];
    }
    
    // g12 and h12
    for (i = 0; i < 3; i++) {
        for (j = 0; j < 3; j++) {
            for (k = 0; k < 3; k++) {
                g[i * 9 + j * 3 + k] = x12g * (e[i] * e[j] - 1.0 / 3.0 * kron[3 * i + j]) * e[k] + y12g * (e[i] * kron[j * 3 + k] + e[j] * kron[i * 3 + k] - 2.0 * e[i] * e[j] * e[k]);
                
                hdot1 = 0.0;
                hdot2 = 0.0;
                for (l = 0; l < 3; l++) {
                    hdot1 += y12h * levi[9 * j + 3 * k + l] * e[l];
                    hdot2 += y12h * levi[9 * i + 3 * k + l] * e[l];
                }
                h[i * 9 + j * 3 + k] = hdot1 * e[i] + hdot2 * e[j];
            }
        }
    }
    
    for (i = 0; i < 3; i++) {
        lgh[0 * 12 + i + 6] = g[9 * 0 + 3 * 0 + i];
        lgh[1 * 12 + i + 6] = g[9 * 0 + 3 * 1 + i];
        lgh[2 * 12 + i + 6] = g[9 * 0 + 3 * 2 + i];
        lgh[3 * 12 + i + 6] = g[9 * 1 + 3 * 2 + i];
        lgh[4 * 12 + i + 6] = g[9 * 1 + 3 * 1 + i];
        
        lgh[0 * 12 + i + 9] = h[9 * 0 + 3 * 0 + i];
        lgh[1 * 12 + i + 9] = h[9 * 0 + 3 * 1 + i];
        lgh[2 * 12 + i + 9] = h[9 * 0 + 3 * 2 + i];
        lgh[3 * 12 + i + 9] = h[9 * 1 + 3 * 2 + i];
        lgh[4 * 12 + i + 9] = h[9 * 1 + 3 * 1 + i];
    }
    
    // g22 = -g11, g21 = -g12, h22 = h11, h21 = h12
    for (i = 0; i < 5; i++) {
        for (j = 0; j < 3; j++) {
            lgh[(i + 5) * 12 + j + 6] = -lgh[i * 12 + j + 0];
            lgh[(i + 5) * 12 + j + 0] = -lgh[i * 12 + j + 6];
            lgh[(i + 5) * 12 + j + 9] = +lgh[i * 12 + j + 3];
            lgh[(i + 5) * 12 + j + 3] = +lgh[i * 12 + j + 9];
        }
    }
    
    free(g);
    free(h);
}

void lub_zm(double xm, double ym, double zm,
            double * e,
            double * lzm)
{
    int i, j, k, l;
    double *m = calloc(81, sizeof(double));
    
    for (i = 0; i < 3; i++) {
        for (j = 0; j < 3; j++) {
            for (k = 0; k < 3; k++) {
                for (l = 0; l < 3; l++) {
                    m[i * 27 + j * 9 + k * 3 + l] = 1.5 * xm * (e[i] * e[j] - 1.0 / 3.0 * kron[i * 3 + j]) * (e[k] * e[l] - 1.0 / 3.0 * kron[k * 3 + l]) + 0.5 * ym * (e[i] * kron[j * 3 + l] * e[k] + e[j] * kron[i * 3 + l] * e[k] + e[i] * kron[j * 3 + k] * e[l] + e[j] * kron[i * 3 + k] * e[l] - 4.0 * e[i] * e[j] * e[k] * e[l]) + 0.5 * zm * (kron[i * 3 + k] * kron[j * 3 + l] + kron[j * 3 + k] * kron[i * 3 + l] - kron[i * 3 + j] * kron[k * 3 + l] + e[i] * e[j] * kron[k * 3 + l] + kron[i * 3 + j] * e[k] * e[l] + e[i] * e[j] * e[k] * e[l] - e[i] * kron[j * 3 + l] * e[k] - e[j] * kron[i * 3 + l] * e[k] - e[i] * kron[j * 3 + k] * e[l] - e[j] * kron[i * 3 + k] * e[l]);
                } // for l
            } // for k
        } // for j
    } // for i
    
    lzm[0 * 10 + 0] = m[0 * 27 + 0 * 9 + 0 * 3 + 0];
    lzm[0 * 10 + 1] = m[0 * 27 + 0 * 9 + 0 * 3 + 1];
    lzm[0 * 10 + 2] = m[0 * 27 + 0 * 9 + 0 * 3 + 2];
    lzm[0 * 10 + 3] = m[0 * 27 + 0 * 9 + 1 * 3 + 2];
    lzm[0 * 10 + 4] = m[0 * 27 + 0 * 9 + 1 * 3 + 1];
    
    lzm[1 * 10 + 0] = m[0 * 27 + 1 * 9 + 0 * 3 + 0];
    lzm[1 * 10 + 1] = m[0 * 27 + 1 * 9 + 0 * 3 + 1];
    lzm[1 * 10 + 2] = m[0 * 27 + 1 * 9 + 0 * 3 + 2];
    lzm[1 * 10 + 3] = m[0 * 27 + 1 * 9 + 1 * 3 + 2];
    lzm[1 * 10 + 4] = m[0 * 27 + 1 * 9 + 1 * 3 + 1];
    
    lzm[2 * 10 + 0] = m[0 * 27 + 2 * 9 + 0 * 3 + 0];
    lzm[2 * 10 + 1] = m[0 * 27 + 2 * 9 + 0 * 3 + 1];
    lzm[2 * 10 + 2] = m[0 * 27 + 2 * 9 + 0 * 3 + 2];
    lzm[2 * 10 + 3] = m[0 * 27 + 2 * 9 + 1 * 3 + 2];
    lzm[2 * 10 + 4] = m[0 * 27 + 2 * 9 + 1 * 3 + 1];
    
    lzm[3 * 10 + 0] = m[1 * 27 + 2 * 9 + 0 * 3 + 0];
    lzm[3 * 10 + 1] = m[1 * 27 + 2 * 9 + 0 * 3 + 1];
    lzm[3 * 10 + 2] = m[1 * 27 + 2 * 9 + 0 * 3 + 2];
    lzm[3 * 10 + 3] = m[1 * 27 + 2 * 9 + 1 * 3 + 2];
    lzm[3 * 10 + 4] = m[1 * 27 + 2 * 9 + 1 * 3 + 1];
    
    lzm[4 * 10 + 0] = m[1 * 27 + 1 * 9 + 0 * 3 + 0];
    lzm[4 * 10 + 1] = m[1 * 27 + 1 * 9 + 0 * 3 + 1];
    lzm[4 * 10 + 2] = m[1 * 27 + 1 * 9 + 0 * 3 + 2];
    lzm[4 * 10 + 3] = m[1 * 27 + 1 * 9 + 1 * 3 + 2];
    lzm[4 * 10 + 4] = m[1 * 27 + 1 * 9 + 1 * 3 + 1];
    
    // m12 = m11
    lzm[0 * 10 + 5] = lzm[0 * 10 + 0];
    lzm[0 * 10 + 6] = lzm[0 * 10 + 1];
    lzm[0 * 10 + 7] = lzm[0 * 10 + 2];
    lzm[0 * 10 + 8] = lzm[0 * 10 + 3];
    lzm[0 * 10 + 9] = lzm[0 * 10 + 4];
    
    lzm[1 * 10 + 5] = lzm[1 * 10 + 0];
    lzm[1 * 10 + 6] = lzm[1 * 10 + 1];
    lzm[1 * 10 + 7] = lzm[1 * 10 + 2];
    lzm[1 * 10 + 8] = lzm[1 * 10 + 3];
    lzm[1 * 10 + 9] = lzm[1 * 10 + 4];
    
    lzm[2 * 10 + 5] = lzm[2 * 10 + 0];
    lzm[2 * 10 + 6] = lzm[2 * 10 + 1];
    lzm[2 * 10 + 7] = lzm[2 * 10 + 2];
    lzm[2 * 10 + 8] = lzm[2 * 10 + 3];
    lzm[2 * 10 + 9] = lzm[2 * 10 + 4];
    
    lzm[3 * 10 + 5] = lzm[3 * 10 + 0];
    lzm[3 * 10 + 6] = lzm[3 * 10 + 1];
    lzm[3 * 10 + 7] = lzm[3 * 10 + 2];
    lzm[3 * 10 + 8] = lzm[3 * 10 + 3];
    lzm[3 * 10 + 9] = lzm[3 * 10 + 4];
    
    lzm[4 * 10 + 5] = lzm[4 * 10 + 0];
    lzm[4 * 10 + 6] = lzm[4 * 10 + 1];
    lzm[4 * 10 + 7] = lzm[4 * 10 + 2];
    lzm[4 * 10 + 8] = lzm[4 * 10 + 3];
    lzm[4 * 10 + 9] = lzm[4 * 10 + 4];
    
    // m21 = m12
    lzm[5 * 10 + 0] = lzm[0 * 10 + 5];
    lzm[5 * 10 + 1] = lzm[0 * 10 + 6];
    lzm[5 * 10 + 2] = lzm[0 * 10 + 7];
    lzm[5 * 10 + 3] = lzm[0 * 10 + 8];
    lzm[5 * 10 + 4] = lzm[0 * 10 + 9];
    
    lzm[6 * 10 + 0] = lzm[1 * 10 + 5];
    lzm[6 * 10 + 1] = lzm[1 * 10 + 6];
    lzm[6 * 10 + 2] = lzm[1 * 10 + 7];
    lzm[6 * 10 + 3] = lzm[1 * 10 + 8];
    lzm[6 * 10 + 4] = lzm[1 * 10 + 9];
    
    lzm[7 * 10 + 0] = lzm[2 * 10 + 5];
    lzm[7 * 10 + 1] = lzm[2 * 10 + 6];
    lzm[7 * 10 + 2] = lzm[2 * 10 + 7];
    lzm[7 * 10 + 3] = lzm[2 * 10 + 8];
    lzm[7 * 10 + 4] = lzm[2 * 10 + 9];
    
    lzm[8 * 10 + 0] = lzm[3 * 10 + 5];
    lzm[8 * 10 + 1] = lzm[3 * 10 + 6];
    lzm[8 * 10 + 2] = lzm[3 * 10 + 7];
    lzm[8 * 10 + 3] = lzm[3 * 10 + 8];
    lzm[8 * 10 + 4] = lzm[3 * 10 + 9];
    
    lzm[9 * 10 + 0] = lzm[4 * 10 + 5];
    lzm[9 * 10 + 1] = lzm[4 * 10 + 6];
    lzm[9 * 10 + 2] = lzm[4 * 10 + 7];
    lzm[9 * 10 + 3] = lzm[4 * 10 + 8];
    lzm[9 * 10 + 4] = lzm[4 * 10 + 9];
    
    // m22 = m11
    lzm[5 * 10 + 5] = lzm[0 * 10 + 0];
    lzm[5 * 10 + 6] = lzm[0 * 10 + 1];
    lzm[5 * 10 + 7] = lzm[0 * 10 + 2];
    lzm[5 * 10 + 8] = lzm[0 * 10 + 3];
    lzm[5 * 10 + 9] = lzm[0 * 10 + 4];
    
    lzm[6 * 10 + 5] = lzm[1 * 10 + 0];
    lzm[6 * 10 + 6] = lzm[1 * 10 + 1];
    lzm[6 * 10 + 7] = lzm[1 * 10 + 2];
    lzm[6 * 10 + 8] = lzm[1 * 10 + 3];
    lzm[6 * 10 + 9] = lzm[1 * 10 + 4];
    
    lzm[7 * 10 + 5] = lzm[2 * 10 + 0];
    lzm[7 * 10 + 6] = lzm[2 * 10 + 1];
    lzm[7 * 10 + 7] = lzm[2 * 10 + 2];
    lzm[7 * 10 + 8] = lzm[2 * 10 + 3];
    lzm[7 * 10 + 9] = lzm[2 * 10 + 4];
    
    lzm[8 * 10 + 5] = lzm[3 * 10 + 0];
    lzm[8 * 10 + 6] = lzm[3 * 10 + 1];
    lzm[8 * 10 + 7] = lzm[3 * 10 + 2];
    lzm[8 * 10 + 8] = lzm[3 * 10 + 3];
    lzm[8 * 10 + 9] = lzm[3 * 10 + 4];
    
    lzm[9 * 10 + 5] = lzm[4 * 10 + 0];
    lzm[9 * 10 + 6] = lzm[4 * 10 + 1];
    lzm[9 * 10 + 7] = lzm[4 * 10 + 2];
    lzm[9 * 10 + 8] = lzm[4 * 10 + 3];
    lzm[9 * 10 + 9] = lzm[4 * 10 + 4];
    
    free(m);
}

/*
 Obtain lubrication tensors
 Input:
 np                  number of particles
 lambda              a parameter that makes Lfu invertible
 */
void lub_mat(struct lub * lub_sys,
             const int np,
             const double lambda)
{
    int nb;
    int ni5, nb5, nb6;
    int ipair;
    int in, jn;
    int npair_ah, pair_acc;
    int pai, npai;
    int i, j;
    int ic, jc, jr;
    
    double x11a, y11a;
    double y11b;
    double x11c, y11c;
    double x11g, y11g;
    double y11h;
    
    double x12a, y12a;
    double y12b;
    double x12c, y12c;
    double x12g, y12g;
    double y12h;
    
    double xm, ym, zm;
    
    double labc[144];
    double lgh[120];
    double lzm[100];
    
    nb = np + lub_sys->tot_pairs;
    ni5 = lub_sys->tot_pairs * 5;
    nb5 = nb * 5;
    nb6 = nb * 6;
    
    lub_sys->lfu = calloc(6 * nb6, sizeof(double));
    lub_sys->lfe = calloc(6 * ni5, sizeof(double));
    lub_sys->lsu = calloc(5 * nb6, sizeof(double));
    lub_sys->lse = calloc(5 * nb5, sizeof(double));
    
    pair_acc = 0;
    
    for (in = 0; in < np; in++) {
        for (i = 0; i < 6; i++) {
            lub_sys->lfu[in * 6 + i * nb6 + i] += lambda;
        }
        
        npair_ah = lub_sys->pair_sys[in].ah_count;
        
        for (ipair = 0; ipair < npair_ah; ipair++) {
            jn = lub_sys->pair_sys[in].ah_list[ipair];
            
            lub_scalars(lub_sys->pair_sys[in].r[ipair],
                        & x11a, & y11a,
                        & y11b,
                        & x11c, & y11c,
                        & x11g, & y11g,
                        & y11h,
                        & x12a, & y12a,
                        & y12b,
                        & x12c, & y12c,
                        & x12g, & y12g,
                        & y12h,
                        & xm, & ym, & zm);
            
            xm *= 0.5;
            ym *= 0.5;
            zm *= 0.5;
            
            lub_abc(x11a, y11a,
                    y11b,
                    x11c, y11c,
                    x12a, y12a,
                    y12b,
                    x12c, y12c,
                    lub_sys->pair_sys[in].e + ipair * 3,
                    labc);
            lub_gh(x11g, y11g,
                   y11h,
                   x12g, y12g,
                   y12h,
                   lub_sys->pair_sys[in].e + ipair * 3,
                   lgh);
            lub_zm(xm, ym, zm,
                   lub_sys->pair_sys[in].e + ipair * 3,
                   lzm);
            
            pai = pair_acc + ipair;
            npai = np + pai;
            
            // Lfu + Lfe
            for (i = 0; i < 6; i++) {
                // Lfu
                for (j = 0; j < 6; j++) {
                    ic = in * 6 + j;
                    jr = jn * 6 + j;
                    jc = npai * 6 + j;
                    
                    lub_sys->lfu[i * nb6 + ic] += labc[(i + 0) * 12 + (j + 0)];
                    lub_sys->lfu[i * nb6 + jc] += labc[(i + 0) * 12 + (j + 6)];
                    lub_sys->lfu[i * nb6 + jr] += labc[(i + 6) * 12 + (j + 6)];
                }
                // Lfe
                for (j = 0; j < 5; j++) {
                    jc = pai * 5 + j;
                    
                    lub_sys->lfe[i * ni5 + jc] += lgh[(j + 5) * 12 + (i + 0)];
                }
            } // for i
            
            // Lsu + Lse
            for (i = 0; i < 5; i++) {
                // Lsu
                for (j = 0; j < 6; j++) {
                    ic = in * 6 + j;
                    jr = jn * 6 + j;
                    jc = npai * 6 + j;
                    
                    lub_sys->lsu[i * nb6 + ic] += lgh[(i + 0) * 12 + (j + 0)];
                    lub_sys->lsu[i * nb6 + jc] += lgh[(i + 0) * 12 + (j + 6)];
                    lub_sys->lsu[i * nb6 + jr] += lgh[(i + 5) * 12 + (j + 6)];
                }
                // Lse
                for (j = 0; j < 5; j++) {
                    ic = in * 5 + j;
                    jr = jn * 5 + j;
                    jc = npai * 5 + j;
                    
                    lub_sys->lse[i * nb5 + ic] += lzm[(i + 0) * 10 + (j + 0)];
                    lub_sys->lse[i * nb5 + jc] += lzm[(i + 0) * 10 + (j + 5)];
                    lub_sys->lse[i * nb5 + jr] += lzm[(i + 5) * 10 + (j + 5)];
                }
            } // for i
            
        } // for ipair
        
        pair_acc += npair_ah;
        
    } // for in
    
}

/*
 Obtain preconditoners for Lfu and \tilde{M}
 1. The preconditioner of Lfu is constructed by incomplete Cholesky decomposition with zero fill-ins
 2. The grand mobility \tilde{M} is
 \tilde{M}  =   | Muf   Mus | + | Lfu^{-1}  0 | - \lambda * | Muf   Mus | . | Lfu^{-1}  0 |
 | Mef   Mes |   | 0         0 |             | Mef   Mes |   | 0         0 |
 where Lfu is approximated as the block diagonal
 M is approxiamted by the diagonal components
 */
void lub_precond(struct lub * lub_sys,
                 const int np, const double lambda)
{
    int nb;
    int n6;
    int nb6;
    int i, j, k;
    int ipair, jpair;
    int in, jn, kn;
    int indent_i, indent_j, indent_k;
    int ind_i, ind_j, ind_k;
    int block_ind_i, block_ind_j, block_ind_k;
    int npair_ah, npair_bh;
    
    double lii, lij;
    double xa, xc;
    
    double * li;
    
    nb = np + lub_sys->tot_pairs;
    n6 = np * 6;
    nb6 = nb * 6;
    
    /* Use self function of non-Ewald case */
    xa = 1.0;
    xc = 0.75;
    
    lub_sys->iccl = calloc(6 * nb6, sizeof(double));
    lub_sys->lfu_l = calloc(6 * n6, sizeof(double));
    lub_sys->lfu_u = calloc(6 * n6, sizeof(double));
    
    li = calloc(6 * 6, sizeof(double));
    
    lub_sys->precond = 1;
    
    /* Obtain the preconditioner of Lfu */
    for (in = 0; in < np; in++) {
        npair_ah = lub_sys->pair_sys[in].ah_count;
        npair_bh = lub_sys->pair_sys[in].bh_count;
        indent_i = in * (in + 1) / 2;
        
        for (i = 0; i < 6; i++) {
            // Calculate Lii
            lii = 0.0;
            for (ipair = 0; ipair < npair_bh; ipair++) {
                jn = lub_sys->pair_sys[in].bh_list[ipair];
                indent_j = jn * (jn + 1) / 2;
                ind_j = jn * np + in - indent_j;
                block_ind_j = lub_sys->index[ind_j];
                
                for (k = 0; k < 6; k++) {
                    lii += lub_sys->iccl[block_ind_j * 6 + i * nb6 + k] * lub_sys->iccl[block_ind_j * 6 + i * nb6 + k];
                } // for k
                
            } // for ipair
            
            for (k = 0; k < i; k++) {
                lii += lub_sys->iccl[in * 6 + i * nb6 + k] * lub_sys->iccl[in * 6 + i * nb6 + k];
            }
            lub_sys->iccl[in * 6 + i * nb6 + i] = sqrt(lub_sys->lfu[in * 6 + i * nb6 + i] - lii);
            if (lub_sys->iccl[in * 6 + i * nb6 + i] != lub_sys->iccl[in * 6 + i * nb6 + i]) {
                lub_sys->precond = 0;
                break;
            }
            
            // Calculate Lij
            for (j = i + 1; j < 6; j++) {
                
                if (lub_sys->lfu[in * 6 + i * nb6 + j] != 0.0) {
                    lij = 0.0;
                    for (ipair = 0; ipair < npair_bh; ipair++) {
                        jn = lub_sys->pair_sys[in].bh_list[ipair];
                        indent_j = jn * (jn + 1) / 2;
                        ind_j = jn * np + in - indent_j;
                        block_ind_j = lub_sys->index[ind_j];
                        
                        for (k = 0; k < 6; k++) {
                            lij += lub_sys->iccl[block_ind_j * 6 + i * nb6 + k] * lub_sys->iccl[block_ind_j * 6 + j * nb6 + k];
                        }
                    } // for ipair
                    
                    for (k = 0; k < i; k++) {
                        lij += lub_sys->iccl[in * 6 + j * nb6 + k] * lub_sys->iccl[in * 6 + i * nb6 + k];
                    }
                    lub_sys->iccl[in * 6 + j * nb6 + i] = 1.0 / lub_sys->iccl[in * 6 + i * nb6 + i] * (lub_sys->lfu[in * 6 + i * nb6 + j] - lij);
                    
                } // if
                
            } // for j
            
            for (ipair = 0; ipair < npair_ah; ipair++) {
                jn = lub_sys->pair_sys[in].ah_list[ipair];
                ind_i = in * np + jn - indent_i;
                block_ind_i = lub_sys->index[ind_i];
                
                
                for (j = 0; j < 6; j++) {
                    
                    if (lub_sys->lfu[block_ind_i * 6 + i * nb6 + j] != 0.0) {
                        lij = 0.0;
                        for (jpair = 0; jpair < npair_bh; jpair++) {
                            kn = lub_sys->pair_sys[in].bh_list[jpair];
                            indent_k = kn * (kn + 1) / 2;
                            ind_j = kn * np + in - indent_k;
                            ind_k = kn * np + jn - indent_k;
                            
                            block_ind_j = lub_sys->index[ind_j];
                            block_ind_k = lub_sys->index[ind_k];
                            
                            if (block_ind_j != 0 && block_ind_k != 0) {
                                for (k = 0; k < 6; k++) {
                                    lij += lub_sys->iccl[block_ind_j * 6 + i * nb6 + k] * lub_sys->iccl[block_ind_k * 6 + j * nb6 + k];
                                }
                                
                            } // if
                            
                        } // for jpair
                        
                        for (k = 0; k < i; k++) {
                            lij += lub_sys->iccl[block_ind_i * 6 + j * nb6 + k] * lub_sys->iccl[in * 6 + i * nb6 + k];
                        }
                        lub_sys->iccl[block_ind_i * 6 + j * nb6 + i] = 1.0 / lub_sys->iccl[in * 6 + i * nb6 + i] * (lub_sys->lfu[block_ind_i * 6 + i * nb6 + j] - lij);
                        
                    } // if != 0.0
                    
                } // for j
                
            } // for ipair
            
        } // for i
        
    } // for in
    
    /* Obtain the preconditioner of \tilde{M} */
    for (i = 0; i < np; i++) {
        // invert block diagonal of Lfu
        cholesky_inverse(6,
                         & lub_sys->lfu[i * 6], nb6,
                         li, 6);
        
        // calculate Li = Li + M - \lambda * M . Li
        for (j = 0; j < 3; j++) {
            
            for (k = 0; k < 3; k++) {
                li[ j * 6      + k    ] = li[ j      * 6 + k    ] - lambda * xa * li[ j      * 6 + k    ] + xa * kron[j * 3 + k];
                li[ j * 6      + k + 3] = li[ j      * 6 + k + 3] - lambda * xa * li[ j      * 6 + k + 3];
                li[(j + 3) * 6 + k    ] = li[(j + 3) * 6 + k    ] - lambda * xc * li[(j + 3) * 6 + k    ];
                li[(j + 3) * 6 + k + 3] = li[(j + 3) * 6 + k + 3] - lambda * xc * li[(j + 3) * 6 + k + 3] + xc * kron[j * 3 + k];
            } // for k
            
        } // for j
        
        // LU decomposition of Li
        lu_decomp(6,
                  li, 6,
                  & lub_sys->lfu_l[i * 6], n6,
                  & lub_sys->lfu_u[i * 6], n6);
    } // for i
    
    free(li);
    
}

