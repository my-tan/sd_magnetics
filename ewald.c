//
//  ewald.c
//  asd_sphere_20200805
//
//  Created by Mingyang Tan on 8/5/20.
//  Copyright © 2020 Mingyang Tan. All rights reserved.
//

#include "ewald.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <fftw3.h>
#include <cblas.h>
#include <float.h>
#include "sd.h"
#include "cell.h"

extern double * kron;
extern double * levi;

struct ewald * ewald_initialize(struct sd * sys)
{
    int i;
    
    struct ewald * ewald_sys = NULL;
    ewald_sys = (struct ewald * ) malloc(sizeof(struct ewald));
    
    ewald_sys->ndx = floor(sys->lx / sys->rcut);
    ewald_sys->ndy = floor(sys->ly / sys->rcut);
    ewald_sys->ndz = floor(sys->lz / sys->rcut);
    ewald_sys->nd = ewald_sys->ndx * ewald_sys->ndy * ewald_sys->ndz;
    ewald_sys->rx = sys->lx / (double) ewald_sys->ndx;
    ewald_sys->ry = sys->ly / (double) ewald_sys->ndy;
    ewald_sys->rz = sys->lz / (double) ewald_sys->ndz;
    
    ewald_sys->tot_pairs = 0;
    ewald_sys->row_id = NULL;
    ewald_sys->col_id = NULL;
    ewald_sys->index = NULL;
    ewald_sys->head = NULL;
    ewald_sys->list = NULL;
    
    ewald_sys->muf = NULL;
    ewald_sys->mus = NULL;
    ewald_sys->mef = NULL;
    ewald_sys->mes = NULL;
    
    ewald_sys->mbd = NULL;
    
    ewald_sys->cell_sys = (struct mcell * ) malloc(ewald_sys->nd * sizeof(struct mcell));
    for (i = 0; i < ewald_sys->nd; i++) {
        ewald_sys->cell_sys[i].tot_boxes = 0;
        ewald_sys->cell_sys[i].inter_box = NULL;
    }
    
    ewald_sys->gridk_sys = (struct gridK * ) malloc(sys->nk * sizeof(struct gridK));
    for (i = 0; i < sys->nk; i++) {
        ewald_sys->gridk_sys[i].kx = 0.0;
        ewald_sys->gridk_sys[i].ky = 0.0;
        ewald_sys->gridk_sys[i].kz = 0.0;
        ewald_sys->gridk_sys[i].k = 0.0;
        
        ewald_sys->gridk_sys[i].ya = 0.0;
        ewald_sys->gridk_sys[i].yb = 0.0;
        ewald_sys->gridk_sys[i].yc = 0.0;
        ewald_sys->gridk_sys[i].yg = 0.0;
        ewald_sys->gridk_sys[i].yh = 0.0;
        ewald_sys->gridk_sys[i].ym = 0.0;
        
        ewald_sys->gridk_sys[i].w = 0.0;
    }
    
    ewald_sys->gridp_sys = (struct gridP * ) malloc(sys->np * sizeof(struct gridP));
    for (i = 0; i < sys->np; i++) {
        ewald_sys->gridp_sys[i].ids = NULL;
        ewald_sys->gridp_sys[i].r2 = NULL;
        ewald_sys->gridp_sys[i].dx = NULL;
        ewald_sys->gridp_sys[i].dy = NULL;
        ewald_sys->gridp_sys[i].dz = NULL;
    }
    
    ewald_sys->pair_sys = (struct mpair * ) malloc(sys->np * sizeof(struct mpair));
    for (i = 0; i < sys->np; i++) {
        ewald_sys->pair_sys[i].bh_count = 0;
        ewald_sys->pair_sys[i].bh_list = NULL;
        
        ewald_sys->pair_sys[i].ah_count = 0;
        ewald_sys->pair_sys[i].ah_list = NULL;
        
        ewald_sys->pair_sys[i].r = NULL;
        ewald_sys->pair_sys[i].e = NULL;
    }
    
    return ewald_sys;
};

/*
 Calculate the reciprocal functions of each grid node.
 Input:
 lx, ly, lz          dimension of simulation box in x,y,z-directions
 nkx, nky, nkz       number of nodes in x,y,z-directions
 nk                  total number of nodes
 xi2                 xi * xi
 eta                 splitting parameter of Gaussian kernel
 strain              strain unit
 */
void gridk_value(struct ewald * ewald_sys,
                 const double lx, const double ly, const double lz,
                 const int nkx, const int nky, const int nkz,
                 const int nk,
                 const double xi2, const double eta,
                 const double strain)
{
    int i, j, k;
    int id;
    
    double kx, ky, kz;
    double k2, ks;
    double kxi, expkxi;
    
    for (i = 0; i < nkx; i++) {
        if (i >= (nkx + 1) / 2) {
            i -= nkx;
        }
        kx = 2.0 * M_PI / lx * (double) i;
        
        for (j = 0; j < nky; j++) {
            if (j >= (nky + 1) / 2) {
                j -= nky;
            }
            ky = 2.0 * M_PI / ly * (double) j;
            ky -= kx * strain;              // deformation of grid due to shear
            
            for (k = 0; k < nkz; k++) {
                if (k >= (nkz + 1) / 2) {
                    k -= nkz;
                }
                kz = 2.0 * M_PI / lz * (double) k;
                id = k * nky * nkx + j * nkx + i;
                
                if (id == 0) {
                    ewald_sys->gridk_sys[id].kx = 0.0;
                    ewald_sys->gridk_sys[id].ky = 0.0;
                    ewald_sys->gridk_sys[id].kz = 0.0;
                    ewald_sys->gridk_sys[id].k = 0.0;
                    
                    ewald_sys->gridk_sys[id].ya = 0.0;
                    ewald_sys->gridk_sys[id].yb = 0.0;
                    ewald_sys->gridk_sys[id].yc = 0.0;
                    ewald_sys->gridk_sys[id].yg = 0.0;
                    ewald_sys->gridk_sys[id].yh = 0.0;
                    ewald_sys->gridk_sys[id].ym = 0.0;
                    
                    ewald_sys->gridk_sys[id].w = 0.0;
                } else {
                    k2 = kx * kx + ky * ky + kz * kz;
                    ks = sqrt(k2);
                    
                    kxi = k2 / (4.0 * xi2);
                    expkxi = 6.0 * M_PI * (1.0 + kxi * (1.0 + 2.0 * kxi)) / k2 * exp( - (1.0 - eta) * kxi) / (double) nk;
                    
                    ewald_sys->gridk_sys[id].kx = kx / ks;
                    ewald_sys->gridk_sys[id].ky = ky / ks;
                    ewald_sys->gridk_sys[id].kz = kz / ks;
                    ewald_sys->gridk_sys[id].k = ks;
                    
                    ewald_sys->gridk_sys[id].ya = expkxi * (1.0 - k2 / 3.0);
                    ewald_sys->gridk_sys[id].yb = expkxi * 0.5 * ks;
                    ewald_sys->gridk_sys[id].yc = expkxi * 0.25 * k2;
                    ewald_sys->gridk_sys[id].yg = expkxi * ks * (1.0 - 4.0 / 15.0 * k2);
                    ewald_sys->gridk_sys[id].yh = expkxi * 0.25 * k2;
                    ewald_sys->gridk_sys[id].ym = expkxi * 0.5 * k2 * (1.0 - k2 / 5.0);
                    
                    ewald_sys->gridk_sys[id].w = 4.0 * M_PI * exp(- (1.0 - eta) * kxi) / (double) nk;
                } // if ... else ...
                
            } // for k
            
        } // for j
        
    } // for i
    
}

/*
 Identify p**3 closest nodes around each particle, and calculate square of distances between particle and nodes
 Input:
 np                  number of particles
 p                   number of closest nodes in one direction
 p3                  p * p * p
 lx, ly, lz          dimension of simulation box in x,y,z-directions
 xlh                 - lx / 2.0
 xrh                 + lx / 2.0
 ylh                 - ly / 2.0
 yrh                 + ly / 2.0
 zlh                 - lz / 2.0
 zrh                 - lz / 2.0
 nkx, nky, nkz       number of nodes in x,y,z-directions
 dkx, dky, dkz       separation between nodes in x,y,z-directions
 strain              unit of strain
 pos                 position of particles
 */
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
                  const double * pos)
{
    int in;
    int i, j, k;
    int ix, iy, iz;
    int x, y, z;
    int id;
    int pd2;
    int count;
    
    double lxd2, lyd2, lzd2;
    double px, py, pz;
    double kx, ky, kz;
    double r2;
    
    double fpos[3];
    double dist[3];
    
    lxd2 = lx / 2.0;
    lyd2 = ly / 2.0;
    lzd2 = lz / 2.0;
    
    pd2 = (int) p / 2;
    
    for (in = 0; in < np; in++) {
        ewald_sys->gridp_sys[in].ids = calloc(p3, sizeof(int));
        ewald_sys->gridp_sys[in].r2 = calloc(p3, sizeof(double));
        ewald_sys->gridp_sys[in].dx = calloc(p3, sizeof(double));
        ewald_sys->gridp_sys[in].dy = calloc(p3, sizeof(double));
        ewald_sys->gridp_sys[in].dz = calloc(p3, sizeof(double));
        
        pos_fraction(lx, ly, lz,
                     xlh, ylh, zlh,
                     strain,
                     & pos[in * 3],
                     fpos);
        
        px = fpos[0] * (double) nkx;
        py = fpos[1] * (double) nky;
        pz = fpos[2] * (double) nkz;
        
        x = (int) px;
        y = (int) py;
        z = (int) pz;
        
        count = 0;
        for (i = 0; i < p; i++) {
            
            for (j = 0; j < p; j++) {
                
                for (k = 0; k < p; k++) {
                    
                    ix = x + i - pd2 + 1 - p % 2 * (px - (double) x < 0.5);
                    iy = y + j - pd2 + 1 - p % 2 * (py - (double) y < 0.5);
                    iz = z + k - pd2 + 1 - p % 2 * (pz - (double) z < 0.5);
                    
                    ix = (ix < 0) ? ix + nkx : ( (ix >= nkx) ? ix - nkx : ix );
                    iy = (iy < 0) ? iy + nky : ( (iy >= nky) ? iy - nky : iy );
                    iz = (iz < 0) ? iz + nkz : ( (iz >= nkz) ? iz - nkz : iz );
                    
                    kx = (double) ix * dkx - lxd2;
                    ky = (double) iy * dky - lyd2;
                    kz = (double) iz * dkz - lzd2;
                    
                    kx = kx + strain * ky;
                    
                    id = iz * nky * nkx + iy * nkx + ix;
                    
                    dist[0] = kx - pos[in * 3 + 0];
                    dist[1] = ky - pos[in * 3 + 1];
                    dist[2] = kz - pos[in * 3 + 2];
                    
                    min_image(lx, ly, lz,
                              xlh, xrh,
                              ylh, yrh,
                              zlh, zrh,
                              strain,
                              dist);
                    r2 = dist[0] * dist[0] + dist[1] * dist[1] + dist[2] * dist[2];
                    
                    ewald_sys->gridp_sys[in].ids[count] = id;
                    ewald_sys->gridp_sys[in].r2[count] = r2;
                    
                    ewald_sys->gridp_sys[in].dx[count] = dist[0];
                    ewald_sys->gridp_sys[in].dy[count] = dist[1];
                    ewald_sys->gridp_sys[in].dz[count] = dist[2];
                    
                    count++;
                } // for k
                
            } // for j
            
        } // for i
        
    } // for in
    
}

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

void near_rcell_neighbor_shear(struct ewald * ewald_sys,
                               const int np,
                               const double strain)
{
    int ndx, ndy, ndz;
    int px;
    int ix, iy, iz;
    int id;
    
    double sx;
    
    ndx = ewald_sys->ndx;
    ndy = ewald_sys->ndy;
    ndz = ewald_sys->ndz;
    
    for (iz = 0; iz < ndz; iz++) {
        
        for (iy = 0; iy < ndy - 1; iy++) {
            
            for (ix = 0; ix < ndx; ix++) {
                id = cell_id(ix, iy, iz,
                             ndx, ndy, ndz);
                
                ewald_sys->cell_sys[id].tot_boxes = 16;
                ewald_sys->cell_sys[id].inter_box = malloc(16 * sizeof(int));
                
                ewald_sys->cell_sys[id].inter_box[0] = cell_id(ix + 1, iy    , iz,
                                                               ndx, ndy, ndz);
                ewald_sys->cell_sys[id].inter_box[1] = cell_id(ix + 1, iy + 1, iz,
                                                               ndx, ndy, ndz);
                ewald_sys->cell_sys[id].inter_box[2] = cell_id(ix    , iy + 1, iz,
                                                               ndx, ndy, ndz);
                if (ewald_sys->ndx == 2) {
                    ewald_sys->cell_sys[id].inter_box[3] = - 1;
                } else {
                    ewald_sys->cell_sys[id].inter_box[3] = cell_id(ix - 1, iy + 1, iz,
                                                                   ndx, ndy, ndz);
                }
                
                ewald_sys->cell_sys[id].inter_box[4] = cell_id(ix + 1, iy    , iz - 1,
                                                               ndx, ndy, ndz);
                ewald_sys->cell_sys[id].inter_box[5] = cell_id(ix + 1, iy + 1, iz - 1,
                                                               ndx, ndy, ndz);
                ewald_sys->cell_sys[id].inter_box[6] = cell_id(ix    , iy + 1, iz - 1,
                                                               ndx, ndy, ndz);
                if (ewald_sys->ndx == 2) {
                    ewald_sys->cell_sys[id].inter_box[7] = - 1;
                } else {
                    ewald_sys->cell_sys[id].inter_box[7] = cell_id(ix - 1, iy + 1, iz - 1,
                                                                   ndx, ndy, ndz);
                }
                
                if (ewald_sys->ndz == 2) {
                    ewald_sys->cell_sys[id].inter_box[8] = - 1;
                    ewald_sys->cell_sys[id].inter_box[9] = - 1;
                    ewald_sys->cell_sys[id].inter_box[10] = - 1;
                    ewald_sys->cell_sys[id].inter_box[11] = - 1;
                } else {
                    ewald_sys->cell_sys[id].inter_box[8] = cell_id(ix + 1, iy    , iz + 1,
                                                                   ndx, ndy, ndz);
                    ewald_sys->cell_sys[id].inter_box[9] = cell_id(ix + 1, iy + 1, iz + 1,
                                                                   ndx, ndy, ndz);
                    ewald_sys->cell_sys[id].inter_box[10] = cell_id(ix    , iy + 1, iz + 1,
                                                                    ndx, ndy, ndz);
                    ewald_sys->cell_sys[id].inter_box[11] = cell_id(ix - 1, iy + 1, iz + 1,
                                                                    ndx, ndy, ndz);
                }
                
                ewald_sys->cell_sys[id].inter_box[12] = cell_id(ix    , iy    , iz + 1,
                                                                ndx, ndy, ndz);
                
                ewald_sys->cell_sys[id].inter_box[13] = - 1;
                ewald_sys->cell_sys[id].inter_box[14] = - 1;
                ewald_sys->cell_sys[id].inter_box[15] = - 1;
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
            
            ewald_sys->cell_sys[id].tot_boxes = 16;
            ewald_sys->cell_sys[id].inter_box = malloc(16 * sizeof(int));
            
            ewald_sys->cell_sys[id].inter_box[0] = cell_id(ix + 1,      iy    , iz,
                                                           ndx, ndy, ndz);
            ewald_sys->cell_sys[id].inter_box[1] = cell_id(ix + 1 - px, iy + 1, iz,
                                                           ndx, ndy, ndz);
            ewald_sys->cell_sys[id].inter_box[2] = cell_id(ix     - px, iy + 1, iz,
                                                           ndx, ndy, ndz);
            if (ewald_sys->ndx == 2) {
                ewald_sys->cell_sys[id].inter_box[3] = - 1;
            } else {
                ewald_sys->cell_sys[id].inter_box[3] = cell_id(ix - 1 - px, iy + 1, iz,
                                                               ndx, ndy, ndz);
            }
            
            ewald_sys->cell_sys[id].inter_box[4] = cell_id(ix + 1,      iy    , iz - 1,
                                                           ndx, ndy, ndz);
            ewald_sys->cell_sys[id].inter_box[5] = cell_id(ix + 1 - px, iy + 1, iz - 1,
                                                           ndx, ndy, ndz);
            ewald_sys->cell_sys[id].inter_box[6] = cell_id(ix     - px, iy + 1, iz - 1,
                                                           ndx, ndy, ndz);
            if (ewald_sys->ndx == 2) {
                ewald_sys->cell_sys[id].inter_box[7] = - 1;
            } else {
                ewald_sys->cell_sys[id].inter_box[7] = cell_id(ix - 1 - px, iy + 1, iz - 1,
                                                               ndx, ndy, ndz);
            }
            
            if (ewald_sys->ndz == 2) {
                ewald_sys->cell_sys[id].inter_box[8] = - 1;
                ewald_sys->cell_sys[id].inter_box[9] = - 1;
                ewald_sys->cell_sys[id].inter_box[10] = - 1;
                ewald_sys->cell_sys[id].inter_box[11] = - 1;
            } else {
                ewald_sys->cell_sys[id].inter_box[8] = cell_id(ix + 1,      iy    , iz + 1,
                                                               ndx, ndy, ndz);
                ewald_sys->cell_sys[id].inter_box[9] = cell_id(ix + 1 - px, iy + 1, iz + 1,
                                                               ndx, ndy, ndz);
                ewald_sys->cell_sys[id].inter_box[10] = cell_id(ix     - px, iy + 1, iz + 1,
                                                                ndx, ndy, ndz);
                ewald_sys->cell_sys[id].inter_box[11] = cell_id(ix - 1 - px, iy + 1, iz + 1,
                                                                ndx, ndy, ndz);
            }
            
            ewald_sys->cell_sys[id].inter_box[12] = cell_id(ix         , iy    , iz + 1,
                                                            ndx, ndy, ndz);
            if (ewald_sys->ndx < 4) {
                ewald_sys->cell_sys[id].inter_box[13] = - 1;
                ewald_sys->cell_sys[id].inter_box[14] = - 1;
                ewald_sys->cell_sys[id].inter_box[15] = - 1;
            } else {
                ewald_sys->cell_sys[id].inter_box[13] = cell_id(ix - 2 - px, iy + 1, iz,
                                                                ndx, ndy, ndz);
                ewald_sys->cell_sys[id].inter_box[14] = cell_id(ix - 2 - px, iy + 1, iz - 1,
                                                                ndx, ndy, ndz);
                ewald_sys->cell_sys[id].inter_box[15] = cell_id(ix - 2 - px, iy + 1, iz + 1,
                                                                ndx, ndy, ndz);
            }
            
        } // for ix
        
    } // for iz
    
}

void rcell_list(struct ewald * ewald_sys,
                const int np,
                const double * pos)
{
    int ndx, ndy, ndz;
    int ix, iy, iz;
    int id;
    int i;
    
    ndx = ewald_sys->ndx;
    ndy = ewald_sys->ndy;
    ndz = ewald_sys->ndz;
    
    ewald_sys->head = calloc(ewald_sys->nd, sizeof(int));
    ewald_sys->list = calloc(np, sizeof(int));
    
    for (i = 0; i < ewald_sys->nd; i++) {
        ewald_sys->head[i] = -1;
    }
    
    for (i = 0; i < np; i++) {
        ix = (int) (pos[i * 3 + 0] / ewald_sys->rx);
        iy = (int) (pos[i * 3 + 1] / ewald_sys->ry);
        iz = (int) (pos[i * 3 + 2] / ewald_sys->rz);
        
        id = cell_id(ix, iy, iz,
                     ndx, ndy, ndz);
        ewald_sys->list[i] = ewald_sys->head[id];
        ewald_sys->head[id] = i;
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


void mob_pairs_list(struct ewald * ewald_sys,
                    const int np,
                    const double xlh, const double xrh,
                    const double ylh, const double yrh,
                    const double zlh, const double zrh,
                    const double lx, const double ly, const double lz,
                    const double rcut2,
                    const double strain,
                    double * pos)
{
    int i, j, k;
    int ib, jb;
    int ipair;
    int npair_ah, pair_acc;
    int indent;
    
    double s2, s;
    
    double dist[3];
    
    int * aux;
    double * aux_r;
    double * aux_e;
    
    aux = NULL;
    aux_r = NULL;
    aux_e = NULL;
    
    ewald_sys->tot_pairs = 0;
    ewald_sys->index = calloc(np * (np + 1) / 2, sizeof(int));
    
    for (i = 0; i < np; i++) {
        ewald_sys->pair_sys[i].bh_count = 0;
        ewald_sys->pair_sys[i].bh_list = malloc(ewald_sys->pair_sys[i].bh_count * sizeof(int));
        
        ewald_sys->pair_sys[i].ah_count = 0;
        ewald_sys->pair_sys[i].ah_list = malloc(ewald_sys->pair_sys[i].ah_count * sizeof(int));
        ewald_sys->pair_sys[i].r = malloc(ewald_sys->pair_sys[i].ah_count * sizeof(double));
        ewald_sys->pair_sys[i].e = malloc(ewald_sys->pair_sys[i].ah_count * 3 * sizeof(double));
    } // for i
    
    for (ib = 0; ib < ewald_sys->nd; ib++) {
        i = ewald_sys->head[ib];
        
        // loop over all particles in current box
        while (i >= 0) {
            
            // loop over all particles below i in the current box
            j = ewald_sys->list[i];
            
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
                
                if (s2 < rcut2) {
                    ewald_sys->tot_pairs += 1;
                    s = sqrt(s2);
                    
                    if (i > j) {
                        // pair afterhand
                        ewald_sys->pair_sys[j].ah_count += 1;
                        ewald_sys->pair_sys[j].ah_list = realloc(ewald_sys->pair_sys[j].ah_list, ewald_sys->pair_sys[j].ah_count * sizeof(int));
                        ewald_sys->pair_sys[j].r = realloc(ewald_sys->pair_sys[j].r, ewald_sys->pair_sys[j].ah_count * sizeof(double));
                        ewald_sys->pair_sys[j].e = realloc(ewald_sys->pair_sys[j].e, ewald_sys->pair_sys[j].ah_count * 3 * sizeof(double));
                        
                        ewald_sys->pair_sys[j].ah_list[ewald_sys->pair_sys[j].ah_count - 1] = i;
                        ewald_sys->pair_sys[j].r[ewald_sys->pair_sys[j].ah_count - 1] = s;
                        ewald_sys->pair_sys[j].e[(ewald_sys->pair_sys[j].ah_count - 1) * 3 + 0] = - dist[0] / s;
                        ewald_sys->pair_sys[j].e[(ewald_sys->pair_sys[j].ah_count - 1) * 3 + 1] = - dist[1] / s;
                        ewald_sys->pair_sys[j].e[(ewald_sys->pair_sys[j].ah_count - 1) * 3 + 2] = - dist[2] / s;
                        
                        // pair beforehand
                        ewald_sys->pair_sys[i].bh_count += 1;
                        ewald_sys->pair_sys[i].bh_list = realloc(ewald_sys->pair_sys[i].bh_list, ewald_sys->pair_sys[i].bh_count * sizeof(int));
                        
                        ewald_sys->pair_sys[i].bh_list[ewald_sys->pair_sys[i].bh_count - 1] = j;
                        
                    } else {
                        // pair afterhand
                        ewald_sys->pair_sys[i].ah_count += 1;
                        ewald_sys->pair_sys[i].ah_list = realloc(ewald_sys->pair_sys[i].ah_list, ewald_sys->pair_sys[i].ah_count * sizeof(int));
                        ewald_sys->pair_sys[i].r = realloc(ewald_sys->pair_sys[i].r, ewald_sys->pair_sys[i].ah_count * sizeof(double));
                        ewald_sys->pair_sys[i].e = realloc(ewald_sys->pair_sys[i].e, ewald_sys->pair_sys[i].ah_count * 3 * sizeof(double));
                        
                        ewald_sys->pair_sys[i].ah_list[ewald_sys->pair_sys[i].ah_count - 1] = j;
                        ewald_sys->pair_sys[i].r[ewald_sys->pair_sys[i].ah_count - 1] = s;
                        ewald_sys->pair_sys[i].e[(ewald_sys->pair_sys[i].ah_count - 1) * 3 + 0] = dist[0] / s;
                        ewald_sys->pair_sys[i].e[(ewald_sys->pair_sys[i].ah_count - 1) * 3 + 1] = dist[1] / s;
                        ewald_sys->pair_sys[i].e[(ewald_sys->pair_sys[i].ah_count - 1) * 3 + 2] = dist[2] / s;
                        
                        // pair beforehand
                        ewald_sys->pair_sys[j].bh_count += 1;
                        ewald_sys->pair_sys[j].bh_list = realloc(ewald_sys->pair_sys[j].bh_list, ewald_sys->pair_sys[j].bh_count * sizeof(int));
                        
                        ewald_sys->pair_sys[j].bh_list[ewald_sys->pair_sys[j].bh_count - 1] = i;
                    } // if i > j ...
                    
                } // if s2 < rcut2
                
                j = ewald_sys->list[j];
            } // while j >= 0
            
            // loop over neighbor boxes
            for (k = 0; k < ewald_sys->cell_sys[ib].tot_boxes; k++) {
                jb = ewald_sys->cell_sys[ib].inter_box[k];
                
                if (jb >= 0) {
                    // loop over all particles in neighbor boxes
                    j = ewald_sys->head[jb];
                    
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
                        
                        
                        if (s2 < rcut2) {
                            ewald_sys->tot_pairs += 1;
                            s = sqrt(s2);
                            if (i > j) {
                                // pair afterhand
                                ewald_sys->pair_sys[j].ah_count += 1;
                                ewald_sys->pair_sys[j].ah_list = realloc(ewald_sys->pair_sys[j].ah_list, ewald_sys->pair_sys[j].ah_count * sizeof(int));
                                ewald_sys->pair_sys[j].r = realloc(ewald_sys->pair_sys[j].r, ewald_sys->pair_sys[j].ah_count * sizeof(double));
                                ewald_sys->pair_sys[j].e = realloc(ewald_sys->pair_sys[j].e, ewald_sys->pair_sys[j].ah_count * 3 * sizeof(double));
                                
                                ewald_sys->pair_sys[j].ah_list[ewald_sys->pair_sys[j].ah_count - 1] = i;
                                ewald_sys->pair_sys[j].r[ewald_sys->pair_sys[j].ah_count - 1] = s;
                                ewald_sys->pair_sys[j].e[(ewald_sys->pair_sys[j].ah_count - 1) * 3 + 0] = - dist[0] / s;
                                ewald_sys->pair_sys[j].e[(ewald_sys->pair_sys[j].ah_count - 1) * 3 + 1] = - dist[1] / s;
                                ewald_sys->pair_sys[j].e[(ewald_sys->pair_sys[j].ah_count - 1) * 3 + 2] = - dist[2] / s;
                                
                                // pair beforehand
                                ewald_sys->pair_sys[i].bh_count += 1;
                                ewald_sys->pair_sys[i].bh_list = realloc(ewald_sys->pair_sys[i].bh_list, ewald_sys->pair_sys[i].bh_count * sizeof(int));
                                
                                ewald_sys->pair_sys[i].bh_list[ewald_sys->pair_sys[i].bh_count - 1] = j;
                                
                            } else {
                                // pair afterhand
                                ewald_sys->pair_sys[i].ah_count += 1;
                                ewald_sys->pair_sys[i].ah_list = realloc(ewald_sys->pair_sys[i].ah_list, ewald_sys->pair_sys[i].ah_count * sizeof(int));
                                ewald_sys->pair_sys[i].r = realloc(ewald_sys->pair_sys[i].r, ewald_sys->pair_sys[i].ah_count * sizeof(double));
                                ewald_sys->pair_sys[i].e = realloc(ewald_sys->pair_sys[i].e, ewald_sys->pair_sys[i].ah_count * 3 * sizeof(double));
                                
                                ewald_sys->pair_sys[i].ah_list[ewald_sys->pair_sys[i].ah_count - 1] = j;
                                ewald_sys->pair_sys[i].r[ewald_sys->pair_sys[i].ah_count - 1] = s;
                                ewald_sys->pair_sys[i].e[(ewald_sys->pair_sys[i].ah_count - 1) * 3 + 0] = dist[0] / s;
                                ewald_sys->pair_sys[i].e[(ewald_sys->pair_sys[i].ah_count - 1) * 3 + 1] = dist[1] / s;
                                ewald_sys->pair_sys[i].e[(ewald_sys->pair_sys[i].ah_count - 1) * 3 + 2] = dist[2] / s;
                                
                                // pair beforehand
                                ewald_sys->pair_sys[j].bh_count += 1;
                                ewald_sys->pair_sys[j].bh_list = realloc(ewald_sys->pair_sys[j].bh_list, ewald_sys->pair_sys[j].bh_count * sizeof(int));
                                
                                ewald_sys->pair_sys[j].bh_list[ewald_sys->pair_sys[j].bh_count - 1] = i;
                            } // if i > j ...
                            
                        } // if s2 < rcut2
                        
                        j = ewald_sys->list[j];
                    } // while j >= 0
                    
                } // jb >= 0
                
            } // for k
            
            i = ewald_sys->list[i];
        } // while i >= 0
        
    } // for ib
    
    /* Sort pairs */
    for (i = 0; i < np; i++) {
        if (ewald_sys->pair_sys[i].ah_count < 2) {
            continue;
        }
        
        aux = calloc(ewald_sys->pair_sys[i].ah_count, sizeof(int));
        aux_r = calloc(ewald_sys->pair_sys[i].ah_count, sizeof(double));
        aux_e = calloc(ewald_sys->pair_sys[i].ah_count * 3, sizeof(double));
        
        merge_sort_pair_all(0, ewald_sys->pair_sys[i].ah_count - 1,
                            ewald_sys->pair_sys[i].ah_list, ewald_sys->pair_sys[i].r, ewald_sys->pair_sys[i].e,
                            aux, aux_r, aux_e);
        
        free(aux);
        free(aux_r);
        free(aux_e);
    } // for i
    
    for (i = 0; i < np; i++) {
        if (ewald_sys->pair_sys[i].bh_count < 2) {
            continue;
        }
        
        aux = calloc(ewald_sys->pair_sys[i].bh_count, sizeof(int));
        
        merge_sort_pair(0, ewald_sys->pair_sys[i].bh_count - 1,
                        ewald_sys->pair_sys[i].bh_list,
                        aux);
        
        free(aux);
    } // for i
    
    /* Make list of interacting pairs */
    ewald_sys->row_id = malloc(ewald_sys->tot_pairs * sizeof(int));
    ewald_sys->col_id = malloc(ewald_sys->tot_pairs * sizeof(int));
    
    pair_acc = 0;
    for (i = 0; i < np; i++) {
        npair_ah = ewald_sys->pair_sys[i].ah_count;
        indent = i * (i + 1) / 2;
        
        ewald_sys->index[i * np + i - indent] = i;
        
        for (ipair = 0; ipair < npair_ah; ipair++) {
            j = ewald_sys->pair_sys[i].ah_list[ipair];
            
            ewald_sys->row_id[ipair + pair_acc] = i;
            ewald_sys->col_id[ipair + pair_acc] = j;
            
            ewald_sys->index[i * np + j - indent] = np + pair_acc + ipair;
        } // for ipair
        
        pair_acc += npair_ah;
    } // for i
    
}

/*
 Calculate real part of mobility functions
 */
void ewald_real_scalars(double xi, double r,
                        double * xa, double * ya,
                        double * yb,
                        double * xc, double * yc,
                        double * xg, double * yg,
                        double * yh,
                        double * xm, double * ym, double * zm)
{
    double r2, r3, r4, r5;
    double xir, xir2;
    double erfcxir;
    double expxir2;
    
    r2 = r * r;
    r3 = r2 * r;
    r4 = r2 * r2;
    r5 = r3 * r2;
    xir = xi * r;
    xir2 = xir * xir;
    erfcxir = erfc (xir);
    expxir2 = xi / sqrt (M_PI) * exp (- xir2);
    
    (* ya) = (0.75 / r + 0.5 / r3) * erfcxir + ((1.0 + xir2 * (14.0 + 4.0 * xir2 * (- 5.0 + xir2))) / r2 - 4.5 + 3.0 * xir2) * expxir2;
    (* xa) = (1.5 / r - 1.0 / r3) * erfcxir + (-2.0 + xir2 * (-12.0 - 4.0 * xir2)) / r2 * expxir2;
    (* yb) = - 0.75 / r2 * erfcxir - 1.5 * (+ 1.0 + xir2 * (- 6.0 + 2.0 * xir2)) / r * expxir2;
    (* yc) = - 0.375 / r3 * erfcxir - 0.75 * (+ 1.0 + xir2 * (+ 14.0 + xir2 * (-20.0 + 4.0 * xir2))) / r2 * expxir2;
    (* xc) = 0.75 / r3 * erfcxir - 0.75 * (- 2.0 + xir2 * (12.0 - 4.0 * xir2)) / r2 * expxir2;
    (* xg) = (2.25 / r2 - 3.6 / s4) * erfcxir + (- 1.5 * (- 3.0 + 6.0 * xir2) - 0.8 * (9.0 + xir2 * (+ 6.0 + xir2 * (- 48.0 + 12.0 * xir2))) / r2) / r * expxir2;
    (* yg) = 1.2 / r4 * erfcxir + (- 3.0 * ( xir2 * (2.0 - xir2)) - 0.8 * (- 3.0 + xir2 * (- 2.0 + xir2 * (- 26.0 + xir2 * (26.0 - 4.0 * xir2)))) / r2) / r * expxir2;
    (* yh) = - 9.0 / 8.0 / r3 * erfcxir + 1.5 * (- 1.5 + xir2 * (- 1.0 + xir2 * (8.0 - 2.0 * xir2))) / r2 * expxir2;
    (* xm) = (- 4.5 / r3 + 10.8 / r5) * erfcxir + (1.5 * (- 6.0 +  xir2 * (- 12.0 + 12.0 * xir2)) + 1.2 * (18.0 + xir2 * (12.0 + xir2 * (30.0 + xir2 * (- 66.0 + 12.0 * xir2 )))) / r2) / r2 * expxir2;
    (* ym) = (2.25 / r3 - 7.2 / r5) * erfcxir + (- 1.5 * (- 3.0 + xir2 * (6.0 + xir2 * (- 12.0 + 4.0 * xir2))) - 1.2 * (12.0 + xir2 * (8.0 + xir2 * (- 22.0 + xir2 * (58.0 + xir2 * (- 34.0 + 4.0 * xir2))))) / r2) / r2 * expxir2;
    (* zm) = 1.8 / r5 * erfcxir + (- 1.5 * (xir2 * (8.0 - 4.0 * xir2)) - 1.2 * (- 3.0 + xir2 * (- 2.0 + xir2 * (- 26.0 + xir2 * (26.0 - 4.0 * xir2)))) / r2) / r2 * expxir2;
}

/*
 Obtain tensor that couples UF
 */
void mob_abc(double xa, double ya,
             double yb,
             double xc, double yc,
             double * e,
             double * mabc)
{
    int i, j, k;
    double bdot;
    
    for (i = 0; i < 3; i++) {
        for (j = 0; j < 3; j++) {
            // a12, b12, bt12, c12
            mabc[(i + 0) * 12 + (j + 6)] = xa * e[i] * e[j] + ya * (kron[i * 3 + j] - e[i] * e[j]);
            bdot = 0.0;
            for (k = 0; k < 3; k++) {
                bdot += yb * levi[i * 9 + j * 3 + k] * e[k];
            }
            mabc[(i + 3) * 12 + (j + 6)] = bdot;
            mabc[(i + 0) * 12 + (j + 9)] = bdot;
            mabc[(i + 3) * 12 + (j + 9)] = xc * e[i] * e[j] + yc * (kron[i * 3 + j] - e[i] * e[j]);
            
            // a21, b21, bt21, c21
            mabc[(i + 6) * 12 + (j + 0)] = + mabc[(i + 0) * 12 + (j + 6)];
            mabc[(i + 9) * 12 + (j + 0)] = - mabc[(i + 3) * 12 + (j + 6)];
            mabc[(i + 6) * 12 + (j + 3)] = - mabc[(i + 0) * 12 + (j + 9)];
            mabc[(i + 9) * 12 + (j + 3)] = + mabc[(i + 3) * 12 + (j + 9)];
        }
    }
    
}

/*
 Obtain tensor that couples EF
 */
void mob_gh(double xg, double yg,
            double yh,
            double * e,
            double * mgh)
{
    int i, j, k, l;
    double hdot1, hdot2;
    double *g = calloc(27, sizeof(double));
    double *h = calloc(27, sizeof(double));
    
    for (i = 0; i < 3; i++) {
        for (j = 0; j < 3; j++) {
            for (k = 0; k < 3; k++) {
                g[i * 9 + j * 3 + k] = xg * (e[i] * e[j] - 1.0 / 3.0 * kron[i * 3 + j]) * e[k] + yg * (e[i] * kron[j * 3 + k] + e[j] * kron[i * 3 + k] - 2.0 * e[i] * e[j] * e[k]);
                hdot1 = 0.0;    hdot2 = 0.0;
                for (l = 0; l < 3; l++) {
                    hdot1 += yh * levi[j * 9 + k * 3 + l] * e[l];
                    hdot2 += yh * levi[i * 9 + k * 3 + l] * e[l];
                } // for l
                h[i * 9 + j * 3 + k] = hdot1 * e[i] + hdot2 * e[j];
                
            } // for k
        } // for j
    } // for i
    
    for (i = 0; i < 3; i++) {
        // g12
        mgh[0 * 12 + (i + 6)] = 2.0 * g[0 * 9 + 0 * 3 + i] + g[1 * 9 + 1 * 3 + i];
        mgh[1 * 12 + (i + 6)] = 2.0 * g[0 * 9 + 1 * 3 + i];
        mgh[2 * 12 + (i + 6)] = 2.0 * g[0 * 9 + 2 * 3 + i];
        mgh[3 * 12 + (i + 6)] = 2.0 * g[1 * 9 + 2 * 3 + i];
        mgh[4 * 12 + (i + 6)] = g[0 * 9 + 0 * 3 + i] + 2.0 * g[1 * 9 + 1 * 3 + i];
        // h12
        mgh[0 * 12 + (i + 9)] = 2.0 * h[0 * 9 + 0 * 3 + i] + h[1 * 9 + 1 * 3 + i];
        mgh[1 * 12 + (i + 9)] = 2.0 * h[0 * 9 + 1 * 3 + i];
        mgh[2 * 12 + (i + 9)] = 2.0 * h[0 * 9 + 2 * 3 + i];
        mgh[3 * 12 + (i + 9)] = 2.0 * h[1 * 9 + 2 * 3 + i];
        mgh[4 * 12 + (i + 9)] = h[0 * 9 + 0 * 3 + i] + 2.0 * h[1 * 9 + 1 * 3 + i];
        // g21 = -g12
        mgh[5 * 12 + (i + 0)] = -mgh[0 * 12 + (i + 6)];
        mgh[6 * 12 + (i + 0)] = -mgh[1 * 12 + (i + 6)];
        mgh[7 * 12 + (i + 0)] = -mgh[2 * 12 + (i + 6)];
        mgh[8 * 12 + (i + 0)] = -mgh[3 * 12 + (i + 6)];
        mgh[9 * 12 + (i + 0)] = -mgh[4 * 12 + (i + 6)];
        // h21 = h12
        mgh[5 * 12 + (i + 3)] = +mgh[0 * 12 + (i + 9)];
        mgh[6 * 12 + (i + 3)] = +mgh[1 * 12 + (i + 9)];
        mgh[7 * 12 + (i + 3)] = +mgh[2 * 12 + (i + 9)];
        mgh[8 * 12 + (i + 3)] = +mgh[3 * 12 + (i + 9)];
        mgh[9 * 12 + (i + 3)] = +mgh[4 * 12 + (i + 9)];
    }
    
    free(g);
    free(h);
}

/*
 Obtain tensor that couples ES
 */
void mob_zm(double xm, double ym, double zm,
            double * e,
            double * mzm)
{
    int i, j, k, l;
    double *m = calloc(81, sizeof(double));
    
    for (i = 0; i < 3; i++) {
        for (j = 0; j < 3; j++) {
            for (k = 0; k < 3; k++) {
                for (l = 0; l < 3; l++) {
                    m[i * 27 + j * 9 + k * 3 + l] = 1.5 * xm * (e[i] * e[j] - 1.0 / 3.0 * kron[i * 3 + j]) * (e[k] * e[l] - 1.0 / 3.0 * kron[k * 3 + l])
                    + 0.5 * ym * (e[i] * kron[j * 3 + l] * e[k] + e[j] * kron[i * 3 + l] * e[k] + e[i] * kron[j * 3 + k] * e[l] + e[j] * kron[i * 3 + k] * e[l] - 4.0 * e[i] * e[j] * e[k] * e[l])
                    + 0.5 * zm * (kron[i * 3 + k] * kron[j * 3 + l] + kron[j * 3 + k] * kron[i * 3 + l] - kron[i * 3 + j] * kron[k * 3 + l] + e[i] * e[j] * kron[k * 3 + l] + kron[i * 3 + j] * e[k] * e[l] + e[i] * e[j] * e[k] * e[l] - e[i] * kron[j * 3 + l] * e[k] - e[j] * kron[i * 3 + l] * e[k] - e[i] * kron[j * 3 + k] * e[l] - e[j] * kron[i * 3 + k] * e[l]);
                } // for l
            } // for k
        } // for j
    } // for i
    
    // m12
    mzm[0 * 10 + 5] = m[0 * 27 + 0 * 9 + 0 * 3 + 0] - 2.0 * m[0 * 27 + 0 * 9 + 2 * 3 + 2] + m[2 * 27 + 2 * 9 + 2 * 3 + 2];
    mzm[0 * 10 + 6] = 2.0 * (m[0 * 27 + 0 * 9 + 0 * 3 + 1] - m[2 * 27 + 2 * 9 + 0 * 3 + 1]);
    mzm[0 * 10 + 7] = 2.0 * (m[0 * 27 + 0 * 9 + 0 * 3 + 2] - m[2 * 27 + 2 * 9 + 0 * 3 + 2]);
    mzm[0 * 10 + 8] = 2.0 * (m[0 * 27 + 0 * 9 + 1 * 3 + 2] - m[2 * 27 + 2 * 9 + 1 * 3 + 2]);
    mzm[0 * 10 + 9] = m[0 * 27 + 0 * 9 + 1 * 3 + 1] - m[0 * 27 + 0 * 9 + 2 * 3 + 2] - m[2 * 27 + 2 * 9 + 1 * 3 + 1] + m[2 * 27 + 2 * 9 + 2 * 3 + 2];
    
    mzm[1 * 10 + 5] = 2.0 * (m[0 * 27 + 1 * 9 + 0 * 3 + 0] - m[0 * 27 + 1 * 9 + 2 * 3 + 2]);
    mzm[1 * 10 + 6] = 4.0 * m[0 * 27 + 1 * 9 + 0 * 3 + 1];
    mzm[1 * 10 + 7] = 4.0 * m[0 * 27 + 1 * 9 + 0 * 3 + 2];
    mzm[1 * 10 + 8] = 4.0 * m[0 * 27 + 1 * 9 + 1 * 3 + 2];
    mzm[1 * 10 + 9] = 2.0 * (m[0 * 27 + 1 * 9 + 1 * 3 + 1] - m[0 * 27 + 1 * 9 + 2 * 3 + 2]);
    
    mzm[2 * 10 + 5] = 2.0 * (m[0 * 27 + 2 * 9 + 0 * 3 + 0] - m[0 * 27 + 2 * 9 + 2 * 3 + 2]);
    mzm[2 * 10 + 6] = 4.0 * m[0 * 27 + 2 * 9 + 0 * 3 + 1];
    mzm[2 * 10 + 7] = 4.0 * m[0 * 27 + 2 * 9 + 0 * 3 + 2];
    mzm[2 * 10 + 8] = 4.0 * m[0 * 27 + 2 * 9 + 1 * 3 + 2];
    mzm[2 * 10 + 9] = 2.0 * (m[0 * 27 + 2 * 9 + 1 * 3 + 1] - m[0 * 27 + 2 * 9 + 2 * 3 + 2]);
    
    mzm[3 * 10 + 5] = 2.0 * (m[1 * 27 + 2 * 9 + 0 * 3 + 0] - m[1 * 27 + 2 * 9 + 2 * 3 + 2]);
    mzm[3 * 10 + 6] = 4.0 * m[1 * 27 + 2 * 9 + 0 * 3 + 1];
    mzm[3 * 10 + 7] = 4.0 * m[1 * 27 + 2 * 9 + 0 * 3 + 2];
    mzm[3 * 10 + 8] = 4.0 * m[1 * 27 + 2 * 9 + 1 * 3 + 2];
    mzm[3 * 10 + 9] = 2.0 * (m[1 * 27 + 2 * 9 + 1 * 3 + 1] - m[1 * 27 + 2 * 9 + 2 * 3 + 2]);
    
    mzm[4 * 10 + 5] = m[1 * 27 + 1 * 9 + 0 * 3 + 0] - m[1 * 27 + 1 * 9 + 2 * 3 + 2] - m[0 * 27 + 0 * 9 + 2 * 3 + 2] + m[2 * 27 + 2 * 9 + 2 * 3 + 2];
    mzm[4 * 10 + 6] = 2.0 * (m[1 * 27 + 1 * 9 + 0 * 3 + 1] - m[2 * 27 + 2 * 9 + 0 * 3 + 1]);
    mzm[4 * 10 + 7] = 2.0 * (m[1 * 27 + 1 * 9 + 0 * 3 + 2] - m[2 * 27 + 2 * 9 + 0 * 3 + 2]);
    mzm[4 * 10 + 8] = 2.0 * (m[1 * 27 + 1 * 9 + 1 * 3 + 2] - m[2 * 27 + 2 * 9 + 1 * 3 + 2]);
    mzm[4 * 10 + 9] = m[1 * 27 + 1 * 9 + 1 * 3 + 1] - 2.0 * m[1 * 27 + 1 * 9 + 2 * 3 + 2] + m[2 * 27 + 2 * 9 + 2 * 3 + 2];
    // m21 = m12
    mzm[5 * 10 + 0] = mzm[0 * 10 + 5];
    mzm[5 * 10 + 1] = mzm[0 * 10 + 6];
    mzm[5 * 10 + 2] = mzm[0 * 10 + 7];
    mzm[5 * 10 + 3] = mzm[0 * 10 + 8];
    mzm[5 * 10 + 4] = mzm[0 * 10 + 9];
    
    mzm[6 * 10 + 0] = mzm[1 * 10 + 5];
    mzm[6 * 10 + 1] = mzm[1 * 10 + 6];
    mzm[6 * 10 + 2] = mzm[1 * 10 + 7];
    mzm[6 * 10 + 3] = mzm[1 * 10 + 8];
    mzm[6 * 10 + 4] = mzm[1 * 10 + 9];
    
    mzm[7 * 10 + 0] = mzm[2 * 10 + 5];
    mzm[7 * 10 + 1] = mzm[2 * 10 + 6];
    mzm[7 * 10 + 2] = mzm[2 * 10 + 7];
    mzm[7 * 10 + 3] = mzm[2 * 10 + 8];
    mzm[7 * 10 + 4] = mzm[2 * 10 + 9];
    
    mzm[8 * 10 + 0] = mzm[3 * 10 + 5];
    mzm[8 * 10 + 1] = mzm[3 * 10 + 6];
    mzm[8 * 10 + 2] = mzm[3 * 10 + 7];
    mzm[8 * 10 + 3] = mzm[3 * 10 + 8];
    mzm[8 * 10 + 4] = mzm[3 * 10 + 9];
    
    mzm[9 * 10 + 0] = mzm[4 * 10 + 5];
    mzm[9 * 10 + 1] = mzm[4 * 10 + 6];
    mzm[9 * 10 + 2] = mzm[4 * 10 + 7];
    mzm[9 * 10 + 3] = mzm[4 * 10 + 8];
    mzm[9 * 10 + 4] = mzm[4 * 10 + 9];
    
    free(m);
}

/*
 Obtain real part of mobility tensors
 Input:
 np                  number of particles
 xi                  xi parameter
 */
void ewald_real_mat(struct ewald * ewald_sys,
                    const int np,
                    const double xi)
{
    int n5, n6;
    int ipair;
    int in, jn;
    int npair_ah, pair_acc;
    int pai;
    int i, j;
    int jc;
    
    double xa, ya;
    double yb;
    double xc, yc;
    double xg, yg;
    double yh;
    double xm, ym, zm;
    
    double mabc[144];
    double mgh[120];
    double mzm[100];
    
    n5 = ewald_sys->tot_pairs * 5;
    n6 = ewald_sys->tot_pairs * 6;
    
    ewald_sys->muf = calloc(6 * n6, sizeof(double));
    ewald_sys->mus = calloc(6 * n5, sizeof(double));
    ewald_sys->mef = calloc(5 * n6, sizeof(double));
    ewald_sys->mes = calloc(5 * n5, sizeof(double));
    
    pair_acc = 0;
    
    for (in = 0; in < np; in++) {
        npair_ah = ewald_sys->pair_sys[in].ah_count;
        
        for (ipair = 0; ipair < npair_ah; ipair++) {
            jn = ewald_sys->pair_sys[in].ah_list[ipair];
            
            ewald_real_scalars(xi, ewald_sys->pair_sys[in].r[ipair],
                               & xa, & ya,
                               & yb,
                               & xc, & yc,
                               & xg, & yg,
                               & yh,
                               & xm, & ym, & zm);
            
            mob_abc(xa, ya,
                    yb,
                    xc, yc,
                    ewald_sys->pair_sys[in].e + ipair * 3,
                    mabc);
            mob_gh(xg, yg,
                   yh,
                   ewald_sys->pair_sys[in].e + ipair * 3,
                   mgh);
            mob_zm(xm, ym, zm,
                   ewald_sys->pair_sys[in].e + ipair * 3,
                   mzm);
            
            pai = pair_acc + ipair;
            
            // Muf + Mus
            for (i = 0; i < 6; i++) {
                // Muf
                for (j = 0; j < 6; j++) {
                    jc = pai * 6 + j;
                    
                    ewald_sys->muf[i * n6 + jc] += mabc[(i + 0) * 12 + (j + 6)];
                }
                // Mus
                for (j = 0; j < 5; j++) {
                    jc = pai * 5 + j;
                    
                    ewald_sys->mus[i * n5 + jc] += mgh[(j + 5) * 12 + (i + 0)];
                }
            } // for i
            
            // Mef + Mes
            for (i = 0; i < 5; i++) {
                // Mef
                for (j = 0; j < 6; j++) {
                    jc = pai * 6 + j;
                    
                    ewald_sys->mef[i * n6 + jc] += mgh[(i + 0) * 12 + (j + 6)];
                }
                // Mes
                for (j = 0; j < 5; j++) {
                    jc = pai * 5 + j;
                    
                    ewald_sys->mes[i * n5 + jc] += mzm[(i + 0) * 10 + (j + 5)];
                }
            } // for i
            
        } // for ipair
        
        pair_acc += npair_ah;
        
    } // for in
    
}

/*
 Spreads forces/torques/stresslet onto p**3 grid nodes that are closest to particles
 Input:
 mode                                            mode = 0, FT version. mode = 1, FTS version
 np                                              number of particles
 p3                                              p3 = p**3, p = number of node in one direction
 gaussCo, gaussExp                               coeffcients of Gaussian kernel
 f                                               force/torque
 s                                               stresslet (sxx, sxy, sxz, syz, syy)
 Output:
 gridFx, gridFy, gridFz                          x,y,z-components of force of grid node
 gridTx, gridTy, gridTz                          x,y,z-components of torque of grid node
 gridSxx, gridSxy, gridSxz, gridSyz, gridSyy     stresslet of grid node
 */
void ewald_recip_spread(struct ewald * ewald_sys,
                        const int np, const int p3,
                        const double gaussCo, const double gaussExp,
                        double * f, double * s,
                        fftw_complex * gridFx, fftw_complex * gridFy, fftw_complex * gridFz,
                        fftw_complex * gridTx, fftw_complex * gridTy, fftw_complex * gridTz,
                        fftw_complex * gridSxx, fftw_complex * gridSxy, fftw_complex * gridSxz, fftw_complex * gridSyz, fftw_complex * gridSyy)
{
    int i, j;
    int id;
    
    double fx, fy, fz;
    double tx, ty, tz;
    double sxx, sxy, sxz, syz, syy;
    double r2;
    double coeff;
    
    for (i = 0; i < np; i++) {
        fx = f[i * 6 + 0];
        fy = f[i * 6 + 1];
        fz = f[i * 6 + 2];
        
        tx = f[i * 6 + 3];
        ty = f[i * 6 + 4];
        tz = f[i * 6 + 5];
        
        sxx = s[i * 5 + 0];
        sxy = s[i * 5 + 1];
        sxz = s[i * 5 + 2];
        syz = s[i * 5 + 3];
        syy = s[i * 5 + 4];
        
        for (j = 0; j < p3; j++) {
            id = ewald_sys->gridp_sys[i].ids[j];
            r2 = ewald_sys->gridp_sys[i].r2[j];
            
            coeff = gaussCo * exp( - gaussExp * r2 );
            
            gridFx[id][0] += fx * coeff;
            gridFy[id][0] += fy * coeff;
            gridFz[id][0] += fz * coeff;
            
            gridTx[id][0] += tx * coeff;
            gridTy[id][0] += ty * coeff;
            gridTz[id][0] += tz * coeff;
            
            gridSxx[id][0] += sxx * coeff;
            gridSxy[id][0] += sxy * coeff;
            gridSxz[id][0] += sxz * coeff;
            gridSyz[id][0] += syz * coeff;
            gridSyy[id][0] += syy * coeff;
            
        } // for j
        
    } // for i
    
}

/*
 Calculate scaled force/torque/stresslet of grid by multiplying Green's function
 Input:
 mode                                            mode = 0, FT version. mode = 1, FTS version
 nk                                              number of grid points
 Input/Output
 gridFx, gridFy, gridFz                          scaled x,y,z-components of force of grid node
 gridTx, gridTy, gridTz                          scaled x,y,z-components of torque of grid node
 gridSxx, gridSxy, gridSxz, gridSyz, gridSyy     scaled stresslet of grid node
 */
void ewald_recip_green(struct ewald * ewald_sys,
                       const int nk,
                       fftw_complex * gridFx, fftw_complex * gridFy, fftw_complex * gridFz,
                       fftw_complex * gridTx, fftw_complex * gridTy, fftw_complex * gridTz,
                       fftw_complex * gridSxx, fftw_complex * gridSxy, fftw_complex * gridSxz, fftw_complex * gridSyz, fftw_complex * gridSyy)
{
    int i;
    
    double kx, ky, kz;
    double kxx, kxy, kxz, kyy, kyz;
    
    fftw_complex fx, fy, fz;
    fftw_complex tx, ty, tz;
    fftw_complex fk, tk;
    fftw_complex fxx, fxy, fxz;
    fftw_complex fyx, fyy, fyz;
    fftw_complex fzx, fzy;
    fftw_complex sxx, sxy, sxz, syz, syy, szz;
    fftw_complex skx, sky, skz;
    fftw_complex ksk;
    
    for (i = 0; i < nk; i++) {
        kx = ewald_sys->gridk_sys[i].kx;
        ky = ewald_sys->gridk_sys[i].ky;
        kz = ewald_sys->gridk_sys[i].kz;
        
        kxx = kx * kx;
        kxy = kx * ky;
        kxz = kx * kz;
        kyy = ky * ky;
        kyz = ky * kz;
        
        fx[0] = gridFx[i][0];       fx[1] = gridFx[i][1];
        fy[0] = gridFy[i][0];       fy[1] = gridFy[i][1];
        fz[0] = gridFz[i][0];       fz[1] = gridFz[i][1];
        
        tx[0] = gridTx[i][0];       tx[1] = gridTx[i][1];
        ty[0] = gridTy[i][0];       ty[1] = gridTy[i][1];
        tz[0] = gridTz[i][0];       tz[1] = gridTz[i][1];
        
        sxx[0] = gridSxx[i][0];         sxx[1] = gridSxx[i][1];
        sxy[0] = gridSxy[i][0];         sxy[1] = gridSxy[i][1];
        sxz[0] = gridSxz[i][0];         sxz[1] = gridSxz[i][1];
        syz[0] = gridSyz[i][0];         syz[1] = gridSyz[i][1];
        syy[0] = gridSyy[i][0];         syy[1] = gridSyy[i][1];
        szz[0] = - sxx[0] - syy[0];     szz[1] = - sxx[1] - syy[1];
        
        /* F . k */
        fk[0] = fx[0] * kx + fy[0] * ky + fz[0] * kz;
        fk[1] = fx[1] * kx + fy[1] * ky + fz[1] * kz;
        
        /* F k */
        fxx[0] = fx[0] * kx;   fxx[1] = fx[1] * kx;
        fxy[0] = fx[0] * ky;   fxy[1] = fx[1] * ky;
        fxz[0] = fx[0] * kz;   fxz[1] = fx[1] * kz;
        
        fyx[0] = fy[0] * kx;   fyx[1] = fy[1] * kx;
        fyy[0] = fy[0] * ky;   fyy[1] = fy[1] * ky;
        fyz[0] = fy[0] * kz;   fyz[1] = fy[1] * kz;
        
        fzx[0] = fz[0] * kx;   fzx[1] = fz[1] * kx;
        fzy[0] = fz[0] * ky;   fzy[1] = fz[1] * ky;
        
        /* T . k */
        tk[0] = tx[0] * kx + ty[0] * ky + tz[0] * kz;
        tk[1] = tx[1] * kx + ty[1] * ky + tz[1] * kz;
        
        /*
         S . k
         */
        skx[0] = sxx[0] * kx + sxy[0] * ky + sxz[0] * kz;       skx[1] = sxx[1] * kx + sxy[1] * ky + sxz[1] * kz;
        sky[0] = sxy[0] * kx + syy[0] * ky + syz[0] * kz;       sky[1] = sxy[1] * kx + syy[1] * ky + syz[1] * kz;
        skz[0] = sxz[0] * kx + syz[0] * ky + szz[0] * kz;       skz[1] = sxz[1] * kx + syz[1] * ky + szz[1] * kz;
        
        /*
         k . S . k
         */
        ksk[0] = kx * skx[0] + ky * sky[0] + kz * skz[0];       ksk[1] = kx * skx[1] + ky * sky[1] + kz * skz[1];
        
        //  UF coupling
        gridFx[i][0] = ewald_sys->gridk_sys[i].ya * (fx[0] - fk[0] * kx);
        gridFy[i][0] = ewald_sys->gridk_sys[i].ya * (fy[0] - fk[0] * ky);
        gridFz[i][0] = ewald_sys->gridk_sys[i].ya * (fz[0] - fk[0] * kz);
        
        gridFx[i][1] = ewald_sys->gridk_sys[i].ya * (fx[1] - fk[1] * kx);
        gridFy[i][1] = ewald_sys->gridk_sys[i].ya * (fy[1] - fk[1] * ky);
        gridFz[i][1] = ewald_sys->gridk_sys[i].ya * (fz[1] - fk[1] * kz);
        
        // UT coupling
        gridFx[i][0] += + ewald_sys->gridk_sys[i].yb * (ty[1] * kz - tz[1] * ky);
        gridFy[i][0] += + ewald_sys->gridk_sys[i].yb * (tz[1] * kx - tx[1] * kz);
        gridFz[i][0] += + ewald_sys->gridk_sys[i].yb * (tx[1] * ky - ty[1] * kx);
        
        gridFx[i][1] += - ewald_sys->gridk_sys[i].yb * (ty[0] * kz - tz[0] * ky);
        gridFy[i][1] += - ewald_sys->gridk_sys[i].yb * (tz[0] * kx - tx[0] * kz);
        gridFz[i][1] += - ewald_sys->gridk_sys[i].yb * (tx[0] * ky - ty[0] * kx);
        
        // US coupling
        gridFx[i][0] += + ewald_sys->gridk_sys[i].yg * (skx[1] - kx * ksk[1]);
        gridFy[i][0] += + ewald_sys->gridk_sys[i].yg * (sky[1] - ky * ksk[1]);
        gridFz[i][0] += + ewald_sys->gridk_sys[i].yg * (skz[1] - kz * ksk[1]);
        
        gridFx[i][1] += - ewald_sys->gridk_sys[i].yg * (skx[0] - kx * ksk[0]);
        gridFy[i][1] += - ewald_sys->gridk_sys[i].yg * (sky[0] - ky * ksk[0]);
        gridFz[i][1] += - ewald_sys->gridk_sys[i].yg * (skz[0] - kz * ksk[0]);
        
        // OF coupling
        gridTx[i][0] = + ewald_sys->gridk_sys[i].yb * (fyz[1] - fzy[1]);
        gridTy[i][0] = + ewald_sys->gridk_sys[i].yb * (fzx[1] - fxz[1]);
        gridTz[i][0] = + ewald_sys->gridk_sys[i].yb * (fxy[1] - fyx[1]);
        
        gridTx[i][1] = - ewald_sys->gridk_sys[i].yb * (fyz[0] - fzy[0]);
        gridTy[i][1] = - ewald_sys->gridk_sys[i].yb * (fzx[0] - fxz[0]);
        gridTz[i][1] = - ewald_sys->gridk_sys[i].yb * (fxy[0] - fyx[0]);
        
        // OT coupling
        gridTx[i][0] += + ewald_sys->gridk_sys[i].yc * (tx[0] - kx * tk[0]);
        gridTy[i][0] += + ewald_sys->gridk_sys[i].yc * (ty[0] - ky * tk[0]);
        gridTz[i][0] += + ewald_sys->gridk_sys[i].yc * (tz[0] - kz * tk[0]);
        
        gridTx[i][1] += + ewald_sys->gridk_sys[i].yc * (tx[1] - kx * tk[1]);
        gridTy[i][1] += + ewald_sys->gridk_sys[i].yc * (ty[1] - ky * tk[1]);
        gridTz[i][1] += + ewald_sys->gridk_sys[i].yc * (tz[1] - kz * tk[1]);
        
        // OS coupling
        gridTx[i][0] += + ewald_sys->gridk_sys[i].yh * 2.0 * (skz[0] * ky - sky[0] * kz);
        gridTy[i][0] += + ewald_sys->gridk_sys[i].yh * 2.0 * (skx[0] * kz - skz[0] * kx);
        gridTz[i][0] += + ewald_sys->gridk_sys[i].yh * 2.0 * (sky[0] * kx - skx[0] * ky);
        
        gridTx[i][1] += + ewald_sys->gridk_sys[i].yh * 2.0 * (skz[1] * ky - sky[1] * kz);
        gridTy[i][1] += + ewald_sys->gridk_sys[i].yh * 2.0 * (skx[1] * kz - skz[1] * kx);
        gridTz[i][1] += + ewald_sys->gridk_sys[i].yh * 2.0 * (sky[1] * kx - skx[1] * ky);
        
        /* Calculate E */
        // EF coupling
        gridSxx[i][0] = - ewald_sys->gridk_sys[i].yg * (fxx[1] - kxx * fk[1]);
        gridSxy[i][0] = - ewald_sys->gridk_sys[i].yg * ( 0.5 * (fxy[1] + fyx[1]) - kxy * fk[1]);
        gridSxz[i][0] = - ewald_sys->gridk_sys[i].yg * ( 0.5 * (fxz[1] + fzx[1]) - kxz * fk[1]);
        gridSyz[i][0] = - ewald_sys->gridk_sys[i].yg * ( 0.5 * (fyz[1] + fzy[1]) - kyz * fk[1]);
        gridSyy[i][0] = - ewald_sys->gridk_sys[i].yg * (fyy[1] - kyy * fk[1]);
        
        gridSxx[i][1] = + ewald_sys->gridk_sys[i].yg * (fxx[0] - kxx * fk[0]);
        gridSxy[i][1] = + ewald_sys->gridk_sys[i].yg * ( 0.5 * (fxy[0] + fyx[0]) - kxy * fk[0]);
        gridSxz[i][1] = + ewald_sys->gridk_sys[i].yg * ( 0.5 * (fxz[0] + fzx[0]) - kxz * fk[0]);
        gridSyz[i][1] = + ewald_sys->gridk_sys[i].yg * ( 0.5 * (fyz[0] + fzy[0]) - kyz * fk[0]);
        gridSyy[i][1] = + ewald_sys->gridk_sys[i].yg * (fyy[0] - kyy * fk[0]);
        
        // ET coupling
        gridSxx[i][0] += - ewald_sys->gridk_sys[i].yh * ( 2.0 * kx * (ky * tz[0] - kz * ty[0]) );
        gridSxy[i][0] += - ewald_sys->gridk_sys[i].yh * ( ky * (ky * tz[0] - kz * ty[0]) + kx * (kz * tx[0] - kx * tz[0]) );
        gridSxz[i][0] += - ewald_sys->gridk_sys[i].yh * ( kz * (ky * tz[0] - kz * ty[0]) + kx * (kx * ty[0] - ky * tx[0]) );
        gridSyz[i][0] += - ewald_sys->gridk_sys[i].yh * ( kz * (kz * tx[0] - kx * tz[0]) + ky * (kx * ty[0] - ky * tx[0]) );
        gridSyy[i][0] += - ewald_sys->gridk_sys[i].yh * ( 2.0 * ky * (kz * tx[0] - kx * tz[0]) );
        
        gridSxx[i][1] += - ewald_sys->gridk_sys[i].yh * ( 2.0 * kx * (ky * tz[1] - kz * ty[1]) );
        gridSxy[i][1] += - ewald_sys->gridk_sys[i].yh * ( ky * (ky * tz[1] - kz * ty[1]) + kx * (kz * tx[1] - kx * tz[1]) );
        gridSxz[i][1] += - ewald_sys->gridk_sys[i].yh * ( kz * (ky * tz[1] - kz * ty[1]) + kx * (kx * ty[1] - ky * tx[1]) );
        gridSyz[i][1] += - ewald_sys->gridk_sys[i].yh * ( kz * (kz * tx[1] - kx * tz[1]) + ky * (kx * ty[1] - ky * tx[1]) );
        gridSyy[i][1] += - ewald_sys->gridk_sys[i].yh * ( 2.0 * ky * (kz * tx[1] - kx * tz[1]) );
        
        // ES coupling
        gridSxx[i][0] += + ewald_sys->gridk_sys[i].ym * (2.0 * kx * skx[0] - 2.0 * kxx * ksk[0]);
        gridSxy[i][0] += + ewald_sys->gridk_sys[i].ym * (ky * skx[0] + kx * sky[0] - 2.0 * kxy * ksk[0]);
        gridSxz[i][0] += + ewald_sys->gridk_sys[i].ym * (kz * skx[0] + kx * skz[0] - 2.0 * kxz * ksk[0]);
        gridSyz[i][0] += + ewald_sys->gridk_sys[i].ym * (kz * sky[0] + ky * skz[0] - 2.0 * kyz * ksk[0]);
        gridSyy[i][0] += + ewald_sys->gridk_sys[i].ym * (2.0 * ky * sky[0] - 2.0 * kyy * ksk[0]);
        
        gridSxx[i][1] += + ewald_sys->gridk_sys[i].ym * (2.0 * kx * skx[1] - 2.0 * kxx * ksk[1]);
        gridSxy[i][1] += + ewald_sys->gridk_sys[i].ym * (ky * skx[1] + kx * sky[1] - 2.0 * kxy * ksk[1]);
        gridSxz[i][1] += + ewald_sys->gridk_sys[i].ym * (kz * skx[1] + kx * skz[1] - 2.0 * kxz * ksk[1]);
        gridSyz[i][1] += + ewald_sys->gridk_sys[i].ym * (kz * sky[1] + ky * skz[1] - 2.0 * kyz * ksk[1]);
        gridSyy[i][1] += + ewald_sys->gridk_sys[i].ym * (2.0 * ky * sky[1] - 2.0 * kyy * ksk[1]);
        //
    } // for i
    
}

/*
 Interpolate velocity/velocity-gradient of grid points back to particles
 Input:
 np                                              number of particles
 p3                                              p3 = p**3, p = number of node in one direction
 gaussCo, gaussExp                               coeffcients of Gaussian kernel
 vk                                              volume of one grid cell
 gridFx, gridFy, gridFz                          scaled x,y,z-components of force of grid node
 gridTx, gridTy, gridTz                          scaled x,y,z-components of torque of grid node
 gridSxx, gridSxy, gridSxz, gridSyz, gridSyy     scaled stresslet of grid node
 Output:
 u                                               translational/rotational velocity of particles
 e                                               rate-of-strain of particle (Exx, Exy, Exz, Eyz, Eyy)
 */
void ewald_recip_interp(struct ewald * ewald_sys,
                        const int np, const int p3,
                        const double gaussCo, const double gaussExp,
                        const double vk,
                        fftw_complex * gridFx, fftw_complex * gridFy, fftw_complex * gridFz,
                        fftw_complex * gridTx, fftw_complex * gridTy, fftw_complex * gridTz,
                        fftw_complex * gridSxx, fftw_complex * gridSxy, fftw_complex * gridSxz, fftw_complex * gridSyz, fftw_complex * gridSyy,
                        double * u, double * e)
{
    int i, j;
    int id;
    
    double r2;
    double coeff;
    
    for (i = 0; i < np; i++) {
        
        for (j = 0; j < p3; j++) {
            id = ewald_sys->gridp_sys[i].ids[j];
            r2 = ewald_sys->gridp_sys[i].r2[j];
            
            coeff = vk * gaussCo * exp( - gaussExp * r2 );
            
            u[i * 6 + 0] += coeff * gridFx[id][0];
            u[i * 6 + 1] += coeff * gridFy[id][0];
            u[i * 6 + 2] += coeff * gridFz[id][0];
            
            u[i * 6 + 3] += coeff * gridTx[id][0];
            u[i * 6 + 4] += coeff * gridTy[id][0];
            u[i * 6 + 5] += coeff * gridTz[id][0];
            
            e[i * 5 + 0] += coeff * gridSxx[id][0];
            e[i * 5 + 1] += coeff * gridSxy[id][0];
            e[i * 5 + 2] += coeff * gridSxz[id][0];
            e[i * 5 + 3] += coeff * gridSyz[id][0];
            e[i * 5 + 4] += coeff * gridSyy[id][0];
        }
        
    } // for i
    
}

/*
 Wrap up all the functions that calculate the reciprocal contribution to velocity/velocity-gradient with input of force/torque/stresslet
 Input:
 mode                                            mode = 0, FT version. mode = 1, FTS version
 np                                              number of particles
 p3                                              p3 = p**3, p = number of node in one direction
 nkx, nky, nkz                                   number of grid points in x,y,z-directions
 nk                                              total number of grid points
 gaussCo, gaussExp                               coeffcients of Gaussian kernel
 vk                                              volume of one grid cell
 f                                               force/torque of particles
 s                                               stresslet of particles
 Output:
 u                                               translational/rotational velocity of particles
 e                                               rate-of-strain of particles (Exx-Ezz, 2Exy, 2Exz, 2Eyz, Eyy-Ezz)
 */
void ewald_recip_wrap(struct ewald * ewald_sys,
                      const int np, const int p3,
                      const int nkx, const int nky, const int nkz,
                      const int nk,
                      const double gaussCo, const double gaussExp,
                      const double vk,
                      double * f, double * s,
                      double * u, double * e)
{
    int i;
    
    double * etemp;
    
    fftw_complex * gridFx;
    fftw_complex * gridFy;
    fftw_complex * gridFz;
    
    fftw_complex * gridTx;
    fftw_complex * gridTy;
    fftw_complex * gridTz;
    
    fftw_complex * gridSxx;
    fftw_complex * gridSxy;
    fftw_complex * gridSxz;
    fftw_complex * gridSyz;
    fftw_complex * gridSyy;
    
    fftw_plan plan;
    
    etemp = calloc(np * 5, sizeof(double));
    
    gridFx = fftw_malloc(nk * sizeof(fftw_complex));
    gridFy = fftw_malloc(nk * sizeof(fftw_complex));
    gridFz = fftw_malloc(nk * sizeof(fftw_complex));
    
    gridTx = fftw_malloc(nk * sizeof(fftw_complex));
    gridTy = fftw_malloc(nk * sizeof(fftw_complex));
    gridTz = fftw_malloc(nk * sizeof(fftw_complex));
    
    gridSxx = fftw_malloc(nk * sizeof(fftw_complex));
    gridSxy = fftw_malloc(nk * sizeof(fftw_complex));
    gridSxz = fftw_malloc(nk * sizeof(fftw_complex));
    gridSyz = fftw_malloc(nk * sizeof(fftw_complex));
    gridSyy = fftw_malloc(nk * sizeof(fftw_complex));
    
    for (i = 0; i < nk; i++) {
        gridFx[i][0] = 0.0;     gridFx[i][1] = 0.0;
        gridFy[i][0] = 0.0;     gridFy[i][1] = 0.0;
        gridFz[i][0] = 0.0;     gridFz[i][1] = 0.0;
        
        gridTx[i][0] = 0.0;     gridTx[i][1] = 0.0;
        gridTy[i][0] = 0.0;     gridTy[i][1] = 0.0;
        gridTz[i][0] = 0.0;     gridTz[i][1] = 0.0;
        
        gridSxx[i][0] = 0.0;    gridSxx[i][1] = 0.0;
        gridSxy[i][0] = 0.0;    gridSxy[i][1] = 0.0;
        gridSxz[i][0] = 0.0;    gridSxz[i][1] = 0.0;
        gridSyz[i][0] = 0.0;    gridSyz[i][1] = 0.0;
        gridSyy[i][0] = 0.0;    gridSyy[i][1] = 0.0;
    }
    
    /*
     Spread force/torque/stresslet of particles to grid
     */
    ewald_recip_spread(ewald_sys,
                       np, p3,
                       gaussCo, gaussExp,
                       f, s,
                       gridFx, gridFy, gridFz,
                       gridTx, gridTy, gridTz,
                       gridSxx, gridSxy, gridSxz, gridSyz, gridSyy);
    /*
     FFT of force/torque/stresslet of grid
     */
    plan = fftw_plan_dft_3d(nkz, nky, nkx, gridFx, gridFx, FFTW_FORWARD, FFTW_ESTIMATE);
    fftw_execute(plan);
    fftw_destroy_plan(plan);
    
    plan = fftw_plan_dft_3d(nkz, nky, nkx, gridFy, gridFy, FFTW_FORWARD, FFTW_ESTIMATE);
    fftw_execute(plan);
    fftw_destroy_plan(plan);
    
    plan = fftw_plan_dft_3d(nkz, nky, nkx, gridFz, gridFz, FFTW_FORWARD, FFTW_ESTIMATE);
    fftw_execute(plan);
    fftw_destroy_plan(plan);
    
    plan = fftw_plan_dft_3d(nkz, nky, nkx, gridTx, gridTx, FFTW_FORWARD, FFTW_ESTIMATE);
    fftw_execute(plan);
    fftw_destroy_plan(plan);
    
    plan = fftw_plan_dft_3d(nkz, nky, nkx, gridTy, gridTy, FFTW_FORWARD, FFTW_ESTIMATE);
    fftw_execute(plan);
    fftw_destroy_plan(plan);
    
    plan = fftw_plan_dft_3d(nkz, nky, nkx, gridTz, gridTz, FFTW_FORWARD, FFTW_ESTIMATE);
    fftw_execute(plan);
    fftw_destroy_plan(plan);
    
    plan = fftw_plan_dft_3d(nkz, nky, nkx, gridSxx, gridSxx, FFTW_FORWARD, FFTW_ESTIMATE);
    fftw_execute(plan);
    fftw_destroy_plan(plan);
    
    plan = fftw_plan_dft_3d(nkz, nky, nkx, gridSxy, gridSxy, FFTW_FORWARD, FFTW_ESTIMATE);
    fftw_execute(plan);
    fftw_destroy_plan(plan);
    
    plan = fftw_plan_dft_3d(nkz, nky, nkx, gridSxz, gridSxz, FFTW_FORWARD, FFTW_ESTIMATE);
    fftw_execute(plan);
    fftw_destroy_plan(plan);
    
    plan = fftw_plan_dft_3d(nkz, nky, nkx, gridSyz, gridSyz, FFTW_FORWARD, FFTW_ESTIMATE);
    fftw_execute(plan);
    fftw_destroy_plan(plan);
    
    plan = fftw_plan_dft_3d(nkz, nky, nkx, gridSyy, gridSyy, FFTW_FORWARD, FFTW_ESTIMATE);
    fftw_execute(plan);
    fftw_destroy_plan(plan);
    
    /*
     Obtain scaled force/torque/stresslet by multiplying Green's function
     */
    ewald_recip_green(ewald_sys,
                      nk,
                      gridFx, gridFy, gridFz,
                      gridTx, gridTy, gridTz,
                      gridSxx, gridSxy, gridSxz, gridSyz, gridSyy);
    
    /*
     IFFT of scaled force/torque/stresslet of grid
     */
    plan = fftw_plan_dft_3d(nkz, nky, nkx, gridFx, gridFx, FFTW_BACKWARD, FFTW_ESTIMATE);
    fftw_execute(plan);
    fftw_destroy_plan(plan);
    
    plan = fftw_plan_dft_3d(nkz, nky, nkx, gridFy, gridFy, FFTW_BACKWARD, FFTW_ESTIMATE);
    fftw_execute(plan);
    fftw_destroy_plan(plan);
    
    plan = fftw_plan_dft_3d(nkz, nky, nkx, gridFz, gridFz, FFTW_BACKWARD, FFTW_ESTIMATE);
    fftw_execute(plan);
    fftw_destroy_plan(plan);
    
    plan = fftw_plan_dft_3d(nkz, nky, nkx, gridTx, gridTx, FFTW_BACKWARD, FFTW_ESTIMATE);
    fftw_execute(plan);
    fftw_destroy_plan(plan);
    
    plan = fftw_plan_dft_3d(nkz, nky, nkx, gridTy, gridTy, FFTW_BACKWARD, FFTW_ESTIMATE);
    fftw_execute(plan);
    fftw_destroy_plan(plan);
    
    plan = fftw_plan_dft_3d(nkz, nky, nkx, gridTz, gridTz, FFTW_BACKWARD, FFTW_ESTIMATE);
    fftw_execute(plan);
    fftw_destroy_plan(plan);
    
    plan = fftw_plan_dft_3d(nkz, nky, nkx, gridSxx, gridSxx, FFTW_BACKWARD, FFTW_ESTIMATE);
    fftw_execute(plan);
    fftw_destroy_plan(plan);
    
    plan = fftw_plan_dft_3d(nkz, nky, nkx, gridSxy, gridSxy, FFTW_BACKWARD, FFTW_ESTIMATE);
    fftw_execute(plan);
    fftw_destroy_plan(plan);
    
    plan = fftw_plan_dft_3d(nkz, nky, nkx, gridSxz, gridSxz, FFTW_BACKWARD, FFTW_ESTIMATE);
    fftw_execute(plan);
    fftw_destroy_plan(plan);
    
    plan = fftw_plan_dft_3d(nkz, nky, nkx, gridSyz, gridSyz, FFTW_BACKWARD, FFTW_ESTIMATE);
    fftw_execute(plan);
    fftw_destroy_plan(plan);
    
    plan = fftw_plan_dft_3d(nkz, nky, nkx, gridSyy, gridSyy, FFTW_BACKWARD, FFTW_ESTIMATE);
    fftw_execute(plan);
    fftw_destroy_plan(plan);
    
    /*
     Interpolate back to particle's center
     */
    ewald_recip_interp(ewald_sys,
                       np, p3,
                       gaussCo, gaussExp,
                       vk,
                       gridFx, gridFy, gridFz,
                       gridTx, gridTy, gridTz,
                       gridSxx, gridSxy, gridSxz, gridSyz, gridSyy,
                       u, etemp);
    
    /*
     Convert Exx, Exy, Exz, Eyz, Eyy to
     Exx - Ezz, 2Exy, 2Exz, 2Eyz, Eyy - Ezz
     */
    for (i = 0; i < np; i++) {
        e[i * 5 + 0] += 2.0 * etemp[i * 5 + 0] + etemp[i * 5 + 4];
        e[i * 5 + 1] += 2.0 * etemp[i * 5 + 1];
        e[i * 5 + 2] += 2.0 * etemp[i * 5 + 2];
        e[i * 5 + 3] += 2.0 * etemp[i * 5 + 3];
        e[i * 5 + 4] += 2.0 * etemp[i * 5 + 4] + etemp[i * 5 + 0];
    }
    
    free(etemp);
    fftw_free(gridFx);
    fftw_free(gridFy);
    fftw_free(gridFz);
    fftw_free(gridTx);
    fftw_free(gridTy);
    fftw_free(gridTz);
    fftw_free(gridSxx);
    fftw_free(gridSxy);
    fftw_free(gridSxz);
    fftw_free(gridSyz);
    fftw_free(gridSyy);
}

/*
 Calculates the mobility problem U = M . F, where M is split into a real part and a reciprocal part, and F = force/torque/stresslet
 Input:
 mode                        mode = 0, FT version. mode = 1, FTS version
 np                          number of particles
 p3                          p3 = p**3, p = number of node in one direction
 nkx, nky, nkz               number of grid points in x,y,z-directions
 nk                          total number of grid points
 gaussCo, gaussExp           coeffcients of Gaussian kernel
 vk                          volume of one grid cell
 self_a, self_c, self_m      self interaction functions, x11a, x11c, x11m
 f                           force/torque of particles
 s                           stresslet of particles
 Output:
 u                           translational/rotational velocity of particles
 e                           rate-of-strain of particles (Exx-Ezz, 2Exy, 2Exz, 2Eyz, Eyy-Ezz)
 */
void ewald_mob_prob_wrap(struct ewald * ewald_sys,
                         const int np, const int p3,
                         const int nkx, const int nky, const int nkz,
                         const int nk,
                         const double gaussCo, const double gaussExp,
                         const double vk,
                         const double self_a, const double self_c, const double self_m,
                         double * f, double * s,
                         double * u, double * e)
{
    int n5, n6;
    int ipair;
    int i, j;
    
    n5 = ewald_sys->tot_pairs * 5;
    n6 = ewald_sys->tot_pairs * 6;
    
    /* Reciprocal part */
    ewald_recip_wrap(ewald_sys,
                     np, p3,
                     nkx, nky, nkz,
                     nk,
                     gaussCo, gaussExp,
                     vk,
                     f, s,
                     u, e);
    
        /* Real part */
        for (ipair = 0; ipair < ewald_sys->tot_pairs; ipair++) {
            i = ewald_sys->row_id[ipair];
            j = ewald_sys->col_id[ipair];
    
            // FT input
            // U = Muf . F + Mus : S
            cblas_dgemv(CblasRowMajor,
                        CblasNoTrans,
                        6, 6,
                        1.0,
                        ewald_sys->muf + ipair * 6, n6,
                        f + j * 6, 1,
                        1.0,
                        u + i * 6, 1);
            cblas_dgemv(CblasRowMajor,
                        CblasTrans,
                        6, 6,
                        1.0,
                        ewald_sys->muf + ipair * 6, n6,
                        f + i * 6, 1,
                        1.0,
                        u + j * 6, 1);
    
            // E = Mef . F + Mes : S
            cblas_dgemv(CblasRowMajor,
                        CblasNoTrans,
                        5, 6,
                        1.0,
                        ewald_sys->mef + ipair * 6, n6,
                        f + j * 6, 1,
                        1.0,
                        e + i * 5, 1);
            cblas_dgemv(CblasRowMajor,
                        CblasTrans,
                        6, 5,
                        1.0,
                        ewald_sys->mus + ipair * 5, n5,
                        f + i * 6, 1,
                        1.0,
                        e + j * 5, 1);
    
            // Mus. S
            cblas_dgemv(CblasRowMajor,
                        CblasNoTrans,
                        6, 5,
                        1.0,
                        ewald_sys->mus + ipair * 5, n5,
                        s + j * 5, 1,
                        1.0,
                        u + i * 6, 1);
            cblas_dgemv(CblasRowMajor,
                        CblasTrans,
                        5, 6,
                        1.0,
                        ewald_sys->mef + ipair * 6, n6,
                        s + i * 5, 1,
                        1.0,
                        u + j * 6, 1);
    
            // Mes : S
            cblas_dgemv(CblasRowMajor,
                        CblasNoTrans,
                        5, 5,
                        1.0,
                        ewald_sys->mes + ipair * 5, n5,
                        s + j * 5, 1,
                        1.0,
                        e + i * 5, 1);
            cblas_dgemv(CblasRowMajor,
                        CblasTrans,
                        5, 5,
                        1.0,
                        ewald_sys->mes + ipair * 5, n5,
                        s + i * 5, 1,
                        1.0,
                        e + j * 5, 1);
    
        } // for ipair
    
        /* Contributino from self part */
        for (i = 0; i < np; i++) {
            u[i * 6 + 0] += self_a * f[i * 6 + 0];
            u[i * 6 + 1] += self_a * f[i * 6 + 1];
            u[i * 6 + 2] += self_a * f[i * 6 + 2];
    
            u[i * 6 + 3] += self_c * f[i * 6 + 3];
            u[i * 6 + 4] += self_c * f[i * 6 + 4];
            u[i * 6 + 5] += self_c * f[i * 6 + 5];
    
            e[i * 5 + 0] += self_m * (2.0 * s[i * 5 + 0] + s[i * 5 + 4]);
            e[i * 5 + 1] += self_m *  2.0 * s[i * 5 + 1];
            e[i * 5 + 2] += self_m *  2.0 * s[i * 5 + 2];
            e[i * 5 + 3] += self_m *  2.0 * s[i * 5 + 3];
            e[i * 5 + 4] += self_m * (2.0 * s[i * 5 + 4] + s[i * 5 + 0]);
        } // for i
    
}

void dipole_real_scalars(double xi ,double r,
                         double * xd, double * yd)
{
    double r2, r3;
    double xi2;
    double xir, xir2;
    double erfcxir, expxir2;
    
    r2 = r * r;
    r3 = r2 * r;
    xi2 = xi * xi;
    xir = xi * r;
    xir2 = xir * xir;
    
    erfcxir = erfc(xir);
    expxir2 = 2.0 * xi / sqrt(M_PI) * exp(- xir2);
    
    (* xd) = - 2.0 / r3 * erfcxir - (2.0 / r2 + 2.0 * xi2) * expxir2;
    (* yd) =  erfcxir / r2 + expxir2 / r2;
}

void dip_a(double xd, double yd,
           double * e,
           double * mb)
{
    int i, j;
    for (i = 0; i < 3; i++) {
        
        for (j = 0; j < 3; j++) {
            mb[(i + 0) * 6 + (j + 3)] = xd * e[i] * e[j] + yd * (kron[i * 3 + j] - e[i] * e[j]);
            mb[(i + 3) * 6 + (j + 0)] = xd * e[i] * e[j] + yd * (kron[i * 3 + j] - e[i] * e[j]);
        } // for j
        
    } // for i
    
}

void dipole_real_mat(struct ewald * ewald_sys,
                     const int np,
                     const double xi)
{
    int n3;
    int ipair;
    int in, jn;
    int npair_ah, pair_acc;
    int pai;
    int i, j;
    int jc;
    
    double xd, yd;
    
    double mb[36];
    
    n3 = ewald_sys->tot_pairs * 3;
    
    ewald_sys->mbd = calloc(3 * n3, sizeof(double));
    
    pair_acc = 0;
    for (in = 0; in < np; in++) {
        npair_ah = ewald_sys->pair_sys[in].ah_count;
        
        for (ipair = 0; ipair < npair_ah; ipair++) {
            jn = ewald_sys->pair_sys[in].ah_list[ipair];
            
            dipole_real_scalars(xi, ewald_sys->pair_sys[in].r[ipair],
                                & xd, & yd);
            
            dip_a(xd, yd,
                  ewald_sys->pair_sys[in].e + ipair * 3,
                  mb);
            
            pai = pair_acc + ipair;
            
            for (i = 0; i < 3; i++) {
                for (j = 0; j < 3; j++) {
                    jc = pai * 3 + j;
                    
                    ewald_sys->mbd[i * n3 + jc] += mb[(i + 0) * 6 + (j + 3)];
                }
            }
            
        } // for ipair
        
        pair_acc += npair_ah;
    } // for in
    
}

/*
 Spreads magnetic dipole onto p**3 grid nodes that are closest to particles
 Input:
 np                                              number of particles
 p3                                              p3 = p**3, p = number of node in one direction
 gaussCo, gaussExp                               coeffcients of Gaussian kernel
 d                                               dipole
 Output:
 gridDx, gridDy, gridDz                          x,y,z-components of dipole of grid node
 */
void dipole_recip_spread(struct ewald * ewald_sys,
                         const int np, const int p3,
                         const double gaussCo, const double gaussExp,
                         double * d,
                         fftw_complex * gridDx, fftw_complex * gridDy, fftw_complex * gridDz)
{
    int i, j;
    int id;
    
    double dx, dy, dz;
    double r2;
    double coeff;
    
    for (i = 0; i < np; i++) {
        dx = d[i * 3 + 0];
        dy = d[i * 3 + 1];
        dz = d[i * 3 + 2];
        
        for (j = 0; j < p3; j++) {
            id = ewald_sys->gridp_sys[i].ids[j];
            r2 = ewald_sys->gridp_sys[i].r2[j];
            
            coeff = gaussCo * exp( - gaussExp * r2 );
            
            gridDx[id][0] += dx * coeff;
            gridDy[id][0] += dy * coeff;
            gridDz[id][0] += dz * coeff;
        } // for j
        
    } // for i
}

/*
 Calculate scaled force/torque of grid by multiplying Green's function
 Input:
 mode                                            mode = 0, calculate field. mode = 1, calculate force
 nk                                              number of grid points
 Input/Output
 gridDx, gridDy, gridDz                          scaled x,y,z-components of dipole of grid node
 */
void dipole_recip_green(struct ewald * ewald_sys,
                        const int nk,
                        fftw_complex * gridDx, fftw_complex * gridDy, fftw_complex * gridDz)
{
    int i;
    
    double kx, ky, kz;
    
    fftw_complex dx, dy, dz;
    fftw_complex dk;
    
    for (i = 0; i < nk; i++) {
        kx = ewald_sys->gridk_sys[i].kx;
        ky = ewald_sys->gridk_sys[i].ky;
        kz = ewald_sys->gridk_sys[i].kz;
        
        dx[0] = gridDx[i][0];       dx[1] = gridDx[i][1];
        dy[0] = gridDy[i][0];       dy[1] = gridDy[i][1];
        dz[0] = gridDz[i][0];       dz[1] = gridDz[i][1];
        
        /* D . k */
        dk[0] = dx[0] * kx + dy[0] * ky + dz[0] * kz;
        dk[1] = dx[1] * kx + dy[1] * ky + dz[1] * kz;
        
        gridDx[i][0] = + ewald_sys->gridk_sys[i].w * kx * dk[0];
        gridDy[i][0] = + ewald_sys->gridk_sys[i].w * ky * dk[0];
        gridDz[i][0] = + ewald_sys->gridk_sys[i].w * kz * dk[0];
        
        gridDx[i][1] = + ewald_sys->gridk_sys[i].w * kx * dk[1];
        gridDy[i][1] = + ewald_sys->gridk_sys[i].w * ky * dk[1];
        gridDz[i][1] = + ewald_sys->gridk_sys[i].w * kz * dk[1];
        
    } // for i
    
}

/*
 Interpolate field of grid points back to particles
 
 Input:
 np                                              number of particles
 p3                                              p3 = p**3, p = number of node in one direction
 gaussCo, gaussExp                               coeffcients of Gaussian kernel
 vk                                              volume of one grid cell
 gridDx, gridDy, gridDz                          scaled x,y,z-components of dipole of grid node
 
 Output:
 b                                               magnetic field
 */
void field_recip_interp(struct ewald * ewald_sys,
                        const int np, const int p3,
                        const double gaussCo, const double gaussExp,
                        const double vk,
                        fftw_complex * gridDx, fftw_complex * gridDy, fftw_complex * gridDz,
                        double * b)
{
    int i, j;
    int id;
    
    double dx, dy, dz;
    double r2;
    double coeff;
    
    for (i = 0; i < np; i++) {
        
        for (j = 0; j < p3; j++) {
            id = ewald_sys->gridp_sys[i].ids[j];
            r2 = ewald_sys->gridp_sys[i].r2[j];
            
            dx = gridDx[id][0];
            dy = gridDy[id][0];
            dz = gridDz[id][0];
            
            coeff = vk * gaussCo * exp( - gaussExp * r2 );
            
            b[i * 3 + 0] += coeff * dx;
            b[i * 3 + 1] += coeff * dy;
            b[i * 3 + 2] += coeff * dz;
            
        } // for j
        
    } // for i
}

/*
 Interpolate force of grid points back to particles
 
 Input:
 np                                              number of particles
 p3                                              p3 = p**3, p = number of node in one direction
 xi2                                             xi**2
 eta                                             Gaussian splitting parameter
 gaussCo, gaussExp                               coeffcients of Gaussian kernel
 vk                                              volume of one grid cell
 gridDx, gridDy, gridDz                          scaled x,y,z-components of dipole of grid node
 
 Output:
 f                                               magnetic force
 */
void force_recip_interp(struct ewald * ewald_sys,
                        const int np, const int p3,
                        const double xi2, const double eta,
                        const double gaussCo, const double gaussExp,
                        const double vk,
                        fftw_complex * gridDx, fftw_complex * gridDy, fftw_complex * gridDz,
                        const double * d,
                        double * f)
{
    int i, j;
    int id;
    
    double dx, dy, dz;
    double dd;
    double x, y, z;
    double r2;
    double coeff;
    
    for (i = 0; i < np; i++) {
        dx = d[i * 3 + 0];
        dy = d[i * 3 + 1];
        dz = d[i * 3 + 2];
        
        for (j = 0; j < p3; j++) {
            id = ewald_sys->gridp_sys[i].ids[j];
            r2 = ewald_sys->gridp_sys[i].r2[j];
            
            x = ewald_sys->gridp_sys[i].dx[j];
            y = ewald_sys->gridp_sys[i].dy[j];
            z = ewald_sys->gridp_sys[i].dz[j];
            
            dd = dx * gridDx[id][0] + dy * gridDy[id][0] + dz * gridDz[id][0];
            
            coeff = 2.0 * gaussCo * vk * (2.0 * xi2 / eta) * exp(- gaussExp * r2);
            f[i * 3 + 0] += - x * coeff * dd;
            f[i * 3 + 1] += - y * coeff * dd;
            f[i * 3 + 2] += - z * coeff * dd;
        } // for j
        
    } // for i
}

/*
 Wraps up all the functions that calculate the reciprocal contribution to magnetic field with input of magnetic dipole
 Input:
 mode                                            mode = 0, calculate field. mode = 1, calculate force
 np                                              number of particles
 p3                                              p3 = p**3, p = number of node in one direction
 nkx, nky, nkz                                   number of grid points in x,y,z-directions
 nk                                              total number of grid points
 gaussCo, gaussExp                               coeffcients of Gaussian kernel
 vk                                              volume of one grid cell
 d                                               dipole of particles
 Output:
 b                                               magnetic field
 */
void field_recip_wrap(struct ewald * ewald_sys,
                      const int np, const int p3,
                      const int nkx, const int nky, const int nkz,
                      const int nk,
                      const double gaussCo, const double gaussExp,
                      const double vk,
                      double * d,
                      double * b)
{
    int i;
    
    fftw_complex * gridDx;
    fftw_complex * gridDy;
    fftw_complex * gridDz;
    
    fftw_plan plan;
    
    gridDx = fftw_malloc(nk * sizeof(fftw_complex));
    gridDy = fftw_malloc(nk * sizeof(fftw_complex));
    gridDz = fftw_malloc(nk * sizeof(fftw_complex));
    
    for (i = 0; i < nk; i++) {
        gridDx[i][0] = 0.0;     gridDx[i][1] = 0.0;
        gridDy[i][0] = 0.0;     gridDy[i][1] = 0.0;
        gridDz[i][0] = 0.0;     gridDz[i][1] = 0.0;
    } // for i
    
    /* Spread dipole of particles to grid */
    dipole_recip_spread(ewald_sys,
                        np, p3,
                        gaussCo, gaussExp,
                        d,
                        gridDx, gridDy, gridDz);
    
    /* FFT of dipole of grid */
    plan = fftw_plan_dft_3d(nkz, nky, nkx, gridDx, gridDx, FFTW_FORWARD, FFTW_ESTIMATE);
    fftw_execute(plan);
    fftw_destroy_plan(plan);
    
    plan = fftw_plan_dft_3d(nkz, nky, nkx, gridDy, gridDy, FFTW_FORWARD, FFTW_ESTIMATE);
    fftw_execute(plan);
    fftw_destroy_plan(plan);
    
    plan = fftw_plan_dft_3d(nkz, nky, nkx, gridDz, gridDz, FFTW_FORWARD, FFTW_ESTIMATE);
    fftw_execute(plan);
    fftw_destroy_plan(plan);
    
    /* Obtain scaled dipole by multiplying Green's function */
    dipole_recip_green(ewald_sys,
                       nk,
                       gridDx, gridDy, gridDz);
    
    /* IFFT of scaled dipole of grid */
    plan = fftw_plan_dft_3d(nkz, nky, nkx, gridDx, gridDx, FFTW_BACKWARD, FFTW_ESTIMATE);
    fftw_execute(plan);
    fftw_destroy_plan(plan);
    
    plan = fftw_plan_dft_3d(nkz, nky, nkx, gridDy, gridDy, FFTW_BACKWARD, FFTW_ESTIMATE);
    fftw_execute(plan);
    fftw_destroy_plan(plan);
    
    plan = fftw_plan_dft_3d(nkz, nky, nkx, gridDz, gridDz, FFTW_BACKWARD, FFTW_ESTIMATE);
    fftw_execute(plan);
    fftw_destroy_plan(plan);
    
    /* Interpolate back to particle's center */
    field_recip_interp(ewald_sys,
                       np, p3,
                       gaussCo, gaussExp,
                       vk,
                       gridDx, gridDy, gridDz,
                       b);
    
    fftw_free(gridDx);
    fftw_free(gridDy);
    fftw_free(gridDz);
}

/*
 Wraps up all functions that calculate the reciprocal part of force
 */
void force_recip_wrap(struct ewald * ewald_sys,
                      const int np, const int p3,
                      const int nkx, const int nky, const int nkz,
                      const int nk,
                      const double xi2, const double eta,
                      const double gaussCo, const double gaussExp,
                      const double vk,
                      double * d,
                      double * f)
{
    int i;
    
    fftw_complex * gridDx;
    fftw_complex * gridDy;
    fftw_complex * gridDz;
    
    fftw_plan plan;
    
    gridDx = fftw_malloc(nk * sizeof(fftw_complex));
    gridDy = fftw_malloc(nk * sizeof(fftw_complex));
    gridDz = fftw_malloc(nk * sizeof(fftw_complex));
    
    for (i = 0; i < nk; i++) {
        gridDx[i][0] = 0.0;     gridDx[i][1] = 0.0;
        gridDy[i][0] = 0.0;     gridDy[i][1] = 0.0;
        gridDz[i][0] = 0.0;     gridDz[i][1] = 0.0;
    }
    
    /* Spread dipole */
    dipole_recip_spread(ewald_sys,
                        np, p3,
                        gaussCo, gaussExp,
                        d,
                        gridDx, gridDy, gridDz);
    
    /* FFT */
    plan = fftw_plan_dft_3d(nkz, nky, nkx, gridDx, gridDx, FFTW_FORWARD, FFTW_ESTIMATE);
    fftw_execute(plan);
    fftw_destroy_plan(plan);
    
    plan = fftw_plan_dft_3d(nkz, nky, nkx, gridDy, gridDy, FFTW_FORWARD, FFTW_ESTIMATE);
    fftw_execute(plan);
    fftw_destroy_plan(plan);
    
    plan = fftw_plan_dft_3d(nkz, nky, nkx, gridDz, gridDz, FFTW_FORWARD, FFTW_ESTIMATE);
    fftw_execute(plan);
    fftw_destroy_plan(plan);
    
    /* Multiply Green's function */
    dipole_recip_green(ewald_sys,
                       nk,
                       gridDx, gridDy, gridDz);
    
    /* IFFT */
    plan = fftw_plan_dft_3d(nkz, nky, nkx, gridDx, gridDx, FFTW_BACKWARD, FFTW_ESTIMATE);
    fftw_execute(plan);
    fftw_destroy_plan(plan);
    
    plan = fftw_plan_dft_3d(nkz, nky, nkx, gridDy, gridDy, FFTW_BACKWARD, FFTW_ESTIMATE);
    fftw_execute(plan);
    fftw_destroy_plan(plan);
    
    plan = fftw_plan_dft_3d(nkz, nky, nkx, gridDz, gridDz, FFTW_BACKWARD, FFTW_ESTIMATE);
    fftw_execute(plan);
    fftw_destroy_plan(plan);
    
    /* Interpolate */
    force_recip_interp(ewald_sys,
                       np, p3,
                       xi2, eta,
                       gaussCo, gaussExp,
                       vk,
                       gridDx, gridDy, gridDz,
                       d,
                       f);
    
    fftw_free(gridDx);
    fftw_free(gridDy);
    fftw_free(gridDz);
}

/*
 Calculate the potential problem b = M . d, where M is split into a real part and a reciprocal part, d is magnetic dipole and b is magnetic field.
 
 Input:
 np                          number of particles
 p3                          p3 = p**3, p = number of node in one direction
 nkx, nky, nkz               number of grid points in x,y,z-directions
 nk                          total number of grid points
 gaussCo, gaussExp           coeffcients of Gaussian kernel
 vk                          volume of one grid cell
 self_d                      self interaction functions, x11d
 d                           magnetic dipole
 
 Output:
 b                           magnetic field
 */
void field_wrap(struct ewald * ewald_sys,
                const int np, const int p3,
                const int nkx, const int nky, const int nkz,
                const int nk,
                const double gaussCo, const double gaussExp,
                const double vk,
                const double self_d,
                double * d,
                double * b)
{
    int n3;
    int ipair;
    int i, j;
    
    n3 = ewald_sys->tot_pairs * 3;
    
    /* Contribution from reciprocal part */
    field_recip_wrap(ewald_sys,
                     np, p3,
                     nkx, nky, nkz,
                     nk,
                     gaussCo, gaussExp,
                     vk,
                     d,
                     b);
    
    /* Contribution from real part */
    for (ipair = 0; ipair < ewald_sys->tot_pairs; ipair++) {
        i = ewald_sys->row_id[ipair];
        j = ewald_sys->col_id[ipair];
        
        cblas_dgemv(CblasRowMajor,
                    CblasNoTrans,
                    3, 3,
                    1.0,
                    ewald_sys->mbd + ipair * 3, n3,
                    d + j * 3, 1,
                    1.0,
                    b + i * 3, 1);
        cblas_dgemv(CblasRowMajor,
                    CblasTrans,
                    3, 3,
                    1.0,
                    ewald_sys->mbd + ipair * 3, n3,
                    d + i * 3, 1,
                    1.0,
                    b + j * 3, 1);
        
    } // for ipair
    
    /* Self contribution */
    for (i = 0; i < np; i++) {
        for (j = 0; j < 3; j++) {
            b[i * 3 + j] += self_d * d[i * 3 + j];
        }
    }
    
}

void force_wrap(struct ewald * ewald_sys,
                const int np,
                const double xi, const double xi2,
                const int p, const int p3,
                const int nkx, const int nky, const int nkz,
                const int nk,
                const double gaussCo, const double gaussExp,
                const double vk, const double eta,
                double * d,
                double * f)
{
    int ipair;
    int i, j;
    int npair;
    
    double ex, ey, ez;
    double r, r2, r3, r4;
    double pisqr;
    double xi4;
    double xir, xir2;
    double erfcxir;
    double expxir2;
    double dae, dbe;
    double dd;
    double xa, ya;
    
    /* Contribution from reciprocal part */
    force_recip_wrap(ewald_sys,
                     np, p3,
                     nkx, nky, nkz,
                     nk,
                     xi2, eta,
                     gaussCo, gaussExp,
                     vk,
                     d,
                     f);
    
    /* Contribution from real part */
    xi4 = xi2 * xi2;
    pisqr = 2.0 * xi / sqrt(M_PI);
    
    for (i = 0; i < np; i++) {
        npair = ewald_sys->pair_sys[i].ah_count;
        
        for (ipair = 0; ipair < npair; ipair++) {
            j = ewald_sys->pair_sys[i].ah_list[ipair];
            
            r = ewald_sys->pair_sys[i].r[ipair];
            ex = - ewald_sys->pair_sys[i].e[ipair * 3    ];
            ey = - ewald_sys->pair_sys[i].e[ipair * 3 + 1];
            ez = - ewald_sys->pair_sys[i].e[ipair * 3 + 2];
            
            r2 = r * r;
            r3 = r2 * r;
            r4 = r2 * r2;
            
            xir = xi * r;
            xir2 = xir * xir;
            erfcxir = erfc(xir);
            expxir2 = pisqr * exp(- xir2);
            
            xa = 3.0 * (erfcxir / r4 + expxir2 * (1 / r3 + xi2 / (3.0 * r)));
            ya = - 15.0 * erfcxir / r4 - expxir2 * (15.0 / r3 + 10.0 * xi2 / r + 4.0 * xi4 * r);
            
            dae = d[i * 3 + 0] * ex + d[i * 3 + 1] * ey + d[i * 3 + 2] * ez;
            dbe = d[j * 3 + 0] * ex + d[j * 3 + 1] * ey + d[j * 3 + 2] * ez;
            dd = d[i * 3 + 0] * d[j * 3 + 0] + d[i * 3 + 1] * d[j * 3 + 1] + d[i * 3 + 2] * d[j * 3 + 2];
            
            f[i * 3 + 0] += xa * (d[i * 3 + 0] * dbe + d[j * 3 + 0] * dae + ex * dd) + ya * (dae * dbe * ex);
            f[i * 3 + 1] += xa * (d[i * 3 + 1] * dbe + d[j * 3 + 1] * dae + ey * dd) + ya * (dae * dbe * ey);
            f[i * 3 + 2] += xa * (d[i * 3 + 2] * dbe + d[j * 3 + 2] * dae + ez * dd) + ya * (dae * dbe * ez);
            
            f[j * 3 + 0] -= xa * (d[i * 3 + 0] * dbe + d[j * 3 + 0] * dae + ex * dd) + ya * (dae * dbe * ex);
            f[j * 3 + 1] -= xa * (d[i * 3 + 1] * dbe + d[j * 3 + 1] * dae + ey * dd) + ya * (dae * dbe * ey);
            f[j * 3 + 2] -= xa * (d[i * 3 + 2] * dbe + d[j * 3 + 2] * dae + ez * dd) + ya * (dae * dbe * ez);
            
        } // for ipair
        
    } // for i
    
}
