//
//  dynamics.c
//  asd_sphere_20200805
//
//  Created by Mingyang Tan on 8/6/20.
//  Copyright © 2020 Mingyang Tan. All rights reserved.
//

#include "dynamics.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <cblas.h>
#include "sd.h"
#include "ewald.h"
#include "lubrication.h"
#include "repulsive.h"
#include "magnetics.h"
#include "iterative.h"
#include "cell.h"

void update_config(struct sd * sys,
                   const double delt,
                   const int ts,
                   const double t,
                   double * d,
                   double * fr, double * fm,
                   double * frx, double * fmx,
                   double * fhn, double * fhf,
                   double * f,
                   double * ulub, double * uoe,
                   double * u, double * s)
{
    int i, j, k;
    /*
     Ewald structure
     1. Initialize structure
     2. Calculate wave value of each grid node
     3. Find p*p*p closest points
     4. Build linked-cell list
     5. Identify pairs within cut-off radius
     6. Build real part of mobility matrix
     */
    struct ewald * ewald_sys = ewald_initialize(sys);
    gridk_value(ewald_sys,
                sys->lx, sys->ly, sys->lz,
                sys->nkx, sys->nky, sys->nkz,
                sys->nk,
                sys->xi2, sys->eta,
                sys->strain);
    id_part_grid(ewald_sys,
                 sys->np,
                 sys->p, sys->p3,
                 sys->xlh, sys->xrh,
                 sys->ylh, sys->yrh,
                 sys->zlh, sys->zrh,
                 sys->lx, sys->ly, sys->lz,
                 sys->nkx, sys->nky, sys->nkz,
                 sys->dkx, sys->dky, sys->dkz,
                 sys->strain, sys->pos);
    
    near_rcell_neighbor_shear(ewald_sys,
                              sys->np,
                              sys->strain);
    rcell_list(ewald_sys,
               sys->np,
               sys->pos);
    mob_pairs_list(ewald_sys,
                   sys->np,
                   sys->xlh, sys->xrh,
                   sys->ylh, sys->yrh,
                   sys->zlh, sys->zrh,
                   sys->lx, sys->ly, sys->lz,
                   sys->rcut2,
                   sys->strain,
                   sys->pos);
    ewald_real_mat(ewald_sys,
                   sys->np,
                   sys->xi);
    dipole_real_mat(ewald_sys,
                    sys->np,
                    sys->xi);
    
    /*
     Lubrication structure
     1. Initialize structure
     2. Build linked-cell list
     3. Identify pairs within cut-off radius
     4. Build lubrication matrix
     5. Build preconditioner
     */
    struct lub * lub_sys = lub_initialize(sys);
    near_lcell_neighbor_shear(lub_sys,
                              sys->np,
                              sys->strain);
    lcell_list(lub_sys,
               sys->np,
               sys->pos);
    lub_pairs_list(lub_sys,
                   sys->np,
                   sys->xlh, sys->xrh,
                   sys->ylh, sys->yrh,
                   sys->zlh, sys->zrh,
                   sys->lx, sys->ly, sys->lz,
                   sys->lcut2,
                   sys->strain,
                   sys->pos);
    lub_mat(lub_sys,
            sys->np,
            sys->lambda);
    lub_precond(lub_sys,
                sys->np,
                sys->lambda);
    
    /* Obtain repulsive force */
    repulsive_hard_potential(lub_sys,
                             sys->np,
                             sys->f0, sys->tau,
                             sys->deltam, sys->delta0,
                             fr);
    
    /* Obtain magnetic force */
    magnetics(sys,
              ewald_sys,
              t,
              d,
              fm);
    
    for (i = 0; i < sys->np; i++) {
        for (j = 0; j < 3; j++) {
            for (k = 0; k < 3; k++) {
                frx[i * 9 + j * 3 + k] = fr[i * 3 + k] * (sys->pos[i * 3 + j]);
                fmx[i * 9 + j * 3 + k] = fm[i * 3 + k] * (sys->pos[i * 3 + j]);
            }
        }
    }
    
    /* Calculate disturbance velocities */
    calc_uos(sys,
             ewald_sys,
             lub_sys,
             fr, fm,
             fhn, fhf,
             f,
             ulub, uoe,
             u, s);
    
    /* Velocities of particles */
    for (i = 0; i < sys->np; i++) {
        u[i * 6 + 0] += (sys->pos[i * 3 + 1]) * (sys->shear_rate);
    }
    
    /* Update positions of particles */
    for (i = 0; i < sys->np; i++) {
        for (j = 0; j < 3; j++) {
            sys->pos[i * 3 + j] += u[i * 6 + j] * delt;
        }
    }
    
    /* Move particles back into simulation box */
    for (i = 0; i < sys->np; i++) {
        return_to_cell(sys->lx, sys->ly, sys->lz,
                       sys->xlh, sys->xrh,
                       sys->ylh, sys->yrh,
                       sys->zlh, sys->zrh,
                       sys->strain,
                       & sys->pos[i * 3]);
    }
    
    // calculate viscosity
    double stot = 0.0;
    for (i = 0; i < sys->np; i++) {
        stot += s[i * 5 + 1];
    }
    sys->viscosity = 1.0 + 6.0 * M_PI * stot / sys->volume / sys->shear_rate;
    
    
    /* Free memories */
    free(ewald_sys->row_id);
    free(ewald_sys->col_id);
    free(ewald_sys->index);
    free(ewald_sys->head);
    free(ewald_sys->list);
    free(ewald_sys->muf);
    free(ewald_sys->mus);
    free(ewald_sys->mef);
    free(ewald_sys->mes);
    free(ewald_sys->mbd);
    for (i = 0; i < ewald_sys->nd; i++) {
        free(ewald_sys->cell_sys[i].inter_box);
    }
    free(ewald_sys->cell_sys);
    free(ewald_sys->gridk_sys);
    for (i = 0; i < sys->np; i++) {
        free(ewald_sys->gridp_sys[i].ids);
        free(ewald_sys->gridp_sys[i].r2);
        free(ewald_sys->gridp_sys[i].dx);
        free(ewald_sys->gridp_sys[i].dy);
        free(ewald_sys->gridp_sys[i].dz);
    }
    free(ewald_sys->gridp_sys);
    
    for (i = 0; i < sys->np; i++) {
        free(ewald_sys->pair_sys[i].bh_list);
        free(ewald_sys->pair_sys[i].ah_list);
        free(ewald_sys->pair_sys[i].r);
        free(ewald_sys->pair_sys[i].e);
    }
    free(ewald_sys->pair_sys);
    free(ewald_sys);
    
    free(lub_sys->row_id);
    free(lub_sys->col_id);
    free(lub_sys->index);
    free(lub_sys->head);
    free(lub_sys->list);
    free(lub_sys->lfu);
    free(lub_sys->lfe);
    free(lub_sys->lsu);
    free(lub_sys->lse);
    free(lub_sys->iccl);
    free(lub_sys->lfu_l);
    free(lub_sys->lfu_u);
    for (i = 0; i < lub_sys->nd; i++) {
        free(lub_sys->cell_sys[i].inter_box);
    }
    free(lub_sys->cell_sys);
    
    for (i = 0; i < sys->np; i++) {
        free(lub_sys->pair_sys[i].bh_list);
        free(lub_sys->pair_sys[i].ah_list);
        free(lub_sys->pair_sys[i].r);
        free(lub_sys->pair_sys[i].e);
    }
    free(lub_sys->pair_sys);
    free(lub_sys);
}

void calc_uos(struct sd * sys,
              struct ewald * ewald_sys,
              struct lub * lub_sys,
              double * fr, double * fm,
              double * fhn, double * fhf,
              double * f,
              double * ulub, double * uoe,
              double * u, double * s)
{
    int n5, n6;
    n5 = sys->np * 5;
    n6 = sys->np * 6;
    
    int n11;
    double * fts;
    n11 = sys->np * 11;
    fts = calloc(n11, sizeof(double));
    
    int ni5, nb5, nb6;
    ni5 = lub_sys->tot_pairs * 5;
    nb5 = (lub_sys->tot_pairs + sys->np) * 5;
    nb6 = (lub_sys->tot_pairs + sys->np) * 6;
    
    int i, j;
    int ipair;
    
    /* F = F + Lfe : E */
    for (i = 0; i < sys->np; i++) {
        cblas_dgemv(CblasRowMajor,
                    CblasTrans,
                    5, 6,
                    1.0,
                    lub_sys->lsu + i * 6, nb6,
                    sys->einf, 1,
                    1.0,
                    fhn + i * 6, 1);
    }
    for (ipair = 0; ipair < lub_sys->tot_pairs; ipair++) {
        i = lub_sys->row_id[ipair];
        j = lub_sys->col_id[ipair];
        
        cblas_dgemv(CblasRowMajor,
                    CblasTrans,
                    5, 6,
                    1.0,
                    lub_sys->lsu + (sys->np + ipair) * 6, nb6,
                    sys->einf, 1,
                    1.0,
                    fhn + j * 6, 1);
        cblas_dgemv(CblasRowMajor,
                    CblasNoTrans,
                    6, 5,
                    1.0,
                    lub_sys->lfe + ipair * 5, ni5,
                    sys->einf, 1,
                    1.0,
                    fhn + i * 6, 1);
    } // for ipair
    
    for (i = 0; i < sys->np; i++) {
        for (j = 0; j < 3; j++) {
            f[i * 6 + j] = fr[i * 3 + j] + fm[i * 3 + j];
        }
    }
    
    for (i = 0; i < sys->np; i++) {
        for (j = 0; j < 6; j++) {
            f[i * 6 + j] += fhn[i * 6 + j];
        }
    }
    if (lub_sys->precond == 1) {
        pcg_iccf(lub_sys,
                 sys->np,
                 sys->itmax, sys->tol,
                 f,
                 ulub);
    } else {
        cg(lub_sys,
           sys->np,
           sys->itmax, sys->tol,
           f,
           ulub);
    }
    
    ewald_mob_prob_wrap(ewald_sys,
                        sys->np, sys->p3,
                        sys->nkx, sys->nky, sys->nkz,
                        sys->nk,
                        sys->gaussCo, sys->gaussExp,
                        sys->vk,
                        sys->self_a, sys->self_c, sys->self_m,
                        ulub, s,
                        uoe, uoe + n6);
    for (i = 0; i < n6; i++) {
        uoe[i] = sys->lambda * uoe[i] - ulub[i];
    }
    for (i = 0; i < sys->np; i++) {
        for (j = 0; j < 5; j++) {
            uoe[n6 + i * 5 + j] = sys->lambda * uoe[n6 + i * 5 + j] + sys->einf[j];
        }
    }
    
    pgmres_mob(sys,
               ewald_sys,
               lub_sys,
               uoe,
               fts);
    for (i = 0; i < n6; i++) {
        f[i] += fts[i];
    }
    
    for (i = 0; i < n6; i++) {
        fhf[i] = fts[i];
    }
    
    // disturbance velocity
    if (lub_sys->precond == 1) {
        pcg_iccf(lub_sys,
                 sys->np,
                 sys->itmax, sys->tol,
                 f,
                 u);
    } else {
        cg(lub_sys,
           sys->np,
           sys->itmax, sys->tol,
           f,
           u);
    }
    
    // calculate stresslet
    for (i = 0; i < sys->np; i++) {
        cblas_dgemv(CblasRowMajor,
                    CblasNoTrans,
                    5, 6,
                    -1.0,
                    lub_sys->lsu + i * 6, nb6,
                    u + i * 6, 1,
                    1.0,
                    s + i * 5, 1);
    } // for i
    for (ipair = 0; ipair < lub_sys->tot_pairs; ipair++) {
        i = lub_sys->row_id[ipair];
        j = lub_sys->col_id[ipair];
        
        cblas_dgemv(CblasRowMajor,
                    CblasNoTrans,
                    5, 6,
                    -1.0,
                    lub_sys->lsu + (sys->np + ipair) * 6, nb6,
                    u + j * 6, 1,
                    1.0,
                    s + i * 5, 1);
        cblas_dgemv(CblasRowMajor,
                    CblasTrans,
                    6, 5,
                    -1.0,
                    lub_sys->lfe + ipair * 5, ni5,
                    u + i * 6, 1,
                    1.0,
                    s + j * 5, 1);
    } // for ipair
    
    for (i = 0; i < sys->np; i++) {
        cblas_dgemv(CblasRowMajor,
                    CblasNoTrans,
                    5, 5,
                    1.0,
                    lub_sys->lse + i * 5, nb5,
                    sys->einf, 1,
                    1.0,
                    s + i * 5, 1);
    }
    for (ipair = 0; ipair < lub_sys->tot_pairs; ipair++) {
        i = lub_sys->row_id[ipair];
        j = lub_sys->col_id[ipair];
        
        cblas_dgemv(CblasRowMajor,
                    CblasNoTrans,
                    5, 5,
                    1.0,
                    lub_sys->lse + (sys->np + ipair) * 5, nb5,
                    sys->einf, 1,
                    1.0,
                    s + i * 5, 1);
        cblas_dgemv(CblasRowMajor,
                    CblasTrans,
                    5, 5,
                    1.0,
                    lub_sys->lse + (sys->np + ipair) * 5, nb5,
                    sys->einf, 1,
                    1.0,
                    s + j * 5, 1);
    }
    for (i = 0; i < n5; i++) {
        s[i] += fts[n6 + i];
    }
    
    free(fts);
}
