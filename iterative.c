//
//  iterative.c
//  asd_sphere_20200508
//
//  Created by Mingyang Tan on 5/11/20.
//  Copyright © 2020 Mingyang Tan. All rights reserved.
//

#include "iterative.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <cblas.h>
#include "lubrication.h"
#include "ewald.h"
#include "sd.h"

static void forw_sub_iccl(struct lub * lub_sys,
                          const int np,
                          const double * iccl, const double * b,
                          double * x)
{
    int nb;
    int nb6;
    int ipair;
    int in, jn;
    int npair_bh;
    int indent_i, ind_i, block_ind_i;
    int i, j;
    
    double s;
    
    nb = np + lub_sys->tot_pairs;
    nb6 = nb * 6;
    
    for (in = 0; in < np; in++) {
        npair_bh = lub_sys->pair_sys[in].bh_count;
        
        for (i = 0; i < 6; i++) {
            s = 0.0;
            
            for (ipair = 0; ipair < npair_bh; ipair++) {
                jn = lub_sys->pair_sys[in].bh_list[ipair];
                indent_i = jn * (jn + 1) / 2;
                ind_i = jn * np + in - indent_i;
                block_ind_i = lub_sys->index[ind_i];
                
                for (j = 0; j < 6; j++) {
                    s += iccl[block_ind_i * 6 + i * nb6 + j] * x[jn * 6 + j];
                } // for j
                
            } // for ipair
            
            for (j = 0; j < i; j++) {
                s += iccl[in * 6 + i * nb6 + j] * x[in * 6 + j];
            } // for j
            
            x[in * 6 + i] = (b[in * 6 + i] - s) / iccl[in * 6 + i * nb6 + i];
            
        } // for i
        
    } // for in
    
}

static void back_sub_iccl(struct lub * lub_sys,
                          const int np,
                          const double * iccl, const double * b,
                          double * x)
{
    int nb;
    int nb6;
    int ipair;
    int in, jn;
    int npair_ah;
    int indent_i, ind_i, block_ind_i;
    int i, j;
    
    double s;
    
    nb = np + lub_sys->tot_pairs;
    nb6 = nb * 6;
    
    for (in = np - 1; in >= 0; in--) {
        npair_ah = lub_sys->pair_sys[in].ah_count;
        indent_i = in * (in + 1) / 2;
        
        for (i = 5; i >= 0; i--) {
            s = 0.0;
            
            for (ipair = 0; ipair < npair_ah; ipair++) {
                jn = lub_sys->pair_sys[in].ah_list[ipair];
                ind_i = in * np + jn - indent_i;
                block_ind_i = lub_sys->index[ind_i];
                
                for (j = 5; j >= 0; j--) {
                    s += iccl[block_ind_i * 6 + j * nb6 + i] * x[jn * 6 + j];
                }
                
            } // for ipair
            
            for (j = 5; j > i; j--) {
                s += iccl[in * 6 + j * nb6 + i] * x[in * 6 + j];
            } // for j
            
            x[in * 6 + i] = (b[in * 6 + i] - s) / iccl[in * 6 + i * nb6 + i];
            
        } // for i
        
    } // for in
    
}

void pcg_iccf(struct lub * lub_sys,
              const int np,
              const int itmax, const double tol,
              const double * b,
              double * x)
{
    int nb;
    int n6;
    int nb6;
    
    int ipair;
    int in, jn;
    
    double * r;
    double * z;
    double * w;
    double rr;
    double rz1, rz2;
    
    double * p;
    double * ap;
    double pap;
    
    double al;
    double be;
    
    double res;
    
    int i, j, k;
    
    res = 0.0;
    
    nb = np + lub_sys->tot_pairs;
    n6 = np * 6;
    nb6 = nb * 6;
    
    r = calloc(n6, sizeof(double));
    z = calloc(n6, sizeof(double));
    p = calloc(n6, sizeof(double));
    ap = calloc(n6, sizeof(double));
    w = calloc(n6, sizeof(double));
    
    /* Initial guess */
    for (i = 0; i < n6; i++) {
        r[i] += b[i];
    }
    forw_sub_iccl(lub_sys,
                  np,
                  lub_sys->iccl, r,
                  w);
    
    back_sub_iccl(lub_sys,
                  np,
                  lub_sys->iccl, w,
                  z);
    
    for (i = 0; i < n6; i++) {
        p[i] = z[i];
    }
    
    for (j = 0; j < itmax; j++) {
        rz1 = cblas_ddot(n6,
                         r, 1,
                         z, 1);
        
        // A . p
        for (i = 0; i < n6; i++) {
            ap[i] = 0.0;
        }
        for (in = 0; in < np; in++) {
            cblas_dgemv(CblasRowMajor,
                        CblasNoTrans,
                        6, 6,
                        1.0,
                        lub_sys->lfu + in * 6, nb6,
                        p + in * 6, 1,
                        1.0,
                        ap + in * 6, 1);
        }
        for (ipair = 0; ipair < lub_sys->tot_pairs; ipair++) {
            in = lub_sys->row_id[ipair];
            jn = lub_sys->col_id[ipair];
            
            cblas_dgemv(CblasRowMajor,
                        CblasNoTrans,
                        6, 6,
                        1.0,
                        lub_sys->lfu + (np + ipair) * 6, nb6,
                        p + jn * 6, 1,
                        1.0,
                        ap + in * 6, 1);
            cblas_dgemv(CblasRowMajor,
                        CblasTrans,
                        6, 6,
                        1.0,
                        lub_sys->lfu + (np + ipair) * 6, nb6,
                        p + in * 6, 1,
                        1.0,
                        ap + jn * 6, 1);
        }  // for ipair
        pap = cblas_ddot(n6,
                         p, 1,
                         ap, 1);
        
        al = rz1 / pap;
        for (i = 0; i < n6; i++) {
            x[i] += al * p[i];
            r[i] -= al * ap[i];
        }
        rr = cblas_ddot(n6,
                        r, 1,
                        r, 1);
        
        res = sqrt(rr);
        if (res <= tol) {
            j++;
            break;
        }
        for (i = 0; i < n6; i++) {
            w[i] = 0.0;
            z[i] = 0.0;
        }
        forw_sub_iccl(lub_sys,
                      np,
                      lub_sys->iccl, r,
                      w);
        back_sub_iccl(lub_sys,
                      np,
                      lub_sys->iccl, w,
                      z);
        rz2 = cblas_ddot(n6,
                         r, 1,
                         z, 1);
        be = rz2 / rz1;
        for (i = 0; i < n6; i++) {
            p[i] = z[i] + be * p[i];
        }
    }
    
    printf("\nNumber of PCG iterations %d.\n", j);
    if (j == itmax) {
        printf("Number of PCG iterations reaches max.\n");
    }
    
    free(z);
    free(r);
    free(p);
    free(ap);
    free(w);
}

void cg(struct lub * lub_sys,
        const int np,
        const int itmax, const double tol,
        const double * b,
        double * x)
{
    int nb;
    int n6;
    int nb6;
    
    int ipair;
    int in, jn;
    
    double *r1;
    double *r2;
    double rr1, rr2;
    
    double *p;
    double *ap;
    double pap;
    
    double al;
    double be;
    
    double res;
    
    int i, j;
    
    res = 0.0;
    
    nb = np + lub_sys->tot_pairs;
    n6 = np * 6;
    nb6 = nb * 6;
    
    r1 = calloc(n6, sizeof(double));
    r2 = calloc(n6, sizeof(double));
    p = calloc(n6, sizeof(double));
    ap = calloc(n6, sizeof(double));
    
    /* Initial guess */
    for (i = 0; i < n6; i++) {
        r1[i] += b[i];
        p[i] = r1[i];
    }
    
    /* Start iteration */
    for (j = 0; j < itmax; j++) {
        // Calculate r1 . r1
        rr1 = cblas_ddot(n6,
                         r1, 1,
                         r1, 1);
        
        // Calculate A . p
        for (i = 0; i < n6; i++) {
            ap[i] = 0.0;
        }
        for (in = 0; in < np; in++) {
            cblas_dgemv(CblasRowMajor,
                        CblasNoTrans,
                        6, 6,
                        1.0,
                        lub_sys->lfu + in * 6, nb6,
                        p + in * 6, 1,
                        1.0,
                        ap + in * 6, 1);
        } // for in
        
        for (ipair = 0; ipair < lub_sys->tot_pairs; ipair++) {
            in = lub_sys->row_id[ipair];
            jn = lub_sys->col_id[ipair];
            
            cblas_dgemv(CblasRowMajor,
                        CblasNoTrans,
                        6, 6,
                        1.0,
                        lub_sys->lfu + (np + ipair) * 6, nb6,
                        p + jn * 6, 1,
                        1.0,
                        ap + in * 6, 1);
            cblas_dgemv(CblasRowMajor,
                        CblasTrans,
                        6, 6,
                        1.0,
                        lub_sys->lfu + (np + ipair) * 6, nb6,
                        p + in * 6, 1,
                        1.0,
                        ap + jn * 6, 1);
        }  // for ipair
        
        // Calculate p . Ap
        pap = cblas_ddot(n6,
                         ap, 1,
                         p, 1);
        al = rr1 / pap;
        
        // Calculate x^{k+1} = x^{k} + al * p^{k}
        // r^{k+1} = r^{k} - al * L^{-1} . A . p^{k}
        for (i = 0; i < n6; i++) {
            x[i] += al * p[i];
            r2[i] = r1[i] - al * ap[i];
        }
        
        rr2 = cblas_ddot(n6,
                         r2, 1,
                         r2, 1);
        res = sqrt(rr2);
        if (res <= tol) {
            j++;
            break;
        }
        
        be = rr2 / rr1;
        
        for (i = 0; i < n6; i++) {
            p[i] = r2[i] + be * p[i];
            r1[i] = r2[i];
        }
        
    } // for j
    
    printf("\nNumber of CG iterations %d.\n", j);
    if (j == itmax) {
        printf("Number of CG iterations reaches max.\n");
    }
    
    free(r1);
    free(r2);
    free(p);
    free(ap);
}

static void forw_sub_l(const int m,
                       const double * L, const int ldl,
                       const double * b,
                       double * x)
{
    int i, j;
    
    double s;
    
    for (i = 0; i < m; i++) {
        s = 0.0;
        
        for (j = 0; j < i; j++) {
            s += L[i * ldl + j] * x[j];
        } // for j
        
        x[i] = (b[i] - s) / L[i * ldl + i];
    } // for i
    
}

static void back_sub_u(const int m,
                       const double * U, const int ldu,
                       const double * b,
                       double * x)
{
    int i, j;
    
    double s;
    
    for (i = m - 1; i >= 0; i--) {
        s = 0.0;
        
        for (j = m - 1; j > i; j--) {
            s += U[i * ldu + j] * x[j];
        } // for j
        
        x[i] = (b[i] - s) / U[i * ldu + i];
    } // for i
    
}

static void minimize_y(int m, int n,
                       const double * r, const double * g,
                       double * y)
{
    int i, j, k;
    
    for (j = m - 1, k = 0; k < m; j--, k++) {
        y[j] = g[j];
        
        for (i = j + 1; i < m; i++) {
            y[j] -= r[j * n + i] * y[i];
        }
        y[j] = y[j] / r[j * n + j];
        
    } // for j, k
}

void pgmres_mob(struct sd * sys,
                struct ewald * ewald_sys,
                struct lub * lub_sys,
                const double * b,
                double * x)
{
    int n5, n6, n11;
    int i, j, k;
    int m;
    
    double xm;
    double hh, tt, hv;
    double t1, t2;
    double g0;
    double res;
    
    double * uff;
    double * w;
    double * v;
    double * z;
    double * h;
    double * g;
    double * cf;
    double * sf;
    
    n5 = sys->np * 5;
    n6 = sys->np * 6;
    n11 = n5 + n6;
    m = sys->itmax;
    
    xm = 0.9;
    
    uff = calloc(n6, sizeof(double));
    w = calloc(n6, sizeof(double));
    v = calloc((m + 1) * n11, sizeof(double));
    z = calloc((m + 1) * n11, sizeof(double));
    h = calloc(m * m, sizeof(double));
    g = calloc(m + 1, sizeof(double));
    cf = calloc(m, sizeof(double));
    sf = calloc(m, sizeof(double));
    
    /* Initial guess */
    res = 0.0;
    for (i = 0; i < n11; i++) {
        v[i] += b[i];
    }
    g[0] = cblas_ddot(n11,
                      v + 0, 1,
                      v + 0, 1);
    g[0] = sqrt(g[0]);
    cblas_dscal(n11,
                1.0 / g[0],
                v + 0, 1);
    
    /* Start iteration */
    for (j = 0; j < m; j++) {
        // multiply preconditioner
        for (i = 0; i < n6; i++) {
            w[i] = 0.0;
        }
        for (i = 0; i < sys->np; i++) {
            forw_sub_l(6,
                       & lub_sys->lfu_l[i * 6], n6,
                       v + j * n11 + i * 6,
                       w + i * 6);
            back_sub_u(6,
                       & lub_sys->lfu_u[i * 6], n6,
                       w + i * 6,
                       z + j * n11 + i * 6);
            
            z[j * n11 + n6 + i * 5 + 0] = (2.0 * v[j * n11 + n6 + i * 5 + 0] - v[j * n11 + n6 + i * 5 + 4]) / (3.0 * xm);
            z[j * n11 + n6 + i * 5 + 1] =        v[j * n11 + n6 + i * 5 + 1] / (2.0 * xm);
            z[j * n11 + n6 + i * 5 + 2] =        v[j * n11 + n6 + i * 5 + 2] / (2.0 * xm);
            z[j * n11 + n6 + i * 5 + 3] =        v[j * n11 + n6 + i * 5 + 3] / (2.0 * xm);
            z[j * n11 + n6 + i * 5 + 4] = (2.0 * v[j * n11 + n6 + i * 5 + 4] - v[j * n11 + n6 + i * 5 + 0]) / (3.0 * xm);
        } // for i
        
        // uff = Lfu^{-1} . F
        for (i = 0; i < n6; i++) {
            uff[i] = 0.0;
        }
        if (lub_sys->precond == 1) {
            pcg_iccf(lub_sys,
                     sys->np,
                     sys->itmax, sys->tol,
                     z + j * n11,
                     uff);
        } else {
            cg(lub_sys,
               sys->np,
               sys->itmax, sys->tol,
               z + j * n11,
               uff);
        }
       
        // -\lambda * M . [uff, 0]^{T}
        ewald_mob_prob_wrap(ewald_sys,
                            sys->np, sys->p3,
                            sys->nkx, sys->nky, sys->nkz,
                            sys->nk,
                            sys->gaussCo, sys->gaussExp,
                            sys->vk,
                            sys->self_a, sys->self_c, sys->self_m,
                            uff, x + n6,
                            v + (j + 1) * n11, v + (j + 1) * n11 + n6);
        
        cblas_dscal(n11,
                    - sys->lambda,
                    v + (j + 1) * n11, 1);
        
        // M . [Fff, Sff]^{T}
        ewald_mob_prob_wrap(ewald_sys,
                            sys->np, sys->p3,
                            sys->nkx, sys->nky, sys->nkz,
                            sys->nk,
                            sys->gaussCo, sys->gaussExp,
                            sys->vk,
                            sys->self_a, sys->self_c, sys->self_m,
                            z + j * n11, z + j * n11 + n6,
                            v + (j + 1) * n11, v + (j + 1) * n11 + n6);
        
        
        for (i = 0; i < n6; i++) {
            v[(j + 1) * n11 + i] += uff[i];
        }
        
        // Calculate h_{ij}
        for (i = 0; i <= j; i++) {
            h[i * m + j] = cblas_ddot(n11,
                                      v + (j + 1) * n11, 1,
                                      v + i * n11, 1);
        }
        
        // Calculate v_{j+1} = Av_{j} - h_{ij} . v_{i}
        // v_{j+1} = v_{j+1} / || v_{j+1} ||
        for (k = 0; k < n11; k++) {
            hv = 0.0;
            for (i = 0; i <= j; i++) {
                hv += h[i * m + j] * v[i * n11 + k];
            }
            v[(j + 1) * n11 + k] -= hv;
        } // for k
        
        hh = cblas_ddot(n11,
                        v + (j + 1) * n11, 1,
                        v + (j + 1) * n11, 1);
        hh = sqrt(hh);
        cblas_dscal(n11,
                    1.0 / hh,
                    v + (j + 1) * n11, 1);
        
        // Rotation
        for (i = 0; i < j; i++) {
            t1 = h[ i      * m + j];
            t2 = h[(i + 1) * m + j];
            h[ i      * m + j] = cf[i] * t1 - sf[i] * t2;
            h[(i + 1) * m + j] = sf[i] * t1 + cf[i] * t2;
        } // for i
        tt = h[j * m + j];
        hv = sqrt(tt * tt + hh * hh);
        cf[j] = + tt / hv;
        sf[j] = - hh / hv;
        h[j * m + j] = hv;
        
        g0 = g[j];
        g[j    ] = cf[j] * g0;
        g[j + 1] = sf[j] * g0;
        
        res = fabs(g[j + 1]);
        if (res <= sys->tol) {
            j++;
            break;
        }
        
    } // for j
    
    if (j == sys->itmax) {
        printf("Number of mob GMRES iterations reaches max.\n");
    }
    printf("Number of mob GMRES iterations %d.\n", j);
    
    minimize_y(j, m,
               h, g,
               cf);
    for (i = 0; i < n11; i++) {
        for (k = 0; k < j; k++) {
            x[i] += z[k * n11 + i] * cf[k];
        }
    }
    
    free(uff);
    free(w);
    free(v);
    free(z);
    free(h);
    free(g);
    free(cf);
    free(sf);
}

void pgmres_mag(struct sd * sys,
                struct ewald * ewald_sys,
                const double * b,
                double * x)
{
    int n3;
    int i, j, k;
    int m;
    
    double hh, tt, hv;
    double t1, t2;
    double g0;
    double res;
    
    double * v;
    double * h;
    double * g;
    double * cf;
    double * sf;
    
    n3 = sys->np * 3;
    m = sys->itmax;
    
    v = calloc((m + 1) * n3, sizeof(double));
    h = calloc(m * m, sizeof(double));
    g = calloc(m + 1, sizeof(double));
    cf = calloc(m, sizeof(double));
    sf = calloc(m, sizeof(double));
    
    /* Initial guess */
    res = 0.0;
    for (i = 0; i < sys->np; i++) {
        for (j = 0; j < 3; j++) {
            v[i * 3 + j] += b[j];
        }
    }
    g[0] = cblas_ddot(n3,
                      v + 0, 1,
                      v + 0, 1);
    g[0] = sqrt(g[0]);
    cblas_dscal(n3,
                1.0 / g[0],
                v + 0, 1);
    
    /* Start iteration */
    for (j = 0; j < m; j++) {
        field_wrap(ewald_sys,
                   sys->np, sys->p3,
                   sys->nkx, sys->nky, sys->nkz,
                   sys->nk,
                   sys->gaussCo, sys->gaussExp,
                   sys->vk,
                   sys->self_d,
                   v + j * n3,
                   v + (j + 1) * n3);
        
        // Calculate h_{ij}
        for (i = 0; i <= j; i++) {
            h[i * m + j] = cblas_ddot(n3,
                                      v + (j + 1) * n3, 1,
                                      v + i * n3, 1);
        }
        
        // Calculate v_{j+1} = Av_{j} - h_{ij} . v_{i}
        // v_{j+1} = v_{j+1} / || v_{j+1} ||
        for (k = 0; k < n3; k++) {
            hv = 0.0;
            for (i = 0; i <= j; i++) {
                hv += h[i * m + j] * v[i * n3 + k];
            }
            v[(j + 1) * n3 + k] -= hv;
        } // for k
        
        hh = cblas_ddot(n3,
                        v + (j + 1) * n3, 1,
                        v + (j + 1) * n3, 1);
        hh = sqrt(hh);
        cblas_dscal(n3,
                    1.0 / hh,
                    v + (j + 1) * n3, 1);
        
        // Rotation
        for (i = 0; i < j; i++) {
            t1 = h[ i      * m + j];
            t2 = h[(i + 1) * m + j];
            h[ i      * m + j] = cf[i] * t1 - sf[i] * t2;
            h[(i + 1) * m + j] = sf[i] * t1 + cf[i] * t2;
        } // for i
        tt = h[j * m + j];
        hv = sqrt(tt * tt + hh * hh);
        cf[j] = + tt / hv;
        sf[j] = - hh / hv;
        h[j * m + j] = hv;
        
        g0 = g[j];
        g[j    ] = cf[j] * g0;
        g[j + 1] = sf[j] * g0;
        
        res = fabs(g[j + 1]);
        if (res <= sys->tol) {
            j++;
            break;
        }
        
    } // for j
    
    if (j == sys->itmax) {
        printf("Number of mag GMRES iterations reaches max.\n");
    }
    printf("\nNumber of mag GMRES iterations %d.\n", j);
    
    minimize_y(j, m,
               h, g,
               cf);
    
    for (i = 0; i < n3; i++) {
        for (k = 0; k < j; k++) {
            x[i] += v[k * n3 + i] * cf[k];
        }
    }
    
    free(v);
    free(h);
    free(g);
    free(cf);
    free(sf);
}
