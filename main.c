//
//  main.c
//  asd_sphere_20200805
//
//  Created by Mingyang Tan on 8/5/20.
//  Copyright © 2020 Mingyang Tan. All rights reserved.
//

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <sys/stat.h>
#include "sd.h"
#include "global.h"
#include "dynamics.h"

int main(int argc, const char * argv[]) {
    glob_variables();
    
    struct sd * sys = sd_initialize();
    int i, j;
    
    // number of particles
    int m;
    int np;
    int n3, n5, n6, n9, n11;
    m = strtod(argv[1], NULL);
    np = m * m * m;
    n3 = np * 3;
    n5 = np * 5;
    n6 = np * 6;
    n9 = np * 9;
    n11 = np * 11;
    sd_np(sys,
          np);
    
    // size of simulation box
    double phi;
    double volume;
    double lx, ly, lz;
    phi = strtod(argv[2], NULL);
    volume = 4.0 / 3.0 * M_PI * (double) np / phi;
    lx = pow(volume, 1.0 / 3.0);
    ly = lx;
    lz = lx;
    sd_cell(sys,
            phi,
            lx, ly, lz);
    
    // position of particles
    double * pos;
    pos = malloc(n3 * sizeof(double));
    FILE * posin = fopen("posini.dat", "r");
    int count = 0;
    while (!feof(posin) && (count < n3)) {
        fscanf(posin, "%lf", &(pos[count++]));
    }
    fclose(posin);
    sd_pos(sys,
           pos);
    
    // parameters of Ewald summation
    double error;
    double xi;
    double rcut;
    int kcut;
    int nkx, nky, nkz;
    error = 1.0e-3;
    sd_set_xi(error,
              lx, ly, lz,
              & xi,
              & rcut, & kcut);
    sd_set_fft_point(kcut,
                     lx, ly, lz,
                     & nkx, & nky, & nkz);
    sd_ewald(sys,
             xi);
    sd_rcut(sys,
            rcut);
    sd_fft(sys,
           nkx, nky, nkz);
    printf("\nxi = %f\nrcut = %f\nnkx = %d, nky = %d, nkz = %d\n", sys->xi, sys->rcut, sys->nkx, sys->nky, sys->nkz);
    
    // parameters of spectral accuracy
    double eta;
    int p;
    sd_set_eta(sys,
               error,
               & eta, & p);
    sd_gauss(sys,
             eta, p);
    printf("\neta = %f\nP = %d\n", sys->eta, sys->p);
    
    // cutoff radius of lubrication interaction
    double lcut;
    lcut = 4.0;
    sd_lcut(sys,
            lcut);
    
    // parameters of iterative solver
    double lambda;
    int itmax;
    double tol;
    lambda = 4.0;
    itmax = 1000;
    //    tol = 1.0e-3;
    tol = error;
    sd_iterative(sys,
                 lambda,
                 itmax, tol);
    
    // external flow field
    double shear_rate;
    double ux, uy, uz;
    double ox, oy, oz;
    shear_rate = 1.0;
    ux = 0.0;   uy = 0.0;   uz = 0.0;
    ox = 0.0;   oy = 0.0;   oz = 0.0;
    sd_uinf(sys,
            shear_rate,
            ux, uy, uz,
            ox, oy, oz);
    
    // magnetic properties
    double Mam, Man;
    double Ma;
    double omem, omen;
    double omega;
    int direct;
    Mam = strtod(argv[3], NULL);
    Man = strtod(argv[4], NULL);
    Ma = Mam * Man;
    omem = strtod(argv[5], NULL);
    omen = strtod(argv[6], NULL);
    omega = omem * omen;
    direct = strtod(argv[7], NULL);
    sd_magnetic(sys,
                Ma,
                direct,
                omega);
    
    // repulsive force
    double f0;
    double deltam, delta0;
    double tau;
    f0 = 1.0;
    deltam = 1.0e-2;
    delta0 = 1.0e-2 * deltam;
    tau = 1250.0;
    sd_repulsive(sys,
                 deltam, delta0,
                 f0,
                 tau);
    
    // forces and velocities
    double * fr = calloc(n3, sizeof(double));
    double * fm = calloc(n3, sizeof(double));
    double * frx = calloc(n9, sizeof(double));
    double * fmx = calloc(n9, sizeof(double));
    double * fhn = calloc(n6, sizeof(double));
    double * fhf = calloc(n6, sizeof(double));
    double * f = calloc(n6, sizeof(double));
    double * u = calloc(n6, sizeof(double));
    double * s = calloc(n5, sizeof(double));
    double * d = calloc(n3, sizeof(double));
    double * ulub = calloc(n6, sizeof(double));
    double * uoe = calloc(n11, sizeof(double));
    
    // time step
    double delt;
    int t0, tf;
    double t;
    int ts;
    if (sys->Ma < 0.01 && sys->Ma > 0.0) {
        delt = 0.1 * sys->Ma;
    } else {
        delt = 0.001;
    }
    t0 = 0;
    tf = (int) (20.0 / delt);
    
    // output files
    char dir_name[80];
    if (Man < 1.0) {
        if (omen == 0.0) {
            sprintf(dir_name, "magnetics_np_%d_phi_0P%d_ma_%dem%d_pi_0_direct_%d_steps_%d_%d", sys->np, (int) (100.0 * sys->phi), (int) Mam, - (int) log10(Man), sys->direct, t0, tf);
        } else if (omen < 1.0) {
            sprintf(dir_name, "magnetics_np_%d_phi_0P%d_ma_%dem%d_pi_%dem%d_direct_%d_steps_%d_%d", sys->np, (int) (100.0 * sys->phi), (int) Mam, - (int) log10(Man), (int) omem, - (int) log10(omen), sys->direct, t0, tf);
        } else {
            sprintf(dir_name, "magnetics_np_%d_phi_0P%d_ma_%dem%d_pi_%de%d_direct_%d_steps_%d_%d", sys->np, (int) (100.0 * sys->phi), (int) Mam, - (int) log10(Man), (int) omem, (int) log10(omen), sys->direct, t0, tf);
        }
    }
    else {
        if (omen == 0.0) {
            sprintf(dir_name, "magnetics_np_%d_phi_0P%d_ma_%de%d_pi_0_direct_%d_steps_%d_%d", sys->np, (int) (100.0 * sys->phi), (int) Mam, (int) log10(Man), sys->direct, t0, tf);
        } else if (omen < 1.0) {
            sprintf(dir_name, "magnetics_np_%d_phi_0P%d_ma_%de%d_pi_%dem%d_direct_%d_steps_%d_%d", sys->np, (int) (100.0 * sys->phi), (int) Mam, (int) log10(Man), (int) omem, - (int) log10(omen), sys->direct, t0, tf);
        } else {
            sprintf(dir_name, "magnetics_np_%d_phi_0P%d_ma_%de%d_pi_%de%d_direct_%d_steps_%d_%d", sys->np, (int) (100.0 * sys->phi), (int) Mam, (int) log10(Man), (int) omem, (int) log10(omen), sys->direct, t0, tf);
        }
    }
    
    mkdir(dir_name, 0755);
    
    char head[100];
    
    head[0] = '\0';
    strcpy(head, "./");
    strcat(head, dir_name);
    FILE * INFO = fopen(strcat(head, "/info.dat"), "w");
    
    head[0] = '\0';
    strcpy(head, "./");
    strcat(head, dir_name);
    FILE * POS = fopen(strcat(head, "/pos.dat"), "w");
    
    head[0] = '\0';
    strcpy(head, "./");
    strcat(head, dir_name);
    FILE * FR = fopen(strcat(head, "/fr.dat"), "w");
    
    head[0] = '\0';
    strcpy(head, "./");
    strcat(head, dir_name);
    FILE * FM = fopen(strcat(head, "/fm.dat"), "w");
    
    head[0] = '\0';
    strcpy(head, "./");
    strcat(head, dir_name);
    FILE * FRX = fopen(strcat(head, "/frx.dat"), "w");
    
    head[0] = '\0';
    strcpy(head, "./");
    strcat(head, dir_name);
    FILE * FMX = fopen(strcat(head, "/fmx.dat"), "w");
    
    head[0] = '\0';
    strcpy(head, "./");
    strcat(head, dir_name);
    FILE * FHN = fopen(strcat(head, "/fhn.dat"), "w");
    
    head[0] = '\0';
    strcpy(head, "./");
    strcat(head, dir_name);
    FILE * FHF = fopen(strcat(head, "/fhf.dat"), "w");
    
    head[0] = '\0';
    strcpy(head, "./");
    strcat(head, dir_name);
    FILE * F = fopen(strcat(head, "/f.dat"), "w");
    
    head[0] = '\0';
    strcpy(head, "./");
    strcat(head, dir_name);
    FILE * STAT = fopen(strcat(head, "/stat.dat"), "w");
    
    head[0] = '\0';
    strcpy(head, "./");
    strcat(head, dir_name);
    FILE * U = fopen(strcat(head, "/u.dat"), "w");
    
    head[0] = '\0';
    strcpy(head, "./");
    strcat(head, dir_name);
    FILE * S = fopen(strcat(head, "/s.dat"), "w");
    
    head[0] = '\0';
    strcpy(head, "./");
    strcat(head, dir_name);
    FILE * ETA = fopen(strcat(head, "/eta.dat"), "w");
    
    head[0] = '\0';
    strcpy(head, "./");
    strcat(head, dir_name);
    FILE * DIPOLE = fopen(strcat(head, "/dipole.dat"), "w");
    
    head[0] = '\0';
    strcpy(head, "./");
    strcat(head, dir_name);
    FILE * ULUB = fopen(strcat(head, "/ulub.dat"), "w");
    
    head[0] = '\0';
    strcpy(head, "./");
    strcat(head, dir_name);
    FILE * UOE = fopen(strcat(head, "/uoe.dat"), "w");
    
    fprintf(INFO, "np = %d\nphi = %f\nlx = %f, ly = %f, lz = %f\n", sys->np, sys->phi, sys->lx, sys->ly, sys->lz);
    fprintf(INFO, "error = %f\nnkx = %d, nky = %d, nkz = %d\nP = %d\n", error, sys->nkx, sys->nky, sys->nkz, sys->p);
    fprintf(INFO, "Ma = %f, Omega = %f, Direction = %d", sys->Ma, sys->omega, sys->direct);
    fprintf(INFO, "gammar = %f\ndt = %f", shear_rate, delt);
    fclose(INFO);
    
    for (ts = t0; ts < tf; ts++) {
        t = (double) ts * delt;
        
        printf("\nStep %d\n", ts);
        
        // strain
        sd_strain(sys,
                  t);
        
        if (ts % 1 == 0) {
            fprintf(STAT, "%d %f %f\n", ts, sys->strain, sys->shift);
            for (i = 0; i < n3; i++) {
                fprintf(POS, "%f ", sys->pos[i]);
            }
            fprintf(POS, "\n");
            fflush(POS);
            fflush(STAT);
        }
        
        // zero arrays
        for (i = 0; i < n3; i++) {
            fr[i] = 0.0;
            fm[i] = 0.0;
            d[i] = 0.0;
        }
        for (i = 0; i < n5; i++) {
            s[i] = 0.0;
        }
        for (i = 0; i < n6; i++) {
            fhn[i] = 0.0;
            fhf[i] = 0.0;
            u[i] = 0.0;
            ulub[i] = 0.0;
            f[i] = 0.0;
        }
        for (i = 0; i < n11; i++) {
            uoe[i] = 0.0;
        }
        
       
        // update particles' positions
        update_config(sys,
                      delt,
                      ts,
                      t,
                      d,
                      fr, fm,
                      frx, fmx,
                      fhn, fhf,
                      f,
                      ulub, uoe,
                      u, s);
        
        if (ts % 1 == 0) {
            for (i = 0; i < n3; i++) {
                fprintf(FR, "%f ", fr[i]);
                fprintf(FM, "%f ", fm[i]);
                fprintf(DIPOLE, "%f ", d[i]);
            }
            fprintf(FR, "\n");
            fprintf(FM, "\n");
            fprintf(DIPOLE, "\n");
            fflush(FR);
            fflush(FM);
            fflush(DIPOLE);
            
            for (i = 0; i < n6; i++) {
                fprintf(FHN, "%f ", fhn[i]);
                fprintf(FHF, "%f ", fhf[i]);
                fprintf(F, "%f ", f[i]);
                fprintf(U, "%f ", u[i]);
                fprintf(ULUB, "%f ", ulub[i]);
            }
            fprintf(FHN, "\n");
            fprintf(FHF, "\n");
            fprintf(F, "\n");
            fprintf(U, "\n");
            fprintf(ULUB, "\n");
            fflush(FHN);
            fflush(FHF);
            fflush(F);
            fflush(U);
            fflush(ULUB);
            
            for (i = 0; i < n9; i++) {
                fprintf(FRX, "%f ", frx[i]);
                fprintf(FMX, "%f ", fmx[i]);
            }
            fprintf(FRX, "\n");
            fprintf(FMX, "\n");
            fflush(FRX);
            fflush(FMX);
     
	    for (i = 0; i < n5; i++) {
		fprintf(S, "%f ", s[i]);
	    }
	    fprintf(S, "\n");
        fflush(S);

	    fprintf(ETA, "%f\n", sys->viscosity);
        fflush(ETA);
            for (i = 0; i < sys->np; i++) {
                for (j = 0; j < 6; j++) {
                    fprintf(UOE, "%f ", uoe[i * 6 + j]);
                }
                for (j = 0; j < 5; j++) {
                    fprintf(UOE, "%f ", uoe[n6 + i * 5 + j]);
                }
            }
            fprintf(UOE, "\n");
            fflush(UOE);
        } // if ts % 10 == 0
        
    }
    
    fclose(STAT);
    fsync(STAT);
    fclose(POS);
    fsync(POS);
    fclose(FR);
    fsync(FR);
    fclose(FM);
    fsync(FM);
    fclose(FRX);
    fsync(FRX);
    fclose(FMX);
    fsync(FMX);
    fclose(F);
    fsync(F);
    fclose(FHN);
    fsync(FHN);
    fclose(FHF);
    fsync(FHF);
    fclose(U);
    fsync(U);
    fclose(S);
    fsync(S);
    fclose(ETA);
    fsync(ETA);
    fclose(FRX);
    fsync(FRX);
    fclose(FMX);
    fsync(FMX);
    fclose(DIPOLE);
    fsync(DIPOLE);
    free(ULUB);
    fsync(ULUB);
    free(UOE);
    fsync(UOE);
    
    free(pos);
    free(fr);
    free(fm);
    free(fhn);
    free(fhf);
    free(f);
    free(u);
    free(s);
    free(d);
    free(ulub);
    free(uoe);
}
