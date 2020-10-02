//
//  global.c
//  asd_trial_20190508
//
//  Created by Mingyang Tan on 2019/5/10.
//  Copyright © 2019 Mingyang Tan. All rights reserved.
//

#include "global.h"
#include <stdlib.h>
#include <stdio.h>

double *kron;
double *levi;

double *rabc;
double *x11as;
double *x12as;
double *y11as;
double *y12as;
double *y11bs;
double *y12bs;
double *x11cs;
double *x12cs;
double *y11cs;
double *y12cs;

double *rghm;
double *x11gs;
double *x12gs;
double *y11gs;
double *y12gs;
double *y11hs;
double *y12hs;
double *xms;
double *yms;
double *zms;

void glob_variables()
{
    kron = calloc(9, sizeof(double));
    levi = calloc(27, sizeof(double));
    
    kron[0] = 1.0;
    kron[4] = 1.0;
    kron[8] = 1.0;
    
    levi[ 5] = +1.0;    levi[ 7] = -1.0;
    levi[11] = -1.0;    levi[15] = +1.0;
    levi[19] = +1.0;    levi[21] = -1.0;
    
    int row, col, size;
    row = 39;
    col = 11;
    size = row * col;
    
    double *temp1 = calloc(size, sizeof(double));
    rabc = calloc(row, sizeof(double));
    x11as = calloc(row, sizeof(double));
    x12as = calloc(row, sizeof(double));
    y11as = calloc(row, sizeof(double));
    y12as = calloc(row, sizeof(double));
    y11bs = calloc(row, sizeof(double));
    y12bs = calloc(row, sizeof(double));
    x11cs = calloc(row, sizeof(double));
    x12cs = calloc(row, sizeof(double));
    y11cs = calloc(row, sizeof(double));
    y12cs = calloc(row, sizeof(double));
    
    
    FILE *datafile1;
    datafile1= fopen("data_abc.txt", "r");
    int count = 0;
    while(!feof(datafile1) && (count < size)){
        fscanf(datafile1, "%lf", &(temp1[count++]));
    }
    int i;
    for (i = 0; i < row; i++) {
        rabc[i] = temp1[i * col + 0];
        x11as[i] = temp1[i * col + 1];
        x12as[i] = temp1[i * col + 2];
        y11as[i] = temp1[i * col + 3];
        y12as[i] = temp1[i * col + 4];
        y11bs[i] = temp1[i * col + 5];
        y12bs[i] = temp1[i * col + 6];
        x11cs[i] = temp1[i * col + 7];
        x12cs[i] = temp1[i * col + 8];
        y11cs[i] = temp1[i * col + 9];
        y12cs[i] = temp1[i * col + 10];
    }
    fclose(datafile1);
    
    row = 47;
    col = 10;
    size = row * col;
    
    double *temp2 = calloc(size, sizeof(double));
    rghm = calloc(row, sizeof(double));
    x11gs = calloc(row, sizeof(double));
    x12gs = calloc(row, sizeof(double));
    y11gs = calloc(row, sizeof(double));
    y12gs = calloc(row, sizeof(double));
    y11hs = calloc(row, sizeof(double));
    y12hs = calloc(row, sizeof(double));
    xms = calloc(row, sizeof(double));
    yms = calloc(row, sizeof(double));
    zms = calloc(row, sizeof(double));
    
    FILE *datafile2;
    datafile2= fopen("data_ghm.txt", "r");
    count = 0;
    while(!feof(datafile2) && (count < size)){
        fscanf(datafile2, "%lf", &(temp2[count++]));
    }
    for (i = 0; i < row; i++) {
        rghm[i] = temp2[i * col + 0];
        x11gs[i] = temp2[i * col + 1];
        x12gs[i] = temp2[i * col + 2];
        y11gs[i] = temp2[i * col + 3];
        y12gs[i] = temp2[i * col + 4];
        y11hs[i] = temp2[i * col + 5];
        y12hs[i] = temp2[i * col + 6];
        xms[i] = temp2[i * col + 7];
        yms[i] = temp2[i * col + 8];
        zms[i] = temp2[i * col + 9];
    }
    fclose(datafile2);
    
    free(temp1);
    free(temp2);
    
}

