//
//  mynla.h
//  asd_sphere_20200508
//
//  Created by Mingyang Tan on 5/10/20.
//  Copyright © 2020 Mingyang Tan. All rights reserved.
//

#ifndef mynla_h
#define mynla_h

#include <stdio.h>

void cholesky_inverse(const int m,
                      const double * A, const int lda,
                      double * B, const int ldb);

void lu_decomp(const int m,
               const double * A, const int lda,
               double * L, const int ldl,
               double * U, const int ldu);

#endif /* mynla_h */
