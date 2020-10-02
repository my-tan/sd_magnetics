//
//  cell.h
//  asd_sphere_20200805
//
//  Created by Mingyang Tan on 8/5/20.
//  Copyright © 2020 Mingyang Tan. All rights reserved.
//

#ifndef cell_h
#define cell_h

#include <stdio.h>

void pos_fraction(const double lx, const double ly, const double lz,
                  const double xlh, const double ylh, const double zlh,
                  const double strain,
                  const double * pos,
                  double * fpos);

void return_to_cell(const double lx, const double ly, const double lz,
                    const double xlh, const double xrh,
                    const double ylh, const double yrh,
                    const double zlh, const double zrh,
                    const double strain,
                    double * pos);

void min_image(const double lx, const double ly, const double lz,
               const double xlh, const double xrh,
               const double ylh, const double yrh,
               const double zlh, const double zrh,
               const double strain,
               double * dist);

#endif /* cell_h */
