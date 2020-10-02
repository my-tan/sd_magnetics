//
//  magnetics.h
//  asd_magnetics_20200622
//
//  Created by Mingyang Tan on 6/22/20.
//  Copyright © 2020 Mingyang Tan. All rights reserved.
//

#ifndef magnetics_h
#define magnetics_h

#include <stdio.h>
#include "sd.h"
#include "ewald.h"

void magnetics(struct sd * sys,
               struct ewald * ewald_sys,
               const double t,
               double * d,
               double * fm);

#endif /* magnetics_h */
