// util.h
//
// Last Modified: 8, Mar 2007
//
// Copyright (c) 2004-2007 Shinsuke Yamada
//
// This file is part of PRIME.
//
// PRIME is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 2 of the License, or
// (at your option) any later version.
//
// PRIME is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with PRIME; if not, write to the Free Software
// Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA

#ifndef __UTIL_H__
#define __UTIL_H__

#include <iosfwd>
#include<cmath>
#ifndef CLIMITS
#include <limits>
#else
#include <cfloat>
#endif

double roundDbl(double);

int iexp(int, size_t);

inline bool isApproxEqual(double a, double b)
{
#ifndef CLIMITS
    return (fabs(a - b) < std::numeric_limits<double>::epsilon());
#else
    return (fabs(a - b) < DBL_EPSILON);
#endif
}

#endif
