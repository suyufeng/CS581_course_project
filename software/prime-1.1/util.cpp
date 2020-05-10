// util.cpp
//
// Last Modified: 15, Dec 2006
//
// Copyright (c) 2004-2006 Shinsuke Yamada
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

#include "util.h"
#include <cmath>

using namespace std;

double roundDbl(double val)
{
    const double upper = ceil(val);
    const double lower = floor(val);
    return (upper-val <= val-lower) ? upper : lower;
}

int iexp(int x, size_t n)
{
    int y = 1, p = x;
    while(1)
    {
	if(n & 1)
	    y *= p;
	n >>= 1;
	if(n == 0)
	    return y;
	p *= p;
    }
}
