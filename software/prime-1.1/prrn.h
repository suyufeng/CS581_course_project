// prrn.h
//
// Last Modified: 20, Mar 2007
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

#ifndef __PRRN_H__
#define __PRRN_H__

#include <cstddef>
#ifndef CLIMITS
#include <limits>
#else
#include <climits>
#include <cfloat>
#endif

namespace prrn
{
    typedef double SCORE;
#ifndef CLIMITS
    const char limitCharMax = std::numeric_limits<char>::max();
    const size_t limitSizeMax = std::numeric_limits<size_t>::max();
    const double limitDoubleMax = std::numeric_limits<double>::max();
    const SCORE limitScoreMax = std::numeric_limits<SCORE>::max();
#else
    const char limitCharMax = CHAR_MAX;
    const size_t limitSizeMax = UINT_MAX;
    const double limitDoubleMax = DBL_MAX;
    const SCORE limitScoreMax = DBL_MAX;
#endif
    const size_t limitSizeZero = 0;
    const double limitDoubleZero = 0.0;
    const SCORE limitScoreZero = 0.0;

    enum psa_alg {SeqAffine, SeqLong};

    enum gsa_alg {ProfAffine, ProfLong};

    enum dist_calc {Psa, Oligo};
}

#endif
