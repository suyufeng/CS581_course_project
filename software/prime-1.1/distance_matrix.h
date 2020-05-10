// distance_matrix.h
//
// Last Modified: 14, Jul 2007
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

#ifndef __DISTANCE_MATRIX_H__
#define __DISTANCE_MATRIX_H__

#include "matrix.hpp"
#include "sequence.h"

class Alignment;
class PairwiseDPAlgorithm;

void calculateDistanceMatrix(const std::vector<Sequence>&, const PairwiseDPAlgorithm&,
	Matrix<size_t>&, bool = false);

void calculateDistanceMatrix(const Alignment&, Matrix<size_t>&, bool = false);

void calculateDistanceMatrix(const std::vector<Sequence>&, Matrix<size_t>&, bool = false);

#endif
