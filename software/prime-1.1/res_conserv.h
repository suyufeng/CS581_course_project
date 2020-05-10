// res_conserv.h
//
// Last Modified: 11, Jul 2007
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

#ifndef __RES_CONSERV_H__
#define __RES_CONSERV_H__

#include <vector>
#include <iosfwd>

class Alignment;
class Sequence;

void calcHenikoffWeight(const Alignment&, std::vector<double>&);

void calcBackground(const std::vector<Sequence>&, size_t);

void calcBackground(const Alignment&, size_t);

void calcREconservation(const Alignment&, const std::vector<double>&, std::vector<double>&);

void showResidueConservation(const Alignment&, const std::vector<double>& rescns, std::ostream&);

#endif
