// anchor.h
//
// Last Modified: 13, Jul 2007
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

#ifndef __ANCHOR_H__
#define __ANCHOR_H__

#include "prrn.h"
#include <vector>

class Alignment;
class PhylogeneticTree;

void eraseConservedBranch(const Alignment&, const PhylogeneticTree&, std::vector<int>&, prrn::SCORE, bool);

enum anchor_mtd {NoAnchor, NonGap, Conservation};

void detectNonGapRegion(const Alignment&, std::vector<std::pair<size_t, size_t> >&);
void detectConservedRegion(const Alignment&, std::vector<std::pair<size_t, size_t> >&);
void calcAnchorPoint(const Alignment&, std::vector<std::pair<size_t, size_t> >&, anchor_mtd);

// UCR: UnChangedRegion, UCS: UnChangedSubmamily, UCRS: both UCR and UCS
enum ucf_strategy {UCR, UCS, UCRS, NoUC};
void detectUnchangedRegion(const Alignment&, const Alignment&, std::vector<std::pair<size_t, size_t> >&);
void eraseUnchangedBranch(const Alignment&, const Alignment&, const PhylogeneticTree&, std::vector<int>&);

void calcNonConservedRegion(size_t, const std::vector<std::pair<size_t, size_t> >&,
	std::vector<std::pair<size_t, size_t> >&);
void decomposeAlignment(const Alignment&, const std::vector<std::pair<size_t, size_t> >&,
	std::vector<Alignment>&);
void replaceAlignment(const std::vector<std::pair<size_t, size_t> >&, const std::vector<Alignment>&,
	Alignment&);

#endif
