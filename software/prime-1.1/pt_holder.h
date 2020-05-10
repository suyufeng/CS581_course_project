// pt_holder.h
//
// Last Modified: 9, May 2007
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

#ifndef __PT_HOLDER_H__
#define __PT_HOLDER_H__

#include "phylogenetic_tree.h"
#include <string>
#include <map>

class PTHolder
{
    public:
	PTHolder();	// do not define
	PTHolder(const std::vector<Sequence>&,
		const Matrix<size_t>&, double = 0.0);

	~PTHolder(){}

	void construct(const Matrix<size_t>&);

	void setDivideBranch(size_t) const;

	size_t getDivideBranch() const;
	const std::vector<Sequence>& getSequenceVector() const;
	const PhylogeneticTree& getPhylogeneticTree() const;
	size_t getIndex(const std::string&) const;

	bool isWeightMode() const;

	double getDivideEdgeWeight() const;
	double getEdgeWeight(const std::string&) const;
	double getPartialWeight(const std::string&) const;
	double getPairWeight(size_t, size_t) const;
	double getPairWeight(const std::string&, const std::string&) const;

	void show() const;
	void showEdge() const;
	void showNode() const;
	void showWeight() const;
	void showPartialWeight() const;
	void showPairWeight() const;
	void showTree() const;

    private:	
	void setDSPartialWeight(size_t, double) const;

	const std::vector<Sequence>& vs_;
	const double eqfactor_;

	mutable std::map<std::string, size_t> seqidx_;
	PhylogeneticTree pt_;
	const size_t nelm_;
	mutable size_t dbn_;	// divided branch number
	mutable std::vector<double> partial_weight_;
};

#endif
