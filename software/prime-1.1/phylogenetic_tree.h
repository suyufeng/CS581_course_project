// phylogenetic_tree.h
//
// Last Modified: 11, May 2007
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

#ifndef __PHYLOGENETIC_TREE_H__
#define __PHYLOGENETIC_TREE_H__

#include "prrn.h"
#include "pt_alg.h"

class PhylogeneticTree
{
    public:
	PhylogeneticTree(){}
	PhylogeneticTree(const Matrix<size_t>&);
	~PhylogeneticTree(){}

	void construct(const Matrix<size_t>&);

	void calculateWeights(double);

	size_t getNumLeaves() const;
	size_t getNumEdges() const;
	size_t getNumNeighbors() const;
	size_t getParent(size_t) const;
	std::pair<size_t, size_t> getNeighbor(size_t) const;
	std::pair<double, double> getNeighborWeights(size_t) const;
	double getEdgeWeight(size_t) const;
	double getPairWeight(size_t, size_t) const;

	void getTreeIndex(size_t, std::vector<size_t>&) const;

	void show() const;
	void showEdge() const;
	void showNode() const;
	void showWeight() const;
	void showPairWeight() const;

    private:
	void   calc_parent_corr_len();
	void calc_pair_weights();
	void calc_downstream_pair_weight(size_t, size_t, double);
	void calc_upstream_pair_weight(size_t, size_t, double);

#ifndef UPGMA
	typedef nj::PTAlg PTAlg;
#else
	typedef upgma::PTAlg PTAlg;
#endif
	typedef struct _Node
	{
	    size_t parent;
	    double len_parent;

	    // corr_child_{left|right} corresponds to branch length
	    // if child_{left|right} is a leaf
	    size_t child_left;
	    size_t child_right;
	    double corr_child_left;
	    double corr_child_right;

	    // corr_parent does NOT correspond to branch length
	    double corr_parent;

	    _Node()
		: parent(prrn::limitSizeMax)
		, len_parent(0.0)
		, child_left(prrn::limitSizeMax)
		, child_right(prrn::limitSizeMax)
		, corr_child_left(0.0)
		, corr_child_right(0.0)
		, corr_parent(0.0)
		{}
	    ~_Node(){}
	} Node;

	std::vector<Node> node_;
	std::vector<double> weight_;
	Matrix<double> pairweight_;
	size_t nelm_;
};

#endif
