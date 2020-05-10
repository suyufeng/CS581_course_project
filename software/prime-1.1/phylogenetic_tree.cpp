// phylogenetic_tree.cpp
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

#include "phylogenetic_tree.h"
#include <iostream>
#include <iomanip>
#ifdef PT_VERBOSE
#include <iterator>
#endif
#include <cmath>

using namespace std;
using namespace prrn;

PhylogeneticTree::PhylogeneticTree(const Matrix<size_t>& m)
{
    construct(m);
}

namespace
{
    inline double correlate(double a, double b)
    {
	return a * b / (a + b);
    }
}

void PhylogeneticTree::construct(const Matrix<size_t>& m)
{
#ifdef DEBUG
    assert(isSymmetric(m));
#endif
    node_.clear();
    nelm_ = m.getRowSize();
    node_.reserve(nelm_-1);
    weight_.assign(getNumEdges(), 1.0);

    PTAlg pta(m);

    vector<size_t> idxs;
    idxs.reserve(nelm_);
    for(size_t i = 0; i < nelm_; ++i)
	idxs.push_back(i);

#ifdef PT_VERBOSE
    copy(idxs.begin(), idxs.end(), ostream_iterator<size_t>(cout, " "));
    pta.show();
#endif

    const size_t nitr = nelm_ - 2;
    for(size_t i = 0; i < nitr; ++i)
    {
	double minLength = limitDoubleMax;
	size_t minJ = 0;
	size_t minK = 0;

	pta.selectOTU(minLength, minJ, minK, idxs, nitr-i+1);

	Node n;
	const double normalizer = 1.0 / ((nitr-i)<<1);
	n.child_left = minJ;
	if(minJ < nelm_)
	    n.corr_child_left = pta.length(minJ, minLength) * normalizer;
	else
	{
	    size_t idx = minJ - nelm_;
	    node_[idx].parent = i + nelm_;
	    node_[idx].len_parent = pta.length(minJ, minLength) * normalizer;

	    n.corr_child_left = node_[idx].len_parent
		+ correlate(node_[idx].corr_child_left, node_[idx].corr_child_right);
	}
	n.child_right = minK;
	if(minK < nelm_)
	    n.corr_child_right = pta.length(minK, minLength) * normalizer;
	else
	{
	    size_t idx = minK - nelm_;
	    node_[idx].parent = i + nelm_;
	    node_[idx].len_parent = pta.length(minK, minLength) * normalizer;

	    n.corr_child_right = node_[idx].len_parent
		+ correlate(node_[idx].corr_child_left, node_[idx].corr_child_right);
	}
	node_.push_back(n);

	pta.update(minJ, minK, idxs);
	idxs.erase(remove(idxs.begin(), idxs.end(), minK), idxs.end());
	*remove(idxs.begin(), idxs.end(), minJ) = i + nelm_;
#ifdef PT_VERBOSE
	copy(idxs.begin(), idxs.end(), ostream_iterator<size_t>(cout, " "));
	pta.show();
#endif
    }
#ifdef DEBUG
    assert(idxs.size() == 2);
    assert(idxs[0] < idxs[1]);
#endif

    // connection of the two rest OTUs
    Node n;
    n.child_left = idxs[0];
    n.child_right = idxs[1];

    // the followings use the condition n.child_left < n.child_right
    const size_t ridx = n.child_right - nelm_;
    node_[ridx].parent = nitr + nelm_;
    node_[ridx].len_parent = pta.distance(n.child_left, n.child_right);
    n.corr_child_right = node_[ridx].len_parent
	+ correlate(node_[ridx].corr_child_left, node_[ridx].corr_child_right);
    // n.child_left can be a leaf
    if(n.child_left < nelm_)
	n.corr_child_left = node_[ridx].len_parent;
    else
    {
	const size_t lidx = n.child_left - nelm_;
	node_[lidx].parent = node_[ridx].parent;
	node_[lidx].len_parent = node_[ridx].len_parent;

	n.corr_child_left = node_[lidx].len_parent
	    + correlate(node_[lidx].corr_child_left, node_[lidx].corr_child_right);
    }
    node_.push_back(n);

    calc_parent_corr_len();
}

void PhylogeneticTree::calculateWeights(double ef)
{
#ifdef DEBUG
    assert(ef >= 0.0);
#endif
    const size_t esize = getNumEdges();
    weight_.assign(esize, 1.0);

    const size_t dnsize = node_.size()-1;
    const double F = 1.0 / ef;
    for(size_t i = 0; i < dnsize; ++i)
    {
	const double prdct_left = node_[i].corr_child_right * node_[i].corr_parent;
	const double prdct_right = node_[i].corr_child_left * node_[i].corr_parent;
	const double prdct_parent = node_[i].corr_child_left * node_[i].corr_child_right;

	const double FS = F / (prdct_left + prdct_right + prdct_parent);

	const double add_rp = node_[i].corr_child_right + node_[i].corr_parent;
	const double add_pl = node_[i].corr_child_left + node_[i].corr_parent;
	const double add_lr = node_[i].corr_child_left + node_[i].corr_child_right;

	weight_[node_[i].child_left]
	    *= sqrt(FS * prdct_left * add_pl * add_lr / (node_[i].corr_child_left * add_rp));

	weight_[node_[i].child_right]
	    *= sqrt(FS * prdct_right * add_rp * add_lr / (node_[i].corr_child_right * add_pl));

	if(i+nelm_ >= esize)
	{
	    if(node_[dnsize].child_left == i+nelm_)
		weight_[node_[dnsize].child_right]
		    *= sqrt(FS * prdct_parent * add_pl * add_rp / (node_[i].corr_parent * add_lr));
	    else
		weight_[node_[dnsize].child_left]
		    *= sqrt(FS * prdct_parent * add_pl * add_rp / (node_[i].corr_parent * add_lr));
	}
	else
	    weight_[i+nelm_]
		*= sqrt(FS * prdct_parent * add_pl * add_rp / (node_[i].corr_parent * add_lr));
#ifdef DEBUG
	cout << "1.0/FS = " << setw(8) << FS << "\n";
	cout << "(" << setw(3) << node_[i].child_left << ", ";
	cout << setw(3) << node_[i].child_right << ", ";
	cout << setw(3) << i+nelm_ << ") =  ";

	cout << "(" << setw(8) << weight_[node_[i].child_left] << ", ";
	cout << setw(8) << weight_[node_[i].child_right] << ", ";
	if(i+nelm_ >= esize && node_[dnsize].child_left == i+nelm_)
	    cout << setw(8) << weight_[node_[dnsize].child_right] << ")\n";
	else if(i+nelm_ >= esize && node_[dnsize].child_right == i+nelm_)
	    cout << setw(8) << weight_[node_[dnsize].child_left] << ")\n";
	else
	    cout << setw(8) << weight_[i+nelm_] << ")\n";
	cout << "\n";
#endif
    }
    calc_pair_weights();
}

size_t PhylogeneticTree::getNumLeaves() const
{
    return nelm_;
}

size_t PhylogeneticTree::getNumEdges() const
{
    return nelm_ + nelm_ - 3;
}

size_t PhylogeneticTree::getNumNeighbors() const
{
    return nelm_ - 1;
}

size_t PhylogeneticTree::getParent(size_t n) const
{
    if(n >= nelm_)
	return node_[n-nelm_].parent;
    for(size_t i = 0; i < nelm_; ++i)
	if(node_[i].child_left == n || node_[i].child_right == n)
	    return i+nelm_;
    return limitSizeMax;
}

pair<size_t, size_t> PhylogeneticTree::getNeighbor(size_t i) const
{
#ifdef DEBUG
    return make_pair(node_.at(i).child_left, node_.at(i).child_right);
#else
    return make_pair(node_[i].child_left, node_[i].child_right);
#endif
}

pair<double, double> PhylogeneticTree::getNeighborWeights(size_t i) const
{
#ifdef DEBUG
    return make_pair(weight_.at(node_.at(i).child_left),
	    weight_.at(node_.at(i).child_right));
#else
    return make_pair(weight_[node_[i].child_left],
	    weight_[node_[i].child_right]);
#endif
}

double PhylogeneticTree::getEdgeWeight(size_t i) const
{
#ifdef DEBUG
    return weight_.at(i);
#else
    return weight_[i];
#endif
}

double PhylogeneticTree::getPairWeight(size_t i, size_t j) const
{
#ifdef DEBUG
    assert(i < nelm_);
    assert(j < nelm_);
    assert(i != j);
    assert(!weight_.empty());
#endif
    return pairweight_(i, j);
}

void PhylogeneticTree::getTreeIndex(size_t i, vector<size_t>& vidx) const
{
    if(node_[i].child_left < nelm_)
	vidx.push_back(node_[i].child_left);
    else
	getTreeIndex(node_[i].child_left-nelm_, vidx);
    if(node_[i].child_right < nelm_)
	vidx.push_back(node_[i].child_right);
    else
	getTreeIndex(node_[i].child_right-nelm_, vidx);
}

void PhylogeneticTree::show() const
{
    static const size_t offset = 3;
    static const size_t offset_len = 8;

    cout << "no. of    leaves:\t" << setw(3) << nelm_ << "\n";
    cout << "no. of     edges:\t" << setw(3) << getNumEdges() << "\n";
    cout << "no. of neighbors:\t" << setw(3) << getNumNeighbors() << "\n";

    const size_t nsize = node_.size()-1;
    for(size_t i = 0; i < nsize; ++i)
    {
	cout << "Node " << setw(offset) << i+nelm_ << " = ";

	if(node_[i].child_left >= nelm_)
	    cout << "Node " << setw(offset) << node_[i].child_left
		<< " (" << setw(offset_len) << node_[node_[i].child_left-nelm_].len_parent << ") + ";
	else
	    cout << "Leaf " << setw(offset) << node_[i].child_left
		<< " (" << setw(offset_len) << node_[i].corr_child_left << ") + ";

	if(node_[i].child_right >= nelm_)
	    cout << "Node " << setw(offset) << node_[i].child_right
		<< " (" << setw(offset_len) << node_[node_[i].child_right-nelm_].len_parent << ")\n";
	else
	    cout << "Leaf " << setw(offset) << node_[i].child_right
		<< " (" << setw(offset_len) << node_[i].corr_child_right << ")\n";
    }
    cout << "FinEdge  = ";
    // child_left can be less than nelm_
    if(node_[nsize].child_left >= nelm_)
	cout << "Node " << setw(offset) << node_[nsize].child_left
	    << " (" << setw(offset_len) << node_[node_[nsize].child_left-nelm_].len_parent << ") + ";
    else
	cout << "Leaf " << setw(offset) << node_[nsize].child_left
	    << " (" << setw(offset_len) << node_[nsize].corr_child_left << ") + ";

    // child_right must be more than or equal to nelm_
    cout << "Node " << setw(offset) << node_[nsize].child_right
	<< " (" << setw(offset_len) << node_[node_[nsize].child_right-nelm_].len_parent << ")\n";
    cout << "\n";
}

void PhylogeneticTree::showEdge() const
{
    const size_t esize = getNumEdges();
    for(size_t i = 0; i <= esize; ++i)
    {
	size_t idx = i >> 1;
	if(i & 1)
	{
	    cout << "(" << setw(3) << node_[idx].child_right << ", "
		<< setw(3) << idx+nelm_ << "; ";
	    if(node_[idx].child_right < nelm_)
		cout << setw(8) << node_[idx].corr_child_right << ")" << "\n";
	    else
		cout << setw(8) << node_[node_[idx].child_right-nelm_].len_parent << ")" << "\n";
	}
	else
	{
	    cout << "(" << setw(3) << node_[idx].child_left << ", "
		<< setw(3) << idx+nelm_ << "; ";
	    if(node_[idx].child_left < nelm_)
		cout << setw(8) << node_[idx].corr_child_left << ")" << "\n";
	    else
		cout << setw(8) << node_[node_[idx].child_left-nelm_].len_parent << ")" << "\n";
	}
    }
}

void PhylogeneticTree::showNode() const
{
    cout << "Node num {pid: parentln, cli, cri,   clcorr,   crcorr,  parcorr}\n";
    const size_t nsize = node_.size();
    for(size_t i = 0; i < nsize; ++i)
    {
	cout << "Node " << setw(3) << i+nelm_ << " {";
	if(node_[i].parent == limitSizeMax)
	    cout << "nul";
	else
	    cout << setw(3) << node_[i].parent;
	cout << ": " << setw(8) << node_[i].len_parent
	    << ", " << setw(3) << node_[i].child_left
	    << ", " << setw(3) << node_[i].child_right
	    << ", " << setw(8) << node_[i].corr_child_left
	    << ", " << setw(8) << node_[i].corr_child_right
	    << "; " << setw(8) << node_[i].corr_parent << "}\n";
    }
}

void PhylogeneticTree::showWeight() const
{
    for(size_t i = 0; i < weight_.size(); ++i)
	cout << setw(3) << i << ":\t" << weight_[i] << "\n";
}

void PhylogeneticTree::showPairWeight() const
{
    for(size_t i = 0; i < nelm_-1; ++i)
	for(size_t j = i+1; j < nelm_; ++j)
	    cout << "(" << setw(3) << i << ", " << setw(3) << j << "):\t" << getPairWeight(i, j) << "\n";
}

void PhylogeneticTree::calc_parent_corr_len()
{
    const int nsize = node_.size();
    for(int i = nsize-2; i >= 0; --i)
    {
	int idx = node_[i].parent - nelm_;
	if(node_[idx].child_left-nelm_ == static_cast<size_t>(i))
	{
	    if(idx == nsize-1)
	    {
		node_[i].corr_parent = node_[idx].corr_child_right;
		continue;
	    }
	    node_[i].corr_parent = node_[i].len_parent;
	    node_[i].corr_parent += correlate(node_[idx].corr_child_right,
		    node_[idx].corr_parent);
	}
	else if(node_[idx].child_right-nelm_ == static_cast<size_t>(i))
	{
	    if(idx == nsize-1)
	    {
		node_[i].corr_parent = node_[idx].corr_child_left;
		continue;
	    }
	    node_[i].corr_parent = node_[i].len_parent;
	    node_[i].corr_parent += correlate(node_[idx].corr_child_left,
		    node_[idx].corr_parent);
	}
	else
	{
	    cerr << "calc_parent_corr_len: unknown error\n";
	    abort();
	}
    }
}

void PhylogeneticTree::calc_pair_weights()
{
    pairweight_.resize(nelm_, nelm_);
    for(size_t i = 0; i < nelm_; ++i)
	pairweight_(i, i) = 0.0;

    const size_t nlidx = node_.size() - 1;
    for(size_t i = 0; i <= nlidx; ++i)
    {
	size_t cl = node_[i].child_left;
	size_t cr = node_[i].child_right;
	if(cl < nelm_)
	{
	    calc_downstream_pair_weight(cl, cr, weight_[cl]);
	    calc_upstream_pair_weight(cl, i+nelm_, weight_[cl]);
	}
	if(cr < nelm_)
	{
	    calc_downstream_pair_weight(cr, cl, weight_[cr]);
	    calc_upstream_pair_weight(cr, i+nelm_, weight_[cr]);
	}
    }
#ifdef DEBUG
    for(size_t i = 0; i < nelm_-1; ++i)
	for(size_t j = i+1; j < nelm_; ++j)
	    assert(isApproxEqual(pairweight_(i, j), pairweight_(j, i)));
#endif

}

void PhylogeneticTree::calc_downstream_pair_weight(size_t source, size_t target, double weight)
{
    if(target < getNumEdges())
	weight *= weight_[target];
    if(target < nelm_)
	pairweight_(source, target) = weight;
    else
    {
	size_t idx = target-nelm_;
	calc_downstream_pair_weight(source, node_[idx].child_left, weight);
	calc_downstream_pair_weight(source, node_[idx].child_right, weight);
    }
}

void PhylogeneticTree::calc_upstream_pair_weight(size_t source, size_t current, double weight)
{
    const size_t lidx = getNumEdges() - 1;
    size_t parent = node_[current-nelm_].parent;
    if(parent == limitSizeMax)
	return;
    while(1)
    {
	size_t pidx = parent - nelm_;
	size_t target = (node_[pidx].child_left == current) ? node_[pidx].child_right : node_[pidx].child_left;

	if(target < lidx)
	    calc_downstream_pair_weight(source, target, weight * weight_[current]);
	else
	{
	    if(target < nelm_)
		pairweight_(source, target) = weight * weight_[target];
	    else
	    {
		size_t tidx = target - nelm_;
		double tw = weight * weight_[lidx];
		calc_downstream_pair_weight(source, node_[tidx].child_left, tw);
		calc_downstream_pair_weight(source, node_[tidx].child_right, tw);
	    }
	}

	weight *= weight_[(current-nelm_ < lidx) ? current : lidx];
	current = parent;
	if((parent = node_[pidx].parent) == limitSizeMax)
	    break;
    }
}
