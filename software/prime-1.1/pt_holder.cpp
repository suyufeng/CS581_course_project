// pt_holder.cpp
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

#include "pt_holder.h"
#include "prrn.h"
#include "sequence.h"
#include "phylogenetic_tree.h"
#include <iostream>
#include <iomanip>

using namespace std;
using namespace prrn;

PTHolder::PTHolder(const vector<Sequence>& vs,
	const Matrix<size_t>& dmat, double eqfactor)
    : vs_(vs)
    , eqfactor_(eqfactor)
    , pt_(dmat)
    , nelm_(vs.size())
    , dbn_(limitSizeMax)
{
    if(isWeightMode())
	pt_.calculateWeights(eqfactor);

    for(size_t i = 0; i < nelm_; ++i)
	seqidx_.insert(make_pair(vs[i].getName(), i));
    partial_weight_.reserve(nelm_);
    partial_weight_.assign(nelm_, 1.0);
}

void PTHolder::construct(const Matrix<size_t>& dmat)
{
    dbn_ = limitSizeMax;
    pt_.construct(dmat);
    if(isWeightMode())
	pt_.calculateWeights(eqfactor_);
    fill(partial_weight_.begin(), partial_weight_.end(), 1.0);
}

void PTHolder::setDivideBranch(size_t branch) const
{
    dbn_ = branch;
    fill(partial_weight_.begin(), partial_weight_.end(), 1.0);
    if(!isWeightMode())
	return;

    if(branch >= nelm_)
	setDSPartialWeight(branch, pt_.getEdgeWeight(branch));
    else
	partial_weight_[branch] = pt_.getEdgeWeight(branch);


    const size_t lastedge = pt_.getNumEdges() - 1;
    double weight = 1.0;
    size_t parent_branch = pt_.getParent(branch);
    while(1)
    {
	pair<size_t, size_t> neighbor = pt_.getNeighbor(parent_branch-nelm_);
	size_t target = (neighbor.first == branch) ? neighbor.second : neighbor.first;
	if(target < nelm_)
	{
	    if(pt_.getParent(parent_branch) == limitSizeMax)
		partial_weight_[target] = weight;
	    else
		partial_weight_[target] = weight * pt_.getEdgeWeight(target);
	}
	else if(target < lastedge)
	    setDSPartialWeight(target, weight * pt_.getEdgeWeight(target));
	else
	    setDSPartialWeight(target, weight);

	branch = parent_branch;
	weight *= pt_.getEdgeWeight((branch < lastedge) ? branch : lastedge);
	if((parent_branch = pt_.getParent(branch)) == limitSizeMax)
	    break;
    }
}

size_t PTHolder::getDivideBranch() const
{
    return dbn_;
}

const vector<Sequence>& PTHolder::getSequenceVector() const
{
    return vs_;
}

const PhylogeneticTree& PTHolder::getPhylogeneticTree() const
{
    return pt_;
}

size_t PTHolder::getIndex(const string& seqname) const
{
    return seqidx_[seqname];
}

bool PTHolder::isWeightMode() const
{
    return (eqfactor_ > 0.0);
}

double PTHolder::getDivideEdgeWeight() const
{
    return (isWeightMode()) ? pt_.getEdgeWeight(dbn_) : 1.0;
}

double PTHolder::getEdgeWeight(const string& name) const
{
    return (isWeightMode()) ? pt_.getEdgeWeight(seqidx_[name]) : 1.0;
}

double PTHolder::getPartialWeight(const string& name) const
{
    return (isWeightMode()) ? partial_weight_[seqidx_[name]] : 1.0;
}

double PTHolder::getPairWeight(size_t a, size_t b) const
{
    return (isWeightMode()) ? pt_.getPairWeight(a, b) : 1.0;
}

double PTHolder::getPairWeight(const string& a, const string& b) const
{
    return (isWeightMode()) ? pt_.getPairWeight(seqidx_[a], seqidx_[b]) : 1.0;
}

void PTHolder::show() const
{
    pt_.show();
}

void PTHolder::showEdge() const
{
    pt_.showEdge();
}

void PTHolder::showNode() const
{
    pt_.showNode();
}

void PTHolder::showWeight() const
{
    if(isWeightMode())
	pt_.showWeight();
}

void PTHolder::showPartialWeight() const
{
    if(!isWeightMode())
	return;
    for(size_t i = 0; i < nelm_; ++i)
	cout << setw(3) << i << ":\t" << partial_weight_[i] << "\n";
}

void PTHolder::showPairWeight() const
{
    if(isWeightMode())
	pt_.showPairWeight();
}

void PTHolder::setDSPartialWeight(size_t nbranch, double weight) const
// DS means DownStream
{
#ifdef DEBUG
    if(nbranch < nelm_)
    {
	partial_weight_[nbranch] = weight * pt_.getEdgeWeight(nbranch);
	abort();
    }
    else
    {
#endif
	pair<size_t, size_t> nidx = pt_.getNeighbor(nbranch-nelm_);

	if(nidx.first < nelm_)
	    partial_weight_[nidx.first] = weight * pt_.getEdgeWeight(nidx.first);
	else
	    setDSPartialWeight(nidx.first, weight * pt_.getEdgeWeight(nidx.first));

	if(nidx.second < nelm_)
	    partial_weight_[nidx.second] = weight * pt_.getEdgeWeight(nidx.second);
	else
	    setDSPartialWeight(nidx.second, weight * pt_.getEdgeWeight(nidx.second));
#ifdef DEBUG
    }
#endif
}
