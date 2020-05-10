// pt_alg.cpp
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

#include "pt_alg.h"
#include "prrn.h"
#include "util_std.hpp"
#include <iostream>
#include <iomanip>

using namespace std;
using namespace prrn;

class Sequence;

namespace nj
{
    PTAlg::iterator PTAlg::get_iter(size_t i, size_t j)
    {
#ifdef DEBUG
	assert(i < j);
	return find_if(dm_.at(i)->begin(), dm_.at(i)->end(),
		bind2nd(equal_to_1st<size_t, double>(), j));
#else
	return find_if(dm_[i]->begin(), dm_[i]->end(),
		bind2nd(equal_to_1st<size_t, double>(), j));
#endif
    }

    PTAlg::const_iterator PTAlg::get_iter(size_t i, size_t j) const
    {
#ifdef DEBUG
	assert(i < j);
	return find_if(dm_.at(i)->begin(), dm_.at(i)->end(),
		bind2nd(equal_to_1st<size_t, double>(), j));
#else
	return find_if(dm_[i]->begin(), dm_[i]->end(),
		bind2nd(equal_to_1st<size_t, double>(), j));
#endif
    }

    PTAlg::PTAlg(const Matrix<size_t>& m)
	: nelm_(m.getColumnSize())
	, insidx_(nelm_)
	, dm_((nelm_<<1) - 3)
	, totdist_((nelm_<<1) - 2)
	{
	    size_t i = 0;
	    size_t ressize = nelm_-1;
	    for(; i < nelm_-1; ++i)
	    {
		dm_[i] = new vector<elm_type>;
		dm_[i]->reserve(ressize--);
		for(size_t j = i+1; j < nelm_; ++j)
		    dm_[i]->push_back(elm_type(j, m(i, j)));
	    }
	    const size_t size = (nelm_<<1) - 3;
	    ressize = nelm_-1;
	    for(; i < size; ++i)
	    {
		dm_[i] = new vector<elm_type>;
		dm_[i]->reserve(ressize--);
	    }

	    for(i = 0; i < nelm_; ++i)
	    {
		for(size_t p = 0, q = i-1; p < i; ++p, --q)
		    totdist_[i] += (*dm_[p])[q].second;

		const size_t size = dm_[i]->size();
		for(size_t j = 0; j < size; ++j)
		    totdist_[i] += (*dm_[i])[j].second;
	    }
	    for(; i <= size; ++i)
		totdist_[i] = 0.0;
	}

    PTAlg::~PTAlg()
    {
	for_each(dm_.begin(), dm_.end(), DeleteObject());
    }

    size_t PTAlg::size() const
    {
	return nelm_;
    }

    void PTAlg::selectOTU(double& minlen, size_t& min_j, size_t& min_k,
	    const vector<size_t>& idxs, size_t idsize) const
    {
	// 'multiplier' is equal to 'n-2' where 'n' is the number of OTUs
	const size_t multiplier = idsize-1;
	for(size_t j = 0; j < idsize; ++j)
	{
	    for(size_t k = j+1; k <= idsize; ++k)
	    {
		double tmp = multiplier * distance(idxs[j], idxs[k]);
		tmp -= totdist_[idxs[j]];
		tmp -= totdist_[idxs[k]];
#ifdef PT_VERBOSE
cout << "(" << idxs[j] << ", " << idxs[k] << "; " << tmp << ")  ";
#endif
		if(minlen > tmp)
		{
		    minlen = tmp;
		    min_j = idxs[j];
		    min_k = idxs[k];
		}
	    }
#ifdef PT_VERBOSE
cout << endl;
#endif
	}
    }

    double PTAlg::distance(size_t i, size_t j) const
    {
#ifdef DEBUG
	assert(i < j);
#endif
	return get_iter(i, j)->second;
    }

    double PTAlg::length(size_t idx, double tmplen) const
    {
	double result = tmplen + 2.0*totdist_[idx];
	return (result > 0.0) ? result : 0.1;
    }

    void PTAlg::update(size_t p, size_t q, const vector<size_t>& vs)
    {
#ifdef DEBUG
	assert(p < q);
#endif
	iterator titer = get_iter(p, q);
	const double tdist = titer->second;
	dm_[p]->erase(titer);

	size_t i = 0;
	for(; vs[i] < p; ++i)
	{
	    const size_t idx = vs[i];
	    iterator m = get_iter(idx, p);
	    iterator n = get_iter(idx, q);
	    double t = (m->second + n->second - tdist) * 0.5;
	    dm_[idx]->erase(n);
	    dm_[idx]->erase(m);
	    dm_[idx]->push_back(elm_type(insidx_, t));
	    totdist_[idx] -= (t+tdist);
	    totdist_[insidx_] += t;
	}
	++i;
	for(; vs[i] < q; ++i)
	{
	    const size_t idx = vs[i];
	    iterator m = get_iter(p, idx);
	    iterator n = get_iter(idx, q);
	    double t = (m->second + n->second - tdist) * 0.5;
	    dm_[p]->erase(m);
	    dm_[idx]->erase(n);
	    dm_[idx]->push_back(elm_type(insidx_, t));
	    totdist_[idx] -= (t+tdist);
	    totdist_[insidx_] += t;
	}
	++i;
	const size_t vssize = vs.size();
	for(; i < vssize; ++i)
	{
	    const size_t idx = vs[i];
	    iterator m = get_iter(p, idx);
	    iterator n = get_iter(q, idx);
	    double t = (m->second + n->second - tdist) * 0.5;
	    dm_[p]->erase(m);
	    dm_[q]->erase(n);
	    dm_[idx]->push_back(elm_type(insidx_, t));
	    totdist_[idx] -= (t+tdist);
	    totdist_[insidx_] += t;
	}
	++insidx_;
#ifdef PT_VERBOSE
	totdist_[p] = 0.0;
	totdist_[q] = 0.0;
#endif
    }

    void PTAlg::show() const
    {
	cout << endl;
	const size_t size = dm_.size();
	for(size_t i = 0; i < size; ++i)
	{
	    cout << "dm_[" << setw(2) << i
		<< "] (" << setw(2) << totdist_[i] << "): ";
	    const size_t dsize = dm_.at(i)->size();
	    for(size_t j = 0; j < dsize; ++j)
	    {
		elm_type t = dm_.at(i)->at(j);
		cout << "(" << setw(2) << t.first
		    << ", " << setw(2) << t.second
		    << ") ";
	    }
	    cout << "\n";
	}
	cout << "\n";
    }
}

