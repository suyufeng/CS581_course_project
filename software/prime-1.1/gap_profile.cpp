// gap_profile.cpp
//
// Last Modified: 4, Jun 2007
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

#include "gap_profile.h"
#include "residue.h"
#include "alignment.h"
#include <iostream>
#include <iomanip>
#ifdef DEBUG
#include <cassert>
#endif

using namespace std;
using namespace prrn;

GapProfile& GapProfile::operator-=(const GapProfile& rg)
{
    const gv_const_iter ie = rg.gv_.end();
    for(gv_const_iter i = rg.gv_.begin(); i != ie; ++i)
    {
	gv_iter j = find_if(gv_.begin(), gv_.end(),
		bind2nd(equal_to_1st<size_t, double>(), i->first));
#ifdef DEBUG
	assert(j != gv_.end());
#endif
	j->second -= i->second;
    }
    return *this;
}

bool GapProfile::operator==(const GapProfile& rg) const
{
    gv_const_iter i = gv_.begin();
    const gv_const_iter ie = gv_.end();
    gv_const_iter j = rg.gv_.begin();
    const gv_const_iter je = rg.gv_.end();

    while(i != ie && j != je)
    {
	if(i->first != j->first
		|| !isApproxEqual(i->second, j->second))
	    return false;
	++i;
	++j;
    }
    return true;
}

void GapProfile::addProfile(const GapProfile& rg)
{
    static vector<gvar> tgp;
    tgp.clear();

    gv_iter i = gv_.begin();
    const gv_iter ie = gv_.end();
    gv_const_iter j = rg.gv_.begin();
    const gv_const_iter je = rg.gv_.end();

    while(i != ie && j != je)
    {
	if(i->first > j->first)
	{
	    tgp.push_back(*j);
	    ++j;
	}
	else if(i->first < j->first)
	{
	    tgp.push_back(*i);
	    ++i;
	}
	else //if(i->first == j->first)
	{
	    tgp.push_back(gvar(i->first, i->second + j->second));
	    ++i, ++j;
	}
    }
    if(i == ie)
	copy(j, je, back_inserter(tgp));
    else //if(j == je)
	copy(i, ie, back_inserter(tgp));
    gv_.swap(tgp);
}

void GapProfile::addEachGapLength(gapstate len)
{
    const gv_iter ie = gv_.end();
    for(gv_iter i = gv_.begin(); i != ie; ++i)
	i->first += len;
}

bool GapProfile::addGapLength(gapstate ngap, size_t len)
{
    const gv_iter ie = gv_.end();
    gv_iter i = find_if(gv_.begin(), ie,
	    bind2nd(equal_to_1st<size_t, double>(), ngap));
    if(i != ie)
    {
	i->first += len;
	return true;
    }
    else
	return false;
}

void GapProfile::substEachGapLength(size_t len)
{
    const gv_iter ie = gv_.end();
    for(gv_iter i = gv_.begin(); i != ie; ++i)
	i->first -= len;
}

bool GapProfile::substGapLength(gapstate ngap, size_t len)
{
    const gv_iter ie = gv_.end();
    gv_iter i = find_if(gv_.begin(), ie,
	    bind2nd(equal_to_1st<size_t, double>(), ngap));
    if(i != ie)
    {
	i->first -= len;
	return true;
    }
    else
	return false;
}

void GapProfile::addState(gapstate ngap, frequency weight)
{
    const gv_iter ie = gv_.end();
    gv_iter i = find_if(gv_.begin(), ie,
	    bind2nd(equal_to_1st<size_t, double>(), ngap));
    if(i != ie)
	i->second += weight;
    else
	gv_.push_back(gvar(ngap, weight));
}

void GapProfile::multiplyWeight(double weight)
{
    const gv_iter ie = gv_.end();
    for(gv_iter i = gv_.begin(); i != ie; ++i)
	i->second *= weight;
}

void GapProfile::cumulate()
{
    for(int i = gv_.size()-1; i > 0; --i)
	gv_[i-1].second += gv_[i].second;
}

void GapProfile::decumulate()
{
    const int size = gv_.size();
    for(int i = 1; i < size; ++i)
	gv_[i-1].second -= gv_[i].second;
}

frequency GapProfile::getGapFreq(gapstate ngap) const
{
    const gv_const_iter ie = gv_.end();
    gv_const_iter i = find_if(gv_.begin(), ie,
	    bind2nd(equal_to_1st<size_t, double>(), ngap));
    return (i != ie) ? i->second : 0.0;
}

frequency GapProfile::getNumGapOpens(const GapProfile &rg) const
{
    double fgap = 0.0;
    const gv_const_iter ie = gv_.end();
    gv_const_iter j = rg.gv_.begin();
    const gv_const_iter je = rg.gv_.end();
    for(gv_const_iter i = gv_.begin(); i != ie; ++i)
    {
	for(; j != je; ++j)
	{
	    if(j->first >= i->first)
	    {
		fgap += j->second * i->second;
		break;
	    }
	}
    }
    return fgap;
}

bool GapProfile::isGapColumn() const
{
    const gv_const_iter ie = gv_.end();
    for(gv_const_iter i = gv_.begin(); i != ie; ++i)
	if(!isApproxEqual(i->second, 0.0))
	    return false;
    return true;
}

#ifdef DEBUG
void GapProfile::show() const
{
    if(gv_.empty())
    {
	cout << "emtpy\n";
	return;
    }

    const gv_const_iter ie = gv_.end();
    for(gv_const_iter i = gv_.begin(); i != ie; ++i)
	cout << "(" << setw(2) << i->first << ", "
	    << setw(2) << i->second << ") "; 
    cout << "\n";
}
#endif

void GapProfileVector::arrange()
{
    for_each(_gvc_res.begin(), _gvc_res.end(), mem_fun_ref(&GapProfile::sort));
    for_each(_gvc_res.begin(), _gvc_res.end(), mem_fun_ref(&GapProfile::cumulate));
    for_each(_gvc_gap.begin(), _gvc_gap.end(), mem_fun_ref(&GapProfile::sort));
}

void GapProfileVector::countGapState(const Alignment& ra, size_t i, double weight)
{
    const size_t len = ra.getLength();
    size_t tg = 0;
    vector<GapProfile>::iterator ri = _gvc_res.begin();
    vector<GapProfile>::iterator gi = _gvc_gap.begin();
    for(size_t k = 0; k < len; ++k, ++ri, ++gi)
    {
	if(isResidue(ra.getResidueCode(i, k)))
	{
	    ri->addState(tg, weight);
	    tg = 0;
	}
	else
	    gi->addState(tg++, weight);
    }
#ifdef DEBUG
    assert(ri == _gvc_res.end());
    assert(gi == _gvc_gap.end());
#endif
}

void GapProfileVector::countGapState(const Alignment& ra, const vector<double>& weight)
{
    const size_t nseq = ra.getNumSequences();
    const size_t len = ra.getLength();
    size_t tgs = 0;
    for(size_t i = 0; i < nseq; ++i)
    {
	tgs = 0;
	for(size_t j = 0; j < len; ++j)
	{
	    if(isResidue(ra.getResidueCode(i, j)))
	    {
		_gvc_res[j].addState(tgs, weight[i]);
		tgs = 0;
	    }
	    else
		_gvc_gap[j].addState(tgs++, weight[i]);
	}
    }
    arrange();
}

void GapProfileVector::getCurGapState(size_t pos, GapProfile& gp) const
{
    gp.gv_.resize(_gvc_gap[pos].gv_.size()+1);
    gv_iter j = gp.gv_.begin();
    *(j++) = gvar(0, _gvc_res[pos].getMinGapFreq());
    gv_const_iter ie = _gvc_gap[pos].gv_.end();
    for(gv_const_iter i = _gvc_gap[pos].gv_.begin(); i != ie; ++i)
	*(j++) = gvar(i->first+1, i->second);
}

GapProfileVector& GapProfileVector::operator+=(const GapProfileVector& rgv)
{
    GapProfile t;
    vecgp::iterator ie = _gvc_res.end();
    vecgp::const_iterator j = rgv._gvc_res.begin();
    for(vecgp::iterator i = _gvc_res.begin(); i != ie; ++i)
    {
	t = *(j++);
	t.decumulate();
	i->decumulate();
	i->addProfile(t);
	i->cumulate();
    }

    ie = _gvc_gap.end();
    j = rgv._gvc_gap.begin();
    for(vecgp::iterator i = _gvc_gap.begin(); i != ie; ++i)
	i->addProfile(*(j++));

    return *this;
}

bool GapProfileVector::operator==(const GapProfileVector& rgv) const
{
    return (_gvc_res == rgv._gvc_res
	    && _gvc_gap == rgv._gvc_gap);
}

void GapProfileVector::multiplyWeight(double weight)
{
    vecgp::iterator ie = _gvc_res.end();
    for(vecgp::iterator i = _gvc_res.begin(); i != ie; ++i)
	i->multiplyWeight(weight);

    ie = _gvc_gap.end();
    for(vecgp::iterator i = _gvc_gap.begin(); i != ie; ++i)
	i->multiplyWeight(weight);
}

double getTotalNumGapOpens(const GapProfileVector& a, const GapProfileVector& b)
{
    double ng = 0.0;
    vector<GapProfile>::const_iterator pe = a._gvc_gap.end();
    for(vector<GapProfile>::const_iterator
	    p = a._gvc_gap.begin(), q = a._gvc_res.begin(),
	    r = b._gvc_gap.begin(), s = b._gvc_res.begin();
	    p != pe;
	    ++p, ++q, ++r, ++s)
    {
	ng += p->getNumGapOpens(*s);
	ng += r->getNumGapOpens(*q);
    }
    return ng;
}

#ifdef DEBUG
void GapProfileVector::show() const
{
    size_t k = 0;
    const vector<GapProfile>::const_iterator ie = _gvc_res.end();
    cerr << "gap state variables\n";
    for(vector<GapProfile>::const_iterator
	    i = _gvc_res.begin(), j = _gvc_gap.begin();
	    i != ie; ++i, ++j)
    {
	cerr << k++ << "\n";

	cerr << "\tS+: ";
	i->show();

	cerr << "\tT : ";
	j->show();
    }
    cerr << "\n";
}
#endif
