// gap_profile.h
//
// Last Modefied: 4, Jun 2007
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

#ifndef __GAP_PROFILE_H__
#define __GAP_PROFILE_H__

#include "util.h"
#include "util_std.hpp"
#include <algorithm>
#include <vector>

typedef std::pair<size_t, double> gvar;
typedef gvar::first_type gapstate;
typedef gvar::second_type frequency;
typedef std::vector<gvar>::iterator gv_iter;
typedef std::vector<gvar>::reverse_iterator gv_rev_iter;
typedef std::vector<gvar>::const_iterator gv_const_iter;
typedef std::vector<gvar>::const_reverse_iterator gv_const_rev_iter;

struct GapProfile
{
	GapProfile(){}
	~GapProfile(){}

	GapProfile(size_t n)
	{
	    gv_.reserve(n);
	}

	GapProfile& operator-=(const GapProfile&);

	bool operator==(const GapProfile&) const;

	void set(const gvar& rgv)
	{
	    gv_.push_back(rgv);
	}

	void clear()
	{
	    gv_.clear();
	}

	void addProfile(const GapProfile&);

	void addEachGapLength(gapstate);
	bool addGapLength(gapstate, size_t);

	void substEachGapLength(size_t);
	bool substGapLength(gapstate, size_t);

	void addState(gapstate, frequency);

	void multiplyWeight(double);

	void cumulate();
	void decumulate();

	frequency getMinGapFreq() const;
	frequency getGapFreq(gapstate) const;

	frequency getNumGapOpens(const GapProfile&) const;

	void sort();

	void eraseZeroFrequency();

	bool isGapColumn() const;
#ifdef DEBUG
	void show() const;
#endif

	std::vector<gvar> gv_;
};

inline frequency GapProfile::getMinGapFreq() const
{
    return gv_.begin()->second;
}

inline void GapProfile::sort()
{
    std::sort(gv_.begin(), gv_.end(), less1st<size_t, double>());
}

inline void GapProfile::eraseZeroFrequency()
{
    gv_.erase(std::remove_if(gv_.begin(), gv_.end(),
		bind2nd(equal_to_2nd<size_t, double>(), 0.0)), gv_.end());
}

class Alignment;

typedef struct _GapProfileVector
{
    typedef std::vector<GapProfile> vecgp;

    vecgp _gvc_res;	// Gap Variable when a residue of Current position is a RESidue
    vecgp _gvc_gap;	// Gap Variable when a residue of Current position is a GAP

    _GapProfileVector(size_t len = 0)
	: _gvc_res(len)
	, _gvc_gap(len)
    {}
    ~_GapProfileVector(){}

    void resize(size_t len)
    {
	_gvc_res.resize(len);
	_gvc_gap.resize(len);
    }
    void arrange();
    void countGapState(const Alignment&, size_t, double);
    void countGapState(const Alignment&, const std::vector<double>&);

    void getCurGapState(size_t pos, GapProfile& gp) const;

    _GapProfileVector& operator+=(const _GapProfileVector&);
    bool operator==(const _GapProfileVector&) const;
    void multiplyWeight(double);
#ifdef DEBUG
    void show() const;
#endif
} GapProfileVector;

double getTotalNumGapOpens(const GapProfileVector&, const GapProfileVector&);

#endif
