// dynamic_gap_state.h
//
// Last Modified: 10, May 2007
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

#ifndef __DYNAMIC_GAP_STATE_H__
#define __DYNAMIC_GAP_STATE_H__

#include "prrn.h"
#include "gap_profile.h"
#include <vector>
#include <functional>

typedef struct _dg	// for one Dynamic Gap
{
    _dg(size_t s = 0, size_t g = 0)
	: static_gap_length_(s)
	, dynamic_gap_length_(g)
	{}

    void incl()
    {
	++dynamic_gap_length_;
    }

    bool isEqualStaticGapLength(size_t len) const
    {
	return (static_gap_length_ == len);
    }

    size_t static_gap_length_;
    size_t dynamic_gap_length_;
} dg;

class GapProfile;
class ProfileDPLong;
class ProfileDPGlobal;

class DynamicGapStates
{
    private:
	friend class ProfileDPLong;
	friend class ProfileDPGlobal;

    public:
	DynamicGapStates(){}
	~DynamicGapStates(){}

	void set(const dg& rdg);

	size_t getDynamicGapLength(size_t) const;

	void incl(const DynamicGapStates&);
	void update(const GapProfile&, const DynamicGapStates&);

	void swap(DynamicGapStates& rdgs)
	{
	    state_.swap(rdgs.state_);
	}

	void clear()
	{
	    state_.clear();
	}

#ifdef DEBUG
	void show() const;
#endif
    private:
	typedef std::vector<dg>::iterator iterator;
	typedef std::vector<dg>::const_iterator const_iterator;
	std::vector<dg> state_;
};

inline size_t DynamicGapStates::getDynamicGapLength(size_t stat_len) const
{
    const const_iterator cie = state_.end();
    const_iterator ci = find_if(state_.begin(), cie,
	    bind2nd(std::mem_fun_ref(&dg::isEqualStaticGapLength), stat_len));
    return  (ci != cie) ? ci->dynamic_gap_length_ : 0;
}

inline void DynamicGapStates::incl(const DynamicGapStates& rdgs)
{
    if(this != &rdgs)
	state_ = rdgs.state_;
    for_each(state_.begin(), state_.end(), std::mem_fun_ref(&dg::incl));
}

typedef struct _cell
{
    prrn::SCORE score;
    DynamicGapStates dgs_a;
    DynamicGapStates dgs_b;
    size_t path;

    _cell(prrn::SCORE = prrn::limitScoreZero);
    ~_cell(){}

    void swap(_cell&);

    void updt_diag(const prrn::SCORE&,
	    const GapProfileVector&, size_t,
	    const GapProfileVector&, size_t, _cell&);
    void updt_vert(const prrn::SCORE&,
	    const GapProfileVector&, size_t, _cell&);
    void updt_hori(const prrn::SCORE&,
	    const GapProfileVector&, size_t, _cell&);
} cell;

inline void cell::swap(_cell& rc)
{
    using std::swap;
    swap(score, rc.score);
    dgs_a.swap(rc.dgs_a);
    dgs_b.swap(rc.dgs_b);
    swap(path, rc.path);
}

inline void cell::updt_diag(const prrn::SCORE& scr,
	const GapProfileVector& gva, size_t posa,
	const GapProfileVector& gvb, size_t posb, _cell& source)
{
    dgs_a.update(gva._gvc_gap[posa], dgs_a);
    dgs_b.update(gvb._gvc_gap[posb], dgs_b);

    swap(source);
    source.score += scr;
}

inline void cell::updt_vert(const prrn::SCORE& scr,
	const GapProfileVector& gv, size_t pos, _cell& source)
{
    if(this != &source)
    {
	score = source.score + scr;
	path = source.path;
    }
    else
	score += scr;

    dgs_a.update(gv._gvc_gap[pos], source.dgs_a);
    dgs_b.incl(source.dgs_b);
}

inline void cell::updt_hori(const prrn::SCORE& scr,
	const GapProfileVector& gv, size_t pos, _cell& source)
{
    if(this != &source)
    {
	score = source.score + scr;
	path = source.path;
    }
    else
	score += scr;

    dgs_a.incl(source.dgs_a);
    dgs_b.update(gv._gvc_gap[pos], source.dgs_b);
}

typedef struct _cell4lng
{
    cell c;
    mutable std::vector<std::pair<int, int> > dgp_a;
    mutable std::vector<std::pair<int, int> > dgp_b;

    _cell4lng(const prrn::SCORE& = prrn::limitScoreZero);
    ~_cell4lng(){}

    _cell4lng(const _cell4lng&);
    _cell4lng& operator=(const _cell4lng&);

    static void updt_dgp(int pos, std::vector<std::pair<int, int> >&);

    void updt_diag(const prrn::SCORE& scr,
	    const GapProfileVector&, size_t,
	    const GapProfileVector&, size_t, _cell4lng&);
    void updt_vert(const prrn::SCORE&, const GapProfileVector&, size_t, _cell4lng&, size_t);
    void updt_hori(const prrn::SCORE&, const GapProfileVector&, size_t, _cell4lng&, size_t);
#ifdef DEBUG
    void show() const;
#endif
} cell4lng;

inline _cell4lng& cell4lng::operator=(const _cell4lng& rc)
{
    if(this != &rc)
    {
	c = rc.c;
	dgp_a = rc.dgp_a;
	dgp_b = rc.dgp_b;
    }
    return *this;
}

inline void cell4lng::updt_dgp(int pos, std::vector<std::pair<int, int> >& rdgp)
{
    if(!rdgp.empty() && rdgp.back().first-1 == pos)
	++rdgp.back().second;
    else
	rdgp.push_back(std::make_pair(pos+1, 1));
}

inline void cell4lng::updt_diag(const prrn::SCORE& scr,
	const GapProfileVector& gva, size_t posa,
	const GapProfileVector& gvb, size_t posb, _cell4lng& source)
{
    c.updt_diag(scr, gva, posa, gvb, posb, source.c);
    dgp_a.swap(source.dgp_a);
    dgp_b.swap(source.dgp_b);
}

inline void cell4lng::updt_vert(const prrn::SCORE& scr,
	const GapProfileVector& gv, size_t pos, _cell4lng& source, size_t vpos)
{
    c.updt_vert(scr, gv, pos, source.c);
    if(this != &source)
    {
	dgp_a = source.dgp_a;
	dgp_b = source.dgp_b;
    }
    updt_dgp(vpos, dgp_b);
}

inline void cell4lng::updt_hori(const prrn::SCORE& scr,
	const GapProfileVector& gv, size_t pos, _cell4lng& source, size_t hpos)
{
    c.updt_hori(scr, gv, pos, source.c);
    if(this != &source)
    {
	dgp_a = source.dgp_a;
	dgp_b = source.dgp_b;
    }
    updt_dgp(hpos, dgp_a);
}

#endif
