// profile_dp_long.cpp
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

#include "profile_dp_long.h"
#include "sequence.h"
#include "alignment.h"
#include "avg_profile.h"
#include "pt_holder.h"
#include "path.h"
#ifdef MSA_MATRIX
#include "params.h"
#endif
#include <iostream>
#include <iomanip>
#include <numeric>

using namespace std;
using namespace prrn;

namespace
{
    prrn::SCORE gfp_;
    size_t begin_, end_;
    const Alignment* pta_;
    const Alignment* ptb_;
    size_t anrow_;
    size_t bnrow_;
    vector<double> aweight_;
    vector<double> bweight_;
    double afreq_;
    double bfreq_;
    typedef char MATCHMODE;
    MATCHMODE mode_;

    void showDGI(const vector<pair<int, int> >& rdgiv)
    {
	if(rdgiv.empty())
	{
	    cout << "empty" << endl;
	    return;
	}
	const size_t size = rdgiv.size();
	for(size_t i = 0; i < size; ++i)
	    cout << "(" << setw(2) << rdgiv[i].first << ", "
		<< setw(2) << rdgiv[i].second << ") ";
	cout << endl;
    }
}

ProfileDPLong::ProfileDPLong()
    : MultipleDPAlgorithm()
{
    gfp_ = go_ + ge_;
    begin_ = Profile::begin_;
    end_ = Profile::end_;
}

void ProfileDPLong::showAlignmentMode() const
{
    cerr << "Profile-to-profile alignment mode is long (terminal gaps are penalized)\n";
    cerr << "Gap penalty function g(x) = min{";

    const size_t L = sm_.getGapSize();
    for(size_t i = 0; i < L; ++i)
    {
	cerr << sm_.getExtensionPenalty(i) << "x + ";
	cerr << sm_.getOpenPenalty(i);
	if(i != L-1)
	    cerr << ", ";
    }
    cerr <<"}\n\n";
}

inline size_t ProfileDPLong::getStartPosition(size_t pos, size_t gstat, size_t& dlen,
	vector<DGI>::const_reverse_iterator& i,
	const vector<DGI>::const_reverse_iterator& ie) const
{
    while(i != ie)
    {
	int tlen = pos - i->first + i->second + dlen;
	if(tlen <= (int)gstat)
	    dlen += (i++)->second;
	else
	{
	    return pos - gstat + dlen
		+ std::max((int)gstat - tlen + i->second, 0); 
	}
    }
    return pos - gstat + dlen;
}

double ProfileDPLong::getGapPenalty(const GapProfile& rgp_res, const DynamicGapStates& rdgs_res,
	const GapProfile& rgp_gap, const DynamicGapStates& rdgs_gap,
	vector<DGI>& vec_dgi, const vector<pair<size_t, double> >& range_prof,
	size_t pos) const
{
    vector<gvar>::const_iterator p = rgp_res.gv_.begin();
    const vector<gvar>::const_iterator pe = rgp_res.gv_.end();
    vector<dg>::const_iterator q = rdgs_res.state_.begin();

    vector<gvar>::const_iterator i = rgp_gap.gv_.begin();
    const vector<gvar>::const_iterator ie = rgp_gap.gv_.end();
    vector<dg>::const_iterator j = rdgs_gap.state_.begin();

    vector<DGI>::const_reverse_iterator u = vec_dgi.rbegin();
    const vector<DGI>::const_reverse_iterator ue = vec_dgi.rend();
    vector<pair<size_t, double> >::const_iterator v = range_prof.begin();

    SCORE gep = 0.0;
    SCORE ngap = 0.0;

    size_t dlen = 0;

    size_t tgsi = 0;
    size_t tlen = 0;
    size_t begin = 0;
    double tfreq = 0.0;

    while(i != ie)
    {
	// gap profiles must be sorted
	// always dynamic gap state is found
	tgsi = i->first;
	while(tgsi != j->static_gap_length_)
	    ++j;
#ifdef DEBUG
	assert(j != rdgs_gap.state_.end());
#endif
	tlen = tgsi + (j++)->dynamic_gap_length_;
	tfreq = (i++)->second;

	begin = getStartPosition(pos, tlen, dlen, u, ue);
	while(v->first > begin)
	    ++v;
#ifdef DEBUG
	assert(v != range_prof.end());
#endif
	gep += tfreq * v->second;

	while(p != pe)
	{
	    while(p->first != q->static_gap_length_)
		++q;
#ifdef DEBUG
	    assert(q != rdgs_res.state_.end());
#endif
	    if(tlen <= p->first + q->dynamic_gap_length_)
	    {
		ngap += tfreq * p->second;
		break;
	    }
	    ++p, ++q;
	}
    }

    if(const int ee = ue - u - 1 > 0)
    {
	vector<DGI>::iterator vb = vec_dgi.begin();
	vec_dgi.erase(vb, vb + ee);
    }
    return - (gep + ngap*go_);
}

SCORE ProfileDPLong::calc_diag(const AverageProfile& A, size_t posa, const CountSequence& csa,
	const AverageProfile& B, size_t posb, const CountSequence& csb, const cell4lng& diag) const
{
    SCORE result = getGapPenalty(A.gvv_._gvc_res[posa], diag.c.dgs_a,
	    B.gvv_._gvc_gap[posb], diag.c.dgs_b,
	    diag.dgp_a, csa.getRangeProfile(posa), posa);
    result += getGapPenalty(B.gvv_._gvc_res[posb], diag.c.dgs_b,
	    A.gvv_._gvc_gap[posa], diag.c.dgs_a,
	    diag.dgp_b, csb.getRangeProfile(posb), posb);

    switch(mode_)
    {
	case 2:
	case 3:
	    for(size_t i = 0; i < anrow_; ++i)
		result += B.getProfile(pta_->getResidueCode(i, posa), posb) * aweight_[i];
	    break;
	case 0:
	    for(size_t k = begin_; k < end_; ++k)
		result += A.getProfile(k, posa) * B.getFrequency(k, posb);
	    break;
	default:	// case 1:
	    for(size_t j = 0; j < bnrow_; ++j)
		result += A.getProfile(ptb_->getResidueCode(j, posb), posa) * bweight_[j];
    }
    return result;
}

inline SCORE ProfileDPLong::calc_vert(const AverageProfile& A, size_t posa,
	const GapProfile& bcgp, const cell4lng& vert,
	const CountSequence& csa) const
{
    return getGapPenalty(A.gvv_._gvc_res[posa], vert.c.dgs_a,
	    bcgp, vert.c.dgs_b,
	    vert.dgp_a, csa.getRangeProfile(posa), posa);
}

inline SCORE ProfileDPLong::calc_hori(const GapProfile& acgp,
	const AverageProfile& B, size_t posb, const cell4lng& hori,
	const CountSequence& csb) const
{
    return getGapPenalty(B.gvv_._gvc_res[posb], hori.c.dgs_b,
	    acgp, hori.c.dgs_a,
	    hori.dgp_b, csb.getRangeProfile(posb), posb);
}

namespace
{
    typedef char pathdir;
    const pathdir fromDiagonal = 0;
    const pathdir fromVertical = 1;
    const pathdir fromHorizontal = 2;
}

Alignment ProfileDPLong::getAlignment(const Alignment& ta, const Alignment& tb,
	const PTHolder* pth) const
{
    pta_ = &ta;
    const size_t I = pta_->getLength();
    anrow_ = pta_->getNumSequences();
    mode_  = (anrow_ <= end_) ? 2 : 0;

    aweight_.resize(anrow_, 1.0);
    if(pth && pth->isWeightMode())
	for(size_t i = 0; i < anrow_; ++i)
	    aweight_[i] = pth->getPartialWeight(pta_->getName(i));

    AverageProfile A(*pta_, aweight_);
    CountSequence cs_a(*pta_, aweight_);

    ptb_ = &tb;
    const size_t J = ptb_->getLength();
    bnrow_ = ptb_->getNumSequences();
    mode_  += (bnrow_ <= end_) ? 1 : 0;

    bweight_.resize(bnrow_, 1.0);
    if(pth && pth->isWeightMode())
	for(size_t j = 0; j < bnrow_; ++j)
	    bweight_[j] = pth->getPartialWeight(ptb_->getName(j));

    AverageProfile B(*ptb_, bweight_);
    CountSequence cs_b(*ptb_, bweight_);

    switch(mode_)
    {
	case 2:
	case 3:
	    B.makeProfile(*ptb_, bweight_);
	    break;
	case 0:
	    B.countFrequency(*ptb_, bweight_);
	    // do not insert `break;' here
	default:	// case 1:
	    A.makeProfile(*pta_, aweight_);
    }

    // initialize hori
    vector<cell4lng> max(J+1, 0.0);
    vector<cell4lng> vert(J+1, - limitScoreMax);

    max.begin()->c.dgs_a.set(dg(0, 0));
    max.begin()->c.dgs_b.set(dg(0, 0));
    for(size_t j = 1; j <= J; ++j)
    {
	max[j].c.dgs_a.incl(max[j-1].c.dgs_a);
	max[j].c.dgs_b.update(B.gvv_._gvc_gap[j-1], max[j-1].c.dgs_b);
	max[j].dgp_a.push_back(make_pair(0, j));
	vert[j].c.dgs_a = max[j].c.dgs_a;
	vert[j].c.dgs_b = max[j].c.dgs_b;
	vert[j].dgp_a = max[j].dgp_a;
    }

    // for weights
    afreq_ = accumulate(aweight_.begin(), aweight_.end(), 0.0);
    bfreq_ = accumulate(bweight_.begin(), bweight_.end(), 0.0);

    // for scores
    GapProfile tgp;
    tgp.set(gvar(0, afreq_));
    for(size_t j = 1; j <= J; ++j)
    {
	max[j].c.score = max[j-1].c.score;
	max[j].c.score -= tgp.getNumGapOpens(B.gvv_._gvc_res[j-1]) * go_;
	max[j].c.score -= cs_b.getRangeProfile(j-1).begin()->second * afreq_;
	++tgp.gv_.begin()->first;
    }
    tgp.gv_.begin()->first = 0;
    tgp.gv_.begin()->second = bfreq_;

    // for path
    PathList path_list(I+J);
    path_list.setPath(0, 0, 0);

    vector<pathdir> path_direction(J+1, fromHorizontal);
    pathdir diag_path_direction = fromDiagonal;

#ifdef MSA_MATRIX
    const RESIDUE* rescode = getResCode(res_type);
    for(size_t n = 0; n < bnrow_; ++n)
    {
	cout << setw(anrow_+11) << " ";
	for(size_t j = 0; j < J; ++j)
	    cout << setw(9) << rescode[tb.getResidueCode(n, j)] << " ";
	cout << endl;
    }
    cout << setw(anrow_+1) << " ";
    for(size_t j = 0; j <= J; ++j)
    {
	cout << setw(8) << max[j].c.score;
	switch(path_direction[j])
	{
	    case fromDiagonal:
		cout << "\\ ";
		break;
	    case fromVertical:
		cout << "| ";
		break;
	    case fromHorizontal:
		cout << "_ ";
		break;
	}
    }
    cout << endl;
#endif

    cell4lng tmp_diag(0.0);
    cell4lng tmp_hori(0.0);

    GapProfile acgp, bcgp;

    //
    // alignment phase
    //
    for(size_t i = 1; i <= I; ++i)
    {
	A.gvv_.getCurGapState(i-1, acgp);

	// initialize vert
	tmp_diag = *max.begin();
	tmp_diag.c.path = 0;

	max.begin()->c.score -= tgp.getNumGapOpens(A.gvv_._gvc_res[i-1]) * go_;
	max.begin()->c.score -= cs_a.getRangeProfile(i-1).begin()->second * bfreq_;
	++tgp.gv_.begin()->first;

	max.begin()->c.dgs_a.update(A.gvv_._gvc_gap[i-1], max.begin()->c.dgs_a);
	max.begin()->c.dgs_b.incl(max.begin()->c.dgs_b);
	max.begin()->dgp_b.assign(1, make_pair(0, i));

	tmp_hori = *max.begin();
	tmp_hori.c.path = 0;

#ifdef MSA_MATRIX
	for(size_t m = 0; m < anrow_; ++m)
	    cout << rescode[ta.getResidueCode(m, i-1)];
	cout << " ";
	cout << setw(8) << max[0].c.score << "| ";
#endif

	// begin alignment
	for(size_t j = 1; j <= J; ++j)
	{
	    B.gvv_.getCurGapState(j-1, bcgp);

	    SCORE t_scr_diag = calc_diag(A, i-1, cs_a, B, j-1, cs_b, tmp_diag);

	    SCORE t_scr_vert_vert = calc_vert(A, i-1, bcgp, vert[j], cs_a);
	    SCORE t_scr_vert_diag = calc_vert(A, i-1, bcgp, max[j], cs_a);

	    SCORE t_scr_hori_hori = calc_hori(acgp, B, j-1, tmp_hori, cs_b);
	    SCORE t_scr_hori_diag = calc_hori(acgp, B, j-1, max[j-1], cs_b);

	    // update score and dynamic gap variables
	    //
	    SCORE t_max = tmp_diag.c.score + t_scr_diag;
	    pathdir t_path_direction = fromDiagonal;

	    // update vert
	    if(vert[j].c.score + t_scr_vert_vert >= max[j].c.score + t_scr_vert_diag)
		vert[j].updt_vert(t_scr_vert_vert, A.gvv_, i-1, vert[j], j-1);
	    else
		vert[j].updt_vert(t_scr_vert_diag, A.gvv_, i-1, max[j], j-1);
	    if(t_max < vert[j].c.score)
	    {
		t_max = vert[j].c.score;
		t_path_direction = fromVertical;
	    }

	    // update hori
	    if(tmp_hori.c.score + t_scr_hori_hori >= max[j-1].c.score + t_scr_hori_diag)
		tmp_hori.updt_hori(t_scr_hori_hori, B.gvv_, j-1, tmp_hori, i-1);
	    else
		tmp_hori.updt_hori(t_scr_hori_diag, B.gvv_, j-1, max[j-1], i-1);
	    if(t_max < tmp_hori.c.score)
	    {
		t_max = tmp_hori.c.score;
		t_path_direction = fromHorizontal;
	    }

	    // update diag and max
	    if(t_path_direction == fromDiagonal)
	    {
		tmp_diag.updt_diag(t_scr_diag, A.gvv_, i-1, B.gvv_, j-1, max[j]);

		// set optimal path candidate
		if(diag_path_direction != fromDiagonal)
		{
		    path_list.setPath(i-1, j-1, max[j].c.path);
		    max[j].c.path = path_list.getSize() - 1;
		}
	    }
	    else
	    {
		tmp_diag = max[j];
		if(t_path_direction == fromVertical)
		    max[j] = vert[j];
		else //if(t_path_direction == fromHorizontal)
		    max[j] = tmp_hori;
	    }
	    diag_path_direction = path_direction[j];
	    path_direction[j] = t_path_direction;
#ifdef MSA_MATRIX
	    cout << setw(8) << max[j].c.score << " ";
	    switch(path_direction[j])
	    {
		case fromDiagonal:
		    cout << "\\ ";
		    break;
		case fromVertical:
		    cout << "| ";
		    break;
		case fromHorizontal:
		    cout << "_ ";
		    break;
	    }
#endif
	}
#ifdef MSA_MATRIX
	cout << endl;
#endif
	diag_path_direction = fromVertical;
    }
    // set path
    path_list.setPath(I, J, max.back().c.path);
#ifdef MSA_MATRIX
    cout << "\n";
#endif

#ifdef MSA_TEST
    score_ = max.back().c.score;
#endif
    return Alignment(ta, tb, path_list, pth);
}

SCORE ProfileDPLong::getPSAscore(const Alignment& a, size_t i , size_t j) const
{
    const size_t len = a.getLength();
    int ig = 0, jg = 0;
    SCORE result = 0.0;
    for(size_t k = 0; k < len; ++k)
    {
	RESCODE irc = a.getResidueCode(i, k);
	RESCODE jrc = a.getResidueCode(j, k);

	// calculate score
	if(isResidue(irc) && isResidue(jrc))
	    result += sm_.getScore(irc, jrc);
	else if(!isResidue(irc) && !isResidue(jrc))
	    continue;
	else if(!isResidue(irc) && isResidue(jrc))
	    result -= (ig <= jg) ? gfp_ : sm_.getDifferencePenalty(ig);
	else if(isResidue(irc) && !isResidue(jrc))
	    result -= (ig >= jg) ? gfp_ : sm_.getDifferencePenalty(jg);

	// update gap variable
	ig = isResidue(irc) ? 0 : ig+1;
	jg = isResidue(jrc) ? 0 : jg+1;
    }
    return result;
}

SCORE ProfileDPLong::getSPScore(const Alignment& a, const PTHolder* pth) const
{
    const int nSeq = a.getNumSequences();
    SCORE sum = 0.0;
    if(pth && pth->isWeightMode())
    {
	for(int i = 0; i < nSeq-1; ++i)
	    for(int j = i+1; j < nSeq; ++j)
		sum += getPSAscore(a, i, j)
		    * pth->getPairWeight(a.getName(i), a.getName(j));
    }
    else
    {
	for(int i = 0; i < nSeq-1; ++i)
	    for(int j = i+1; j < nSeq; ++j)
		sum += getPSAscore(a, i, j);
    }
    return sum;
}
