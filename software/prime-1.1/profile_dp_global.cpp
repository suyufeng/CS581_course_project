// profile_dp_global.cpp
//
// Last Modified: 5, Jun 2007
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

#include "profile_dp_global.h"
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
	size_t begin_, end_;
	const Alignment* pta_;
	const Alignment* ptb_;
	size_t anrow_;
	size_t bnrow_;
	std::vector<double> aweight_;
	std::vector<double> bweight_;
	double afreq_;
	double bfreq_;
	double age_;
	double bge_;
	typedef char MATCHMODE;
	MATCHMODE mode_;
}

ProfileDPGlobal::ProfileDPGlobal()
    : MultipleDPAlgorithm()
{
    begin_ = Profile::begin_;
    end_ = Profile::end_;
}

void ProfileDPGlobal::showAlignmentMode() const
{
    cerr << "Profile-to-profile alignment mode is global (terminal gaps are penalized)\n";
    cerr << "Gap penalty function g(x) = " << ge_ << "x + " << go_ << "\n\n";
}

double ProfileDPGlobal::getNumGaps(const GapProfile& rgp_res, const DynamicGapStates& rdgs_res,
	const GapProfile& rgp_gap, const DynamicGapStates& rdgs_gap) const
{
    double fgap = 0.0;

    size_t tgsi = 0, tgsp = 0;
    size_t tlen = 0;
    double tfreq = 0.0;

    vector<gvar>::const_iterator p = rgp_gap.gv_.begin();
    const vector<gvar>::const_iterator pe = rgp_gap.gv_.end();
    vector<dg>::const_iterator q = rdgs_gap.state_.begin();

    vector<dg>::const_iterator j = rdgs_res.state_.begin();
    vector<gvar>::const_iterator i = rgp_res.gv_.begin();
    const vector<gvar>::const_iterator ie = rgp_res.gv_.end();
    while(i != ie)
    {
	// gap profiles must be sorted
	// always dynamic gap state is found
	tgsi = i->first;
	while(tgsi != j->static_gap_length_)
	    ++j;
#ifdef DEBUG
	assert(j != rdgs_res.state_.end());
#endif
	tlen = tgsi + (j++)->dynamic_gap_length_;
	tfreq = (i++)->second;

	while(p != pe)
	{
	    tgsp = p->first;
	    while(tgsp != q->static_gap_length_)
		++q;
#ifdef DEBUG
	    assert(q != rdgs_gap.state_.end());
#endif
	    if(tlen >= tgsp + q->dynamic_gap_length_)
		fgap += tfreq * (p++)->second;
	    else
		break;
	}
    }
    return fgap;
}

SCORE ProfileDPGlobal::calc_diag(const AverageProfile& A, size_t posa,
	const AverageProfile& B, size_t posb, const cell& diag) const
{
    SCORE result = - getNumGaps(A.gvv_._gvc_res[posa], diag.dgs_a,
	    B.gvv_._gvc_gap[posb], diag.dgs_b);
    result -= getNumGaps(B.gvv_._gvc_res[posb], diag.dgs_b,
	    A.gvv_._gvc_gap[posa], diag.dgs_a);
    result *= go_;

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

inline SCORE ProfileDPGlobal::calc_vert(const GapProfile& acgp,
	const GapProfile& bcgp, const cell& vert) const
{
    return - (acgp.getMinGapFreq()*bge_
	    + getNumGaps(acgp, vert.dgs_a, bcgp, vert.dgs_b)*go_);
}

inline SCORE ProfileDPGlobal::calc_hori(const GapProfile& acgp,
	const GapProfile& bcgp, const cell& hori) const
{
    return - (bcgp.getMinGapFreq()*age_
	    + getNumGaps(bcgp, hori.dgs_b, acgp, hori.dgs_a)*go_);
}

namespace
{
    typedef char pathdir;
    const pathdir fromDiagonal = 0;
    const pathdir fromVertical = 1;
    const pathdir fromHorizontal = 2;
}

Alignment ProfileDPGlobal::getAlignment(const Alignment& ta, const Alignment& tb,
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

    ptb_ = &tb;
    const size_t J = ptb_->getLength();
    bnrow_ = ptb_->getNumSequences();
    mode_  += (bnrow_ <= end_) ? 1 : 0;

    bweight_.resize(bnrow_, 1.0);
    if(pth && pth->isWeightMode())
	for(size_t j = 0; j < bnrow_; ++j)
	    bweight_[j] = pth->getPartialWeight(ptb_->getName(j));

    AverageProfile B(*ptb_, bweight_);

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
    vector<cell> max(J+1, 0.0);
    vector<cell> vert(J+1, - limitScoreMax);

    max.begin()->dgs_a.set(dg(0, 0));
    max.begin()->dgs_b.set(dg(0, 0));
    for(size_t j = 1; j <= J; ++j)
    {
	max[j].dgs_a.incl(max[j-1].dgs_a);
	max[j].dgs_b.update(B.gvv_._gvc_gap[j-1], max[j-1].dgs_b);
	vert[j].dgs_a = max[j].dgs_a;
	vert[j].dgs_b = max[j].dgs_b;
    }

    // for weights
    afreq_ = accumulate(aweight_.begin(), aweight_.end(), 0.0);
    bfreq_ = accumulate(bweight_.begin(), bweight_.end(), 0.0);
    age_ = afreq_ * ge_;
    bge_ = bfreq_ * ge_;

    // for scores
    GapProfile tgp;
    tgp.set(gvar(0, afreq_));
    for(size_t j = 1; j <= J; ++j)
    {
	max[j].score = max[j-1].score;
	max[j].score -= tgp.getNumGapOpens(B.gvv_._gvc_res[j-1]) * go_;
	max[j].score -= B.gvv_._gvc_res[j-1].getMinGapFreq() * age_;
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
	cout << setw(8) << max[j].score;
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

    cell tmp_diag(0.0);
    cell tmp_hori(0.0);

    GapProfile acgp, bcgp;

    //
    // alignment phase
    //
    for(size_t i = 1; i <= I; ++i)
    {
	A.gvv_.getCurGapState(i-1, acgp);

	// initialize vert
	tmp_diag = *max.begin();
	tmp_diag.path = 0;

	max.begin()->score -= tgp.getNumGapOpens(A.gvv_._gvc_res[i-1]) * go_;
	max.begin()->score -= A.gvv_._gvc_res[i-1].getMinGapFreq() * bge_;
	++tgp.gv_.begin()->first;

	max.begin()->dgs_a.update(A.gvv_._gvc_gap[i-1], max.begin()->dgs_a);
	max.begin()->dgs_b.incl(max.begin()->dgs_b);

	tmp_hori = *max.begin();
	tmp_hori.path = 0;

#ifdef MSA_MATRIX
	for(size_t m = 0; m < anrow_; ++m)
	    cout << rescode[ta.getResidueCode(m, i-1)];
	cout << " ";
	cout << setw(8) << max[0].score << "| ";
#endif

	// begin alignment
	for(size_t j = 1; j <= J; ++j)
	{
	    B.gvv_.getCurGapState(j-1, bcgp);

	    SCORE t_scr_diag = calc_diag(A, i-1, B, j-1, tmp_diag);

	    SCORE t_scr_vert_vert = calc_vert(A.gvv_._gvc_res[i-1], bcgp, vert[j]);
	    SCORE t_scr_vert_diag = calc_vert(A.gvv_._gvc_res[i-1], bcgp, max[j]);

	    SCORE t_scr_hori_hori = calc_hori(acgp, B.gvv_._gvc_res[j-1], tmp_hori);
	    SCORE t_scr_hori_diag = calc_hori(acgp, B.gvv_._gvc_res[j-1], max[j-1]);

	    // update score and dynamic gap variables
	    //
	    SCORE t_max = tmp_diag.score + t_scr_diag;
	    pathdir t_path_direction = fromDiagonal;

	    // update vert
	    if(vert[j].score + t_scr_vert_vert >= max[j].score + t_scr_vert_diag)
		vert[j].updt_vert(t_scr_vert_vert, A.gvv_, i-1, vert[j]);
	    else
		vert[j].updt_vert(t_scr_vert_diag, A.gvv_, i-1, max[j]);
	    if(t_max < vert[j].score)
	    {
		t_max = vert[j].score;
		t_path_direction = fromVertical;
	    }

	    // update hori
	    if(tmp_hori.score + t_scr_hori_hori >= max[j-1].score + t_scr_hori_diag)
		tmp_hori.updt_hori(t_scr_hori_hori, B.gvv_, j-1, tmp_hori);
	    else
		tmp_hori.updt_hori(t_scr_hori_diag, B.gvv_, j-1, max[j-1]);
	    if(t_max < tmp_hori.score)
	    {
		t_max = tmp_hori.score;
		t_path_direction = fromHorizontal;
	    }

	    // update diag and max
	    if(t_path_direction == fromDiagonal)
	    {
		tmp_diag.updt_diag(t_scr_diag, A.gvv_, i-1, B.gvv_, j-1, max[j]);

		// set optimal path candidate
		if(diag_path_direction != fromDiagonal)
		{
		    path_list.setPath(i-1, j-1, max[j].path);
		    max[j].path = path_list.getSize() - 1;
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
	    cout << setw(8) << max[j].score << " ";
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
    path_list.setPath(I, J, max.back().path);
#ifdef MSA_MATRIX
    cout << "\n";
#endif

#ifdef MSA_TEST
    score_ = max.back().score;
#endif
    return Alignment(ta, tb, path_list, pth);
}

SCORE ProfileDPGlobal::getPSAscore(const Alignment& a, size_t i, size_t j) const
{
    const size_t len = a.getLength();
    int ig = 0, jg = 0;
    SCORE result = 0.0;
    for(size_t k = 0; k < len; ++k)
    {
	RESCODE irc = a.getResidueCode(i, k);
	RESCODE jrc = a.getResidueCode(j, k);

	// calculate score
	result += sm_.getScore(irc, jrc);
	if((ig <= jg && !isResidue(irc) && isResidue(jrc)) ||
		(ig >= jg && isResidue(irc) && !isResidue(jrc)))
	    result -= go_;

	// update gap variable
	ig = isResidue(irc) ? 0 : ig+1;
	jg = isResidue(jrc) ? 0 : jg+1;
    }
    return result;
}

SCORE ProfileDPGlobal::getInterProfileScore(const AverageProfile& a, const AverageProfile& b) const
{
    const size_t len = a.getLength();
#ifdef DEBUG
    assert(len == b.getLength());
#endif

    SCORE sum = - getTotalNumGapOpens(a.gvv_, b.gvv_) * go_;
    for(size_t i = begin_; i < end_; ++i)
	for(size_t j = 0; j < len; ++j)
	    sum += a.getProfile(i, j) * b.getFrequency(i, j);

    return sum;
}

#ifdef MSA_TEST
SCORE ProfileDPGlobal::getSPScore(const Alignment& a, const PTHolder* pth) const
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

#else

SCORE ProfileDPGlobal::getSPScore(const Alignment& ra, const PTHolder* pth) const
{
    const size_t nseq =ra.getNumSequences();
    const size_t len = ra.getLength();
    if(!pth)
    {
	double sum = 0.0;
	for(size_t i = 0; i < nseq-1; ++i)
	    for(size_t j = i+1; j < nseq; ++j)
		sum += getPSAscore(ra, i, j);
	return sum;
    }

    const PhylogeneticTree& pt = pth->getPhylogeneticTree();
    const size_t nsz = pt.getNumNeighbors();
    const size_t nedge = pt.getNumEdges();
    vector<AverageProfile*> vap(nsz, static_cast<AverageProfile*>(0));
    pair<size_t, size_t> p;
    SCORE sum = 0.0;
    for(size_t i = 0; i < nsz; ++i)
    {
	p = pt.getNeighbor(i);
	if(p.first < nseq && p.second < nseq)
	{
	    sum += getPSAscore(ra, p.first, p.second)
		* pt.getEdgeWeight(p.first) * pt.getEdgeWeight(p.second);
	    vap[i] = new AverageProfile(ra, p.first, p.second,
		    pt.getEdgeWeight(p.first), pt.getEdgeWeight(p.second));
	}
	else if(p.first >= nseq && p.second >= nseq)
	{
	    if(i+nseq < nedge+1)
		sum += getInterProfileScore(*vap[p.first-nseq], *vap[p.second-nseq]);
	    else
		sum += pt.getEdgeWeight(nedge-1) * getInterProfileScore(*vap[p.first-nseq], *vap[p.second-nseq]);
	    vap[i] = vap[p.first-nseq];
	    vap[p.first-nseq] = 0;
	    *vap[i] += *vap[p.second-nseq];
	    delete vap[p.second-nseq];
	    vap[p.second-nseq] = 0;
	}
	else
	{
	    const size_t si = (p.first < nseq) ? p.first : p.second;
	    const size_t gi = (p.first < nseq) ? p.second-nseq : p.first-nseq;
	    const double weight = pt.getEdgeWeight(si);
	    vap[i] = new AverageProfile(ra, si, weight);
	    sum -= getTotalNumGapOpens(vap[i]->getGapProfileVector(), vap[gi]->getGapProfileVector()) * go_;
	    for(size_t k = 0; k < len; ++k)
		sum += vap[gi]->getProfile(ra.getResidueCode(si, k), k) * weight;
	    *vap[i] += *vap[gi];
	    delete vap[gi];
	    vap[gi] = 0;
	}
	if(i+nseq < nedge-1)
	    vap[i]->multiplyWeight(pt.getEdgeWeight(i+nseq));
    }
    for_each(vap.begin(), vap.end(), DeleteObject());
    return sum;
}

#endif
