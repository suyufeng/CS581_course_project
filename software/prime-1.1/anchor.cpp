// anchor.cpp
//
// Last Modified: 13, Jul 2007
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

#include "res_conserv.h"
#include "anchor.h"
#include "phylogenetic_tree.h"
#include "avg_profile.h"
#include "profile_dp_global.h"
#include <iostream>
#include <iomanip>
#ifdef DEBUG
#include <cassert>
#endif
#include <iterator>

using namespace std;
using namespace prrn;

// affine gap cost without weight is employed
// because of a speculation that conserved subfamilies are of similar length
void eraseConservedBranch(const Alignment& a, const PhylogeneticTree& pt, vector<int>& branch_order, SCORE pair_score_thr, bool isVerbose)
{
    if(pair_score_thr <= 0.0)
	return;

    if(isVerbose)
	cerr << "\n===== bundling subfamily =====\n";
    else
	cerr << "bundling subfamily " << flush;

    const ProfileDPGlobal pdg;
    const double go = SubstitutionMatrix::getInstance().getOpenPenalty();

    const size_t nseq = a.getNumSequences();
    const size_t len = a.getLength();

    const int nsz = pt.getNumNeighbors();
    vector<AverageProfile*> vap(nsz, static_cast<AverageProfile*>(0));
    pair<size_t, size_t> inp;
    if(isVerbose)
	cerr << setprecision(3);
    for(int i = 0; i < nsz; ++i)
    {
	inp = pt.getNeighbor(i);
	size_t efflen = 0;
	double ps = 0.0;
	bool isConserved = false;

	if(isVerbose)
	{
	    cerr << "MSA(" << setw(3) << i << ") <- ";
	    if(inp.first < nseq)
		cerr << "SEQ(" << setw(3) << inp.first;
	    else
		cerr << "MSA(" << setw(3) << inp.first-nseq;
	    cerr << ") + ";
	    if(inp.second < nseq)
		cerr << "SEQ(" << setw(3) << inp.second;
	    else
		cerr << "MSA(" << setw(3) << inp.second-nseq;
	    cerr << "): ";
	}
	if(inp.first < nseq && inp.second < nseq)
	{
	    for(size_t k = 0; k < len; ++k)
		if(isResidue(a.getResidueCode(inp.first, k))
			|| isResidue(a.getResidueCode(inp.second, k)))
		    ++efflen;
	    ps = pdg.getPSAscore(a, inp.first, inp.second) / (double)efflen;
	    if(ps > pair_score_thr)
	    {
		isConserved = true;
		vap[i] = new AverageProfile(a, inp.first, inp.second);
	    }
	}
	else if(inp.first >= nseq && inp.second >= nseq)
	{
	    if(vap[inp.first-nseq] && vap[inp.second-nseq])
	    {
		const size_t fi = inp.first-nseq;
		const size_t fns = vap[fi]->getNumSequences();
		const size_t si = inp.second-nseq;
		const size_t sns = vap[si]->getNumSequences();
		for(size_t k = 0; k < len; ++k)
		    if(!isApproxEqual(vap[fi]->getFrequency(PenaltyGap, k) + vap[fi]->getFrequency(nonPenaltyGap, k), fns)
			    || !isApproxEqual(vap[si]->getFrequency(PenaltyGap, k) + vap[si]->getFrequency(nonPenaltyGap, k), sns))
			++efflen;
		ps = pdg.getInterProfileScore(*vap[fi], *vap[si]) / (double)(fns * sns * efflen);
		if(ps > pair_score_thr)
		{
		    isConserved = true;
		    vap[i] = vap[fi];
		    vap[fi] = 0;
		    *vap[i] += *vap[si];
		    delete vap[si];
		    vap[si] = 0;
		}
	    }
	}
	// either element of `inp' is less than `nseq'
	// in the current PhylogeneticTree implementaion, `inp.first' is always less than `nseq'
	else
	{
	    const size_t si = (inp.first < nseq) ? inp.first : inp.second;
	    const size_t gi = (inp.second >= nseq) ? inp.second-nseq : inp.first-nseq;
	    if(vap[gi])
	    {
		vap[i] = new AverageProfile(a, si);
		ps = - getTotalNumGapOpens(vap[i]->getGapProfileVector(), vap[gi]->getGapProfileVector()) * go;
		const size_t gns = vap[gi]->getNumSequences();
		for(size_t k = 0; k < len; ++k)
		    if(isResidue(a.getResidueCode(si, k))
			    || !isApproxEqual(vap[gi]->getFrequency(PenaltyGap, k) + vap[gi]->getFrequency(nonPenaltyGap, k), gns))
		    {
			++efflen;
			ps += vap[gi]->getProfile(a.getResidueCode(si, k), k);
		    }
		ps /= (double)(gns * efflen);
		if(ps > pair_score_thr)
		{
		    isConserved = true;
		    *vap[i] += *vap[gi];
		    delete vap[gi];
		    vap[gi] = 0;
		}
		else
		{
		    delete vap[i];
		    vap[i] = 0;
		}
	    }
	}
	if(isConserved)
	    branch_order.erase(remove(branch_order.begin(), 
			remove(branch_order.begin(), branch_order.end(), (int)inp.first),
			(int)inp.second), branch_order.end());

	if(isVerbose)
	{
	    cerr << setw(5) << ps;
	    if(isConserved)
		cerr << " conserved";
	    cerr << endl;
	}
    }
    if(isVerbose)
	cerr << "\n";
    else
	cerr << "... done.\n";
    for_each(vap.begin(), vap.end(), DeleteObject());
}

namespace
{
    const size_t wordsz = 20;
}

void detectNonGapRegion(const Alignment& a, vector<pair<size_t, size_t> >& cr_idx)
{
    const size_t nseq = a.getNumSequences();
    const size_t len = a.getLength();

    size_t start = 0;
    bool isStart = false;
    for(size_t i = 0; i < len; ++i)
    {
	bool isResColumn = true;
	for(size_t j = 0; j < nseq; ++j)
	    if(!isResidue(a.getResidueCode(j, i)))
	    {
		isResColumn = false;
		break;
	    }

	if(!isStart && isResColumn)
	{
	    isStart = true;
	    start = i;
	    continue;
	}
	if(isStart && !isResColumn)
	{
	    isStart = false;
	    if(i - start >= wordsz)
		cr_idx.push_back(make_pair(start, i-1));
	}
    }
    if(isStart && len - start >= wordsz)
	cr_idx.push_back(make_pair(start, len-1));
}

// currently, bad implementation
void detectConservedRegion(const Alignment& a, vector<pair<size_t, size_t> >& cr_idx)
{
    const double mincons = 0.8;

    const size_t nseq = a.getNumSequences();
    const size_t len = a.getLength();

    vector<double> weight(nseq, 1.0/(double)nseq);
    vector<double> conservation(len, 0.0);
    calcREconservation(a, weight, conservation);

    size_t start = 0;
    double accum_cons = 0.0;
    bool isStart = false;
    for(size_t i = 0; i < len; ++i)
    {
	if(isStart)
	{
	    accum_cons += conservation[i];
	    if(conservation[i] <= mincons)
	    {
		isStart = false;
		if(i - start >= wordsz)
		    cr_idx.push_back(make_pair(start, i-1));
	    }
	}
	else
	{
	    if(conservation[i] > mincons)
	    {
		isStart = true;
		start = i;
		accum_cons = conservation[i];
	    }
	}
    }
    if(isStart && len - start >= wordsz)
	cr_idx.push_back(make_pair(start, len-1));
}

void calcAnchorPoint(const Alignment& a, vector<pair<size_t, size_t> >& anchor_idx, anchor_mtd method)
{
    switch(method)
    {
	case NonGap:
	    detectNonGapRegion(a, anchor_idx);
	    return;
	case Conservation:
	    detectConservedRegion(a, anchor_idx);
	    return;
	case NoAnchor:
	default:
	    // nothing to do
	    return;
    }
}

namespace
{
    void calcResidueIndexMatrix(const Alignment& a, size_t nseq, size_t len,
	    vector<vector<size_t> >& idm)
    {
	vector<size_t> tmpidx(nseq, 0);
	for(size_t i = 0; i < len; ++i)
	{
	    idm[i].reserve(nseq);
	    for(size_t j = 0; j < nseq; ++j)
		idm[i].push_back(isResidue(a.getResidueCode(j, i)) ? ++tmpidx[j] : 0);
#ifdef DEBUG
	    assert(idm[i].size() == nseq);
#endif
	}
    }
}

void detectUnchangedRegion(const Alignment& refr, const Alignment& trgt,
	vector<pair<size_t, size_t> >& ur_idx)
{
    const size_t nseq = refr.getNumSequences();
    if(nseq != trgt.getNumSequences())
	return;

    const size_t rlen = refr.getLength();
    vector<vector<size_t> > ridx(rlen);
    calcResidueIndexMatrix(refr, nseq, rlen, ridx);

    const size_t tlen = trgt.getLength();
    vector<vector<size_t> > tidx(tlen);
    calcResidueIndexMatrix(trgt, nseq, tlen, tidx);

    size_t i = 0, j = 0;
    size_t ri = 0, tj = 0;
    size_t nres = 0, nnull = 0;
    size_t bpos = 0;
    bool isStart = false;
    while(i < rlen && j < tlen)
    {
	nres = nnull = 0;
	ri = tj = 0;
	for(size_t k = 0; k < nseq; ++k)
	{
	    if(ridx[i][k] == tidx[j][k])
	    {
		if(ridx[i][k])
		    ++nres;
		else
		    ++nnull;
	    }
	    else if(tidx[j][k] && tidx[j][k] < ridx[i][k])
		++ri;
	    else if(ridx[i][k] && ridx[i][k] < tidx[j][k])
		++tj;
	}

	// both columns are same
	if(nres+nnull == nseq)
	{
	    if(!isStart)
	    {
		isStart = true;
		bpos = j;
	    }
	    ++i, ++j;
	}
	else
	{
	    if(isStart)
	    {
		isStart = false;
		ur_idx.push_back(make_pair(bpos, j-1));
	    }
	    if(nres || (ri && tj))// || (!ri && !tj))
		++i, ++j;
	    else if(ri)
		++j;
	    else //if(tj)
		++i;
	}
    }
    if(isStart)
	ur_idx.push_back(make_pair(bpos, tlen-1));
}

// this implementation is purposely redundant
void eraseUnchangedBranch(const Alignment& refr, const Alignment& trgt, const PhylogeneticTree& pt,
	vector<int>& branch)
{
    const int nseq = refr.getNumSequences();
#ifdef DEBUG
    if(nseq != (int)trgt.getNumSequences())
	return;
#endif

    const int rlen = refr.getLength();
    vector<vector<size_t> > ridx;
    ridx.reserve(rlen);

    const int tlen = trgt.getLength();
    vector<vector<size_t> > tidx;
    tidx.reserve(tlen);

    const int nsz = pt.getNumNeighbors();
    pair<int, int> p;

    // for index specifying vector index; -1 denotes changed partial alignment
    vector<int> unchangedNeighborIndex(nsz, -1);

    // for length operation
    int ri = 0, ti = 0;

    // for residue indeces
    int rlp = 0, tlp = 0;
    size_t trf = 0, trs = 0;
    size_t ttf = 0, tts = 0;
    vector<size_t> tmp_ridx;
    tmp_ridx.reserve(rlen);
    vector<size_t> tmp_tidx;
    tmp_tidx.reserve(tlen);
    for(int i = 0; i < nsz; ++i)
    {
	p = pt.getNeighbor(i);
	ri = ti = 0;
	rlp = rlen - 1;
	tlp = tlen - 1;
	trf = trs = ttf = tts = 0;
	tmp_ridx.clear(), tmp_tidx.clear();
	if(p.first < nseq && p.second < nseq)
	{
	    // ignore terminal null columns
	    while(!(isResidue(refr.getResidueCode(p.first, rlp))
			|| isResidue(refr.getResidueCode(p.second, rlp))))
		--rlp;
	    while(!(isResidue(trgt.getResidueCode(p.first, tlp))
			|| isResidue(trgt.getResidueCode(p.second, tlp))))
		--tlp;

	    while(ri <= rlp && ti <= tlp)
	    {
		// skip the null columns
		if(!(isResidue(refr.getResidueCode(p.first, ri))
			    || isResidue(refr.getResidueCode(p.second, ri))))
		{
		    ++ri;
		    tmp_ridx.push_back(0);
		    continue;
		}
		if(!(isResidue(trgt.getResidueCode(p.first, ti))
			    || isResidue(trgt.getResidueCode(p.second, ti))))
		{
		    ++ti;
		    tmp_tidx.push_back(0);
		    continue;
		}

		size_t pos_rf = isResidue(refr.getResidueCode(p.first, ri)) ? ++trf : 0;
		size_t pos_rs = isResidue(refr.getResidueCode(p.second, ri)) ? ++trs : 0;
		size_t pos_tf = isResidue(trgt.getResidueCode(p.first, ti)) ? ++ttf : 0;
		size_t pos_ts = isResidue(trgt.getResidueCode(p.second, ti)) ? ++tts : 0;
		if(pos_rf == pos_tf && pos_rs == pos_ts)
		{
		    tmp_ridx.push_back(max(pos_rf, pos_rs));
		    tmp_tidx.push_back(max(pos_tf, pos_ts));
		    ++ri, ++ti;
		}
		else
		    goto next_neighbor;
	    }
	}
	else if(p.first < nseq) // (p.second >= nseq) also holds
	{
	    const int ucis = unchangedNeighborIndex[p.second-nseq];
	    if(ucis == -1)
		goto next_neighbor;

	    // ignore terminal null columns
	    while(!isResidue(refr.getResidueCode(p.first, rlp))
			&& ridx[ucis][rlp] == 0)
		--rlp;
	    while(!isResidue(trgt.getResidueCode(p.first, tlp))
			&& tidx[ucis][tlp] == 0)
		--tlp;

	    while(ri <= rlp && ti <= tlp)
	    {
		// skip the null columns
		if(!isResidue(refr.getResidueCode(p.first, ri))
			&& ridx[ucis][ri] == 0)
		{
		    ++ri;
		    tmp_ridx.push_back(0);
		    continue;
		}
		if(!isResidue(trgt.getResidueCode(p.first, ti))
			&& tidx[ucis][ti] == 0)
		{
		    ++ti;
		    tmp_tidx.push_back(0);
		    continue;
		}

		size_t pos_rf = isResidue(refr.getResidueCode(p.first, ri)) ? ++trf : 0;
		size_t pos_tf = isResidue(trgt.getResidueCode(p.first, ti)) ? ++ttf : 0;
		if(pos_rf == pos_tf && ridx[ucis][ri] == tidx[ucis][ti])
		{
		    tmp_ridx.push_back(max(pos_rf, ridx[ucis][ri++]));
		    tmp_tidx.push_back(max(pos_tf, tidx[ucis][ti++]));
		}
		else
		    goto next_neighbor;
	    }
	}
	else if(p.second < nseq) // (p.first >= nseq) also holds
	{
	    const int ucif = unchangedNeighborIndex[p.first-nseq];
	    if(ucif == -1)
		goto next_neighbor;

	    // ignore terminal null columns
	    while(!isResidue(refr.getResidueCode(p.second, rlp))
			&& ridx[ucif][rlp] == 0)
		--rlp;
	    while(!isResidue(trgt.getResidueCode(p.second, tlp))
			&& tidx[ucif][tlp] == 0)
		--tlp;

	    while(ri <= rlp && ti <= tlp)
	    {
		// skip the null columns
		if(!isResidue(refr.getResidueCode(p.second, ri))
			&& ridx[ucif][ri] == 0)
		{
		    ++ri;
		    tmp_ridx.push_back(0);
		    continue;
		}
		if(!isResidue(trgt.getResidueCode(p.second, ti))
			&& tidx[ucif][ti] == 0)
		{
		    ++ti;
		    tmp_tidx.push_back(0);
		    continue;
		}

		size_t pos_rs = isResidue(refr.getResidueCode(p.second, ri)) ? ++trs : 0;
		size_t pos_ts = isResidue(trgt.getResidueCode(p.second, ti)) ? ++tts : 0;
		if(pos_rs == pos_ts && ridx[ucif][ri] == tidx[ucif][ti])
		{
		    tmp_ridx.push_back(max(pos_rs, ridx[ucif][ri++]));
		    tmp_tidx.push_back(max(pos_ts, tidx[ucif][ti++]));
		}
		else
		    goto next_neighbor;
	    }
	}
	else // (p.first >= nseq) and (p.second >= nseq)
	{
	    const int ucif = unchangedNeighborIndex[p.first-nseq];
	    const int ucis = unchangedNeighborIndex[p.second-nseq];
	    if(ucif == -1 || ucis == -1)
		goto next_neighbor;

	    // ignore terminal null columns
	    while(ridx[ucif][rlp] == 0 && ridx[ucis][rlp] == 0)
		--rlp;
	    while(tidx[ucif][tlp] == 0 && tidx[ucis][tlp] == 0)
		--tlp;

	    while(ri <= rlp && ti <= tlp)
	    {
		// skip the null columns
		if(ridx[ucif][ri] == 0 && ridx[ucis][ri] == 0)
		{
		    ++ri;
		    tmp_ridx.push_back(0);
		    continue;
		}
		if(tidx[ucif][ti] == 0 && tidx[ucis][ti] == 0)
		{
		    ++ti;
		    tmp_tidx.push_back(0);
		    continue;
		}

		if(ridx[ucif][ri] == tidx[ucif][ti]
			&& ridx[ucis][ri] == tidx[ucis][ti])
		{
		    tmp_ridx.push_back(max(ridx[ucif][ri], ridx[ucis][ri]));
		    tmp_tidx.push_back(max(tidx[ucif][ti], tidx[ucis][ti]));
		    ++ri, ++ti;
		}
		else
		    goto next_neighbor;
	    }
	}
#ifdef DEBUG
	assert(ri == rlp+1);
	assert(ti == tlp+1);
#endif
	// pad `0' for terminal null columns
	while(ri++ < rlen)
	    tmp_ridx.push_back(0);
	while(ti++ < tlen)
	    tmp_tidx.push_back(0);
#ifdef DEBUG
	assert(ridx.empty() || ridx[0].size() == tmp_ridx.size() && tmp_ridx.size() == (size_t)rlen);
	assert(tidx.empty() || tidx[0].size() == tmp_tidx.size() && tmp_tidx.size() == (size_t)tlen);
#endif
	branch.erase(remove(branch.begin(), 
		    remove(branch.begin(), branch.end(), p.first),
		    p.second), branch.end());
	ridx.push_back(tmp_ridx);
	tidx.push_back(tmp_tidx);
#ifdef DEBUG
	tmp_ridx.erase(remove(tmp_ridx.begin(), tmp_ridx.end(), 0U), tmp_ridx.end());
	tmp_tidx.erase(remove(tmp_tidx.begin(), tmp_tidx.end(), 0U), tmp_tidx.end());
	assert(tmp_ridx == tmp_tidx);
	assert(ridx.size() == tidx.size());
#endif
	unchangedNeighborIndex[i] = (int)ridx.size() - 1;
next_neighbor:
	;
    }
}

namespace
{
    typedef vector<pair<size_t, size_t> >::const_iterator ci_idx_type;
    typedef vector<pair<size_t, size_t> >::const_reverse_iterator cri_idx_type;
}

void calcNonConservedRegion(size_t len, const vector<pair<size_t, size_t> >& cr_idx,
	vector<pair<size_t, size_t> >& ncr_idx)
{
    if(cr_idx.empty())
    {
	ncr_idx.push_back(make_pair(0, len-1));
	return;
    }
    ncr_idx.reserve(cr_idx.size()+1);
    ci_idx_type j = cr_idx.begin();
    const ci_idx_type je = cr_idx.end();

    ///*
    // fixing all conserved region
    const size_t r = 3;
    while(j != je && j->second - j->first + 1 < r)
	++j;
    if(j == je)
    {
	ncr_idx.push_back(make_pair(0, len-1));
	return;
    }
    if(j->first != 0)
	ncr_idx.push_back(make_pair(0, j->first-1));
    size_t start = (j++)->second + 1;
    while(j != je)
    {
	while(j != je && j->second - j->first + 1 < r)
	    ++j;
	if(j == je)
	{
	    ncr_idx.push_back(make_pair(start, len-1));
	    return;
	}
	ncr_idx.push_back(make_pair(start, j->first-1));
	start = (j++)->second + 1;
    }
    if(start < len)
	ncr_idx.push_back(make_pair(start, len-1));
    //*/

    /*
    // truncating both terminals of conserved region (MUSCLE's strategy)
    // ignoring short conserved region
    const size_t r = 5;
    while(j != je && j->second - j->first <= r+r)
	++j;
    if(j == je)
    {
	ncr_idx.push_back(make_pair(0, len-1));
	return;
    }
    if(j->first != 0)
	ncr_idx.push_back(make_pair(0, j->first+r));
    size_t start = (j++)->second - r;
    while(j != je)
    {
	while(j != je && j->second - j->first <= r+r)
	    ++j;
	if(j == je)
	{
	    ncr_idx.push_back(make_pair(start, len-1));
	    return;
	}
	ncr_idx.push_back(make_pair(start, j->first+r));
	start = (j++)->second - r;
    }
    if(start+r+1 < len)
	ncr_idx.push_back(make_pair(start, len-1));
    */

    /*
    // dividing DP matrix using conserved region (MAFFT's strategy)
    const size_t r = 3;
    while(j != je && j->second - j->first + 1 < r)
	++j;
    if(j == je)
    {
	ncr_idx.push_back(make_pair(0, len-1));
	return;
    }

    size_t midpoint = (j->first + j->second) >> 1;
    ncr_idx.push_back(make_pair(0, midpoint));
    size_t start = midpoint+1;
    ++j;
    while(j != je)
    {
	while(j != je && j->second - j->first + 1 < r)
	    ++j;
	if(j == je)
	{
	    ncr_idx.push_back(make_pair(start, len-1));
	    return;
	}
	midpoint = (j->first + j->second) >> 1;
	ncr_idx.push_back(make_pair(start, midpoint));
	start = midpoint+1;
	++j;
    }
    ncr_idx.push_back(make_pair(start, len-1));
    */
}

void decomposeAlignment(const Alignment& a, const vector<pair<size_t, size_t> >& ncr_idx,
	vector<Alignment>& vpmsa)
{
    vpmsa.reserve(ncr_idx.size());
    for(ci_idx_type i = ncr_idx.begin(), ie = ncr_idx.end(); i != ie; ++i)
	vpmsa.push_back(Alignment(a, i->first, i->second));
#ifdef DEBUG
    assert(ncr_idx.size() == vpmsa.size());
#endif
}

void replaceAlignment(const vector<pair<size_t, size_t> >& ncr_idx, const vector<Alignment>& vpmsa,
	Alignment& a)
{
#ifdef DEBUG
    assert(ncr_idx.size() == vpmsa.size());
#endif
    vector<Alignment>::const_reverse_iterator j = vpmsa.rbegin();
    for(cri_idx_type i = ncr_idx.rbegin(), ie = ncr_idx.rend(); i != ie; ++i)
	a.replace(*(j++), i->first, i->second);
}
