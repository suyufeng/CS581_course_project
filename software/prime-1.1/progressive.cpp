// progressive.cpp
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

#include "multiple_dp_algorithm.h"
#include "pairwise_dp_algorithm.h"
#include "sequence.h"
#include "alignment.h"
#include "pt_holder.h"
#include "util_std.hpp"
#include "util.h"
#include <iostream>
#include <iomanip>
#include <cassert>

using namespace std;
using namespace prrn;

void getProgressiveAlignment(const vector<Sequence>& seqs,
	const PairwiseDPAlgorithm& pda, const MultipleDPAlgorithm& mda,
	const PTHolder& pth, Alignment& A, bool isVerbose)
{
#ifdef DEBUG
    assert(!seqs.empty());
#endif

    const PhylogeneticTree& pt = pth.getPhylogeneticTree();

    int nseq = seqs.size();

    int name_length = seqs[0].getName().size();
    for(int i = 1; i < nseq; ++i)
	if(name_length < static_cast<int>(seqs[i].getName().size()))
	    name_length = seqs[i].getName().size();

    if(isVerbose)
	cerr << "\n===== constructing multiple alignment =====\n";
    else
	cerr << "constructing multiple alignment " << flush;

    Alignment* a = 0;
    Alignment* b = 0;

    vector<Alignment*> va;
    va.reserve(2*nseq-3);

    size_t vasize = 0;
    const int nsize = pt.getNumNeighbors();

    const int window = max(name_length, 10);
    for(int i = 0; i < nsize; ++i)
    {
	pair<int, int> p = pt.getNeighbor(i);
	if(pth.isWeightMode())
	    pth.setDivideBranch(p.first);

	if(isVerbose)
	{
	    cerr << " alignment(" << setw(3) << vasize++ << ") = ";
	    if(p.first < nseq)
		cerr << setw(window) << seqs[p.first].getName() << "("
		    << setw(3) << p.first << ") + ";
	    else
		cerr << setw(window) << " alignment("
		    << setw(3) << p.first-nseq << ") + ";
	    if(p.second < nseq)
		cerr << setw(window) << seqs[p.second].getName() << "("
		    << setw(3) << p.second << ") ";
	    else
		cerr << setw(window) << " alignment("
		    << setw(3) << p.second-nseq << ") ";
	}
	cerr << flush;
#ifdef MSA_MATRIX
	cout << "\n";
#endif
#ifndef MSA_TEST
	if(p.first < nseq && p.second < nseq)
	    va.push_back(new Alignment(pda.getAlignment(seqs[p.first], seqs[p.second])));
	else
	{
#endif
	    if(p.first < nseq)
		a = new Alignment(seqs[p.first]);
	    else
		a = va[p.first-nseq];
	    if(p.second < nseq)
		b = new Alignment(seqs[p.second]);
	    else
		b = va[p.second-nseq];

#ifndef MSA_TEST
	    va.push_back(new Alignment(mda.getAlignment(*a, *b, &pth)));
	}
#else
	va.push_back(new Alignment(mda.getAlignment(*a, *b, &pth)));
	SCORE sum = mda.getSPScore(*va.back(), &pth);
	if(p.first >= nseq)
	    sum -= mda.getSPScore(*a, &pth);
	if(p.second >= nseq)
	    sum -= mda.getSPScore(*b, &pth);
	cerr << "score test:\t";
	if(fabs(mda.getDPScore() - sum) < 1e-8)
	    cerr << "ok";
	else
	{
	    cerr << "failed ("
		<< setw(8) << mda.getDPScore() << ", "
		<< setw(8) << sum << ")";
	}
#endif
	if(isVerbose)
	    cerr << "\n";
	if(p.first < nseq && p.second < nseq)
	    continue;
	if(p.first < nseq)
	    delete a;
	if(p.second < nseq)
	    delete b;
	//cerr << "\n";
    }
    A = *va.back();
    if(isVerbose)
	cerr << "\n";
    else
	cerr << "... done.\n";
    for_each(va.begin(), va.end(), DeleteObject());
}
