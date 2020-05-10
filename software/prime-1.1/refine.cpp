// refine.cpp
//
// Last Modified: 7, Jun 2007
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
#include "sequence.h"
#include "alignment.h"
#include "pt_holder.h"
#include <iostream>
#include <iomanip>

using namespace std;
using namespace prrn;

void getBranchOrder(int nbrch, vector<int>& branch_order)
{
    for(int i = 0; i < nbrch; ++i)
	branch_order.push_back(i);
}

void refineAlignment(const MultipleDPAlgorithm& mda, size_t inner_itr,
	const PTHolder& pth, vector<int>& branch_order, Alignment& ra)
{
    const int nseq = ra.getNumSequences();
    vector<size_t> div_indx;
    div_indx.reserve(nseq);
    Alignment a(ra);
    Alignment b;
    Alignment result;

    SCORE scoremax = mda.getSPScore(ra, &pth);

    const PhylogeneticTree& pt = pth.getPhylogeneticTree();
    const size_t nbrch = branch_order.size();
    size_t limp = nbrch;
    size_t nitr = 0;

    for(size_t j = 0; j < inner_itr; ++j)
    {
	for(size_t i = 0; i < nbrch; ++i)
	{
	    if(i == limp)
		return;
	    // decide which sequences are separated
	    pth.setDivideBranch(branch_order[i]);
	    if(branch_order[i] >= nseq)
	    {
		div_indx.clear();
		pt.getTreeIndex(branch_order[i]-nseq, div_indx);
		sort(div_indx.begin(), div_indx.end());
	    }
	    else
		div_indx.assign(1, branch_order[i]);

	    a.divide(div_indx, b);

	    result = mda.getAlignment(a, b, &pth);

	    SCORE tscore = mda.getSPScore(result, &pth);
	    cerr << "iteration(" << setw(5) << ++nitr << ")\t"
		<< branch_order[i] << "\tscore = "
		<< setprecision(8) << tscore;

	    if(tscore > scoremax)
	    {
		limp = i;
		scoremax = tscore;
		ra = result;
		cerr << "\taccepted";
	    }
	    else
		cerr << "\trejected";
	    cerr << "\n";
	    a = ra;
	}
	if(limp == nbrch)
	    return;
    }
}
