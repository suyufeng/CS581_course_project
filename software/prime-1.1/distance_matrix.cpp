// distance_matrix.cpp
//
// Last Modified: 7, Mar 2007
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

#include "distance_matrix.h"
#include "util.h"
#include "alignment.h"
#include "pairwise_dp_algorithm.h"
#include <iostream>
#include <iomanip>
#include <cassert>

using namespace std;

void calculateDistanceMatrix(const vector<Sequence>& seqs, const PairwiseDPAlgorithm& pda,
	Matrix<size_t>& dist_mat, bool isVerbose)
{
#ifdef DEBUG
    assert(!seqs.empty());
#endif
    const size_t nseq = seqs.size();

    size_t name_length = 0;
    if(isVerbose)
    {
	cerr << "\n===== calculating distance matrix =====\n";
	for(size_t i = 0; i < nseq; ++i)
	    name_length = max(name_length, seqs[i].getName().size());
    }
    else
	cerr << "\ncalculating distance matrix based on PSA " << flush;

    for(size_t i = 0; i < nseq-1; ++i)
    {
	dist_mat(i, i) = 0;
	for(size_t j = i+1; j < nseq; ++j)
	{
	    dist_mat(i, j) = (size_t)roundDbl(100.0 - 100.0 * getSequenceIdentity(pda.getAlignment(seqs[i], seqs[j]), 0, 1));
	    dist_mat(j, i) = dist_mat(i, j);

	    if(isVerbose)
	    {
		cerr << "distance ";
		cerr << setw(name_length) << seqs[i].getName() << "(" << setw(3) << i << ") and ";
		cerr << setw(name_length) << seqs[j].getName() << "(" << setw(3) << j << ")";
		cerr << ":\t" << setw(2) << dist_mat(i, j) << "\n";
	    }
	}
    }
    if(!isVerbose)
	cerr << "... done.\n";
}

void calculateDistanceMatrix(const Alignment& a, Matrix<size_t>& dist_mat, bool isVerbose)
{
#ifdef DEBUG
    assert(isSymmetric(dist_mat)
	    && dist_mat.getRowSize() == a.getNumSequences());
#endif
    const size_t nseq = a.getNumSequences();

    size_t name_length = 0;

    cerr << "\ncalculating distance matrix based on alignment " << flush;
    if(isVerbose)
    {
	cerr << "\n";
	for(size_t i = 0; i < nseq; ++i)
	    if(name_length < a.getName(i).size())
		name_length = a.getName(i).size();
    }
    for(size_t i = 0; i < nseq-1; ++i)
    {
	dist_mat(i, i) = 0;
	for(size_t j = i+1; j < nseq; ++j)
	{
	    dist_mat(i, j) = (size_t)roundDbl(100.0 - 100.0 * getSequenceIdentity(a, i, j));
	    dist_mat(j, i) = dist_mat(i, j);
	    if(isVerbose)
	    {
		cerr << "distance ";
		cerr << setw(name_length) << a.getName(i) << "(" << setw(3) << i << ") and ";
		cerr << setw(name_length) << a.getName(j) << "(" << setw(3) << j << ")";
		cerr << ":\t" << setw(2) << dist_mat(i, j) << "\n";
	    }
	}
    }
    cerr << "... done.\n";
}
