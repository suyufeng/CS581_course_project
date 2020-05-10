// oligo_count.cpp
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

#include "oligo_count.h"
#include "prrn.h"
#include "matrix.hpp"
#include "sequence.h"
#include "util.h"
#include "util_std.hpp"
#include <iostream>
#include <iomanip>
#include <numeric>
#include <cassert>

using namespace std;
using namespace prrn;

size_t Oligomer::kmer_ = 6;
const char* Oligomer::pRC_ = se_b;
size_t Oligomer::nrc_ = nrc_se_b;

Oligomer::Oligomer(const Sequence& s)
    : nsgmt_(s.getLength()-kmer_+1)
    , oligo_()
{
    oligo_.reserve(nsgmt_);
    size_t idx = pRC_[s.getResidueCode(0)];
    for(size_t i = 1; i < kmer_; ++i)
    {
	idx *= nrc_;
	idx += pRC_[s.getResidueCode(i)];
    }
    oligo_.push_back(idx);

    const size_t dim = kmer_-1;
    for(size_t i = 1; i < nsgmt_; ++i)
    {
	size_t j = i-1;
	idx -= pRC_[s.getResidueCode(j)] * iexp(nrc_, dim);
	idx *= nrc_;
	idx += pRC_[s.getResidueCode(j+kmer_)];
	oligo_.push_back(idx);
    }
    sort(oligo_.begin(), oligo_.end());
}

void Oligomer::setResidueClass(const char* prc, size_t nrc)
{
    assert(static_cast<double>(iexp(nrc, kmer_)) <= limitDoubleMax);
    pRC_ = prc;
    nrc_ = nrc;
}

void Oligomer::setOligomerLength(size_t k)
{
    assert(static_cast<double>(iexp(nrc_, k)) <= limitDoubleMax);
    kmer_ = k;
}

double Oligomer::getApproxSeqID(const Oligomer& a, const Oligomer& b)
{
    size_t ncommon = 0;

    vector<size_t>::const_iterator i = a.oligo_.begin();
    vector<size_t>::const_iterator j = b.oligo_.begin();

    const vector<size_t>::const_iterator ie = a.oligo_.end();
    const vector<size_t>::const_iterator je = b.oligo_.end();

    while(i != ie && j != je)
    {
	if(*i < *j)
	    ++i;
	else if(*i > *j)
	    ++j;
	else //if(*i == *j)
	{
	    ++ncommon;
	    ++i, ++j;
	}
    }

    return (double)ncommon / (double)min(a.nsgmt_, b.nsgmt_);
}
 
Oligomer* makeOligomer(const Sequence& s)
{
    return new Oligomer(s);
}

void calculateDistanceMatrix(const vector<Sequence>& seqs, 
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
	cerr << "\ncalculating distance matrix based on oligomer counting " << flush;

    vector<Oligomer*> oligo(nseq);
    transform(seqs.begin(), seqs.end(), oligo.begin(), makeOligomer);
    for(size_t i = 0; i < nseq-1; ++i)
    {
	dist_mat(i, i) = 0;
	for(size_t j = i+1; j < nseq; ++j)
	{
	    dist_mat(i, j) = (size_t)roundDbl(100.0 - 100.0 * Oligomer::getApproxSeqID(*oligo[i], *oligo[j]));
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
    for_each(oligo.begin(), oligo.end(), DeleteObject());
}
