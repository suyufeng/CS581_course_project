// avg_profile.cpp
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

#include "avg_profile.h"
#include "prrn.h"
#include "sequence.h"
#include "alignment.h"
#include "substitution_matrix.h"
#include <algorithm>

using namespace std;
using namespace prrn;

AverageProfile::AverageProfile(const Alignment& ra, const vector<double>& weight)
    : Profile(ra.getNumSequences(), ra.getLength())
    , gvv_(length_)
{
    gvv_.countGapState(ra, weight);
}

AverageProfile::AverageProfile(const Alignment& ra, size_t i, double w)
    : Profile(2, ra.getLength())
    , gvv_(length_)
{
    // construct frequency or profile vectors
    frequency_ = new Matrix<double>(end_, length_);
    for(size_t k = 0; k < length_; ++k)
	(*frequency_)(ra.getResidueCode(i, k), k) += w;

    static const SubstitutionMatrix& sm = SubstitutionMatrix::getInstance();
    profile_ = new Matrix<SCORE>(end_, length_);
    for(size_t k = 0; k < length_; ++k)
	for(RESCODE r = begin_; r < end_; ++r)
	    (*profile_)(r, k) = sm.getScore(ra.getResidueCode(i, k), r) * w;

    gvv_.countGapState(ra, i, w);
}

AverageProfile::AverageProfile(const Alignment& ra, size_t i, size_t j,
	double wi, double wj)
    : Profile(2, ra.getLength())
    , gvv_(length_)
{
    // construct frequency or profile vectors
    frequency_ = new Matrix<double>(end_, length_);
    for(size_t k = 0; k < length_; ++k)
    {
	(*frequency_)(ra.getResidueCode(i, k), k) += wi;
	(*frequency_)(ra.getResidueCode(j, k), k) += wj;
    }

    static const SubstitutionMatrix& sm = SubstitutionMatrix::getInstance();
    profile_ = new Matrix<SCORE>(end_, length_);
    for(size_t k = 0; k < length_; ++k)
	for(RESCODE r = begin_; r < end_; ++r)
	    (*profile_)(r, k) = sm.getScore(ra.getResidueCode(i, k), r) * wi
		+ sm.getScore(ra.getResidueCode(j, k), r) * wj;

    gvv_.countGapState(ra, i, wi);
    gvv_.countGapState(ra, j, wj);
    gvv_.arrange();
}

AverageProfile::~AverageProfile()
{
    delete frequency_;
    delete profile_;
}

AverageProfile& AverageProfile::operator=(const AverageProfile& rap)
{
    if(this != &rap)
    {
	Profile::operator=(rap);
	gvv_._gvc_res = rap.gvv_._gvc_res;
	gvv_._gvc_gap = rap.gvv_._gvc_gap;
    }
    return *this;
}

AverageProfile& AverageProfile::operator+=(const AverageProfile& rap)
{
    Profile::operator+=(rap);
    gvv_ += rap.gvv_;
    return *this;
}

void AverageProfile::construct(const Sequence& rs, double weight)
{
    nsqn_ = 1;
    length_ = rs.getLength();

    Matrix<double>*pfreq = new Matrix<double>(end_, length_);
    for(size_t i = 0; i < length_; ++i)
	(*pfreq)(rs.getResidueCode(i), i) = weight;

    static const SubstitutionMatrix& sm = SubstitutionMatrix::getInstance();
    Matrix<SCORE>*pprof = new Matrix<SCORE>(end_, length_);
    for(size_t i = 0; i < length_; ++i)
    {
	const RESCODE rc = rs.getResidueCode(i);
	for(RESCODE j = begin_; j < end_; ++j)
	    (*pprof)(j, i) = sm.getScore(j, rc) * weight;
    }

    swap(frequency_, pfreq);
    swap(profile_, pprof);
    delete pfreq;
    delete pprof;

    gvv_.resize(length_);
    size_t gstat = 0;
    Sequence::const_iterator i = rs.begin();
    Sequence::const_iterator ie = rs.end();
    vector<GapProfile>::iterator g = gvv_._gvc_gap.begin();
    vector<GapProfile>::iterator r = gvv_._gvc_res.begin();
    while(i != ie)
    {
	if(isResidue(*(i++)))
	{
	    (r++)->set(gvar(gstat, weight));
	    ++g;
	    gstat = 0;
	}
	else
	{
	    (g++)->set(gvar(gstat, weight));
	    ++r;
	    ++gstat;
	}
    }
}

void AverageProfile::makeProfile(const Alignment& a, const vector<double>& weight)
{
#ifdef DEBUG
    assert(begin_ < end_);
#endif
    profile_ = new Matrix<SCORE>(end_, length_);

    static const SubstitutionMatrix& sm = SubstitutionMatrix::getInstance();

    const size_t nrow = a.getNumSequences();
    vector<SCORE> tmp_freq(end_, limitScoreZero);
    for(size_t i = 0; i < length_; ++i)
    {
	for(size_t j = 0; j < nrow; ++j)
	    tmp_freq[a.getResidueCode(j, i)] += weight[j];

	for(RESCODE j = begin_; j < end_; ++j)
	    for(RESCODE k = begin_; k < end_; ++k)
		(*profile_)(j, i) += sm.getScore(j, k) * tmp_freq[k];
	fill(tmp_freq.begin(), tmp_freq.end(), limitScoreZero);
    }
}

void AverageProfile::multiplyWeight(double weight)
{
    *profile_ *= weight;
    *frequency_ *= weight;
    gvv_.multiplyWeight(weight);
}

void AverageProfile::mergeProfile(const AverageProfile& a, const AverageProfile& b, double weight)
{
    *this = a;
    *this += b;
    multiplyWeight(weight);
}

bool AverageProfile::isEqualProfile(AverageProfile& a, AverageProfile& b)
{
    return (Profile::isEqualFrequency(a, b)
	    && Profile::isEqualProfile(a, b)
	    && a.gvv_ == b.gvv_);
}

void AverageProfile::countFrequency(const Alignment& ra, const vector<double>& weight)
{
#ifdef DEBUG
    assert(begin_ < end_);
#endif
    frequency_ = new Matrix<double>(end_, length_);
    for(size_t i = 0; i < length_; ++i)
	for(size_t j = 0; j < nsqn_; ++j)
	    (*frequency_)(ra.getResidueCode(j, i), i) += weight[j];
}
