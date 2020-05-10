// substitution_matrix.cpp
//
// Last Modified 7, Mar 2007
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

#include "substitution_matrix.h"
#include "matrix.hpp"
#include "util.h"
#ifdef DEBUG
#include "params.h"
#endif
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <cassert>

using namespace std;
using namespace prrn;

SubstitutionMatrix* SubstitutionMatrix::pinstance_ = 0;

SubstitutionMatrix& SubstitutionMatrix::getInstance(string mat_name, size_t nrestype, string gap_name)
{
    if(!pinstance_)
	pinstance_ = new SubstitutionMatrix(mat_name, nrestype, gap_name);
    return *pinstance_;
}

SCORE SubstitutionMatrix::getOpenPenalty(size_t i) const
{
#ifdef DEBUG
    return go_.at(i);
#else
    return go_[i];
#endif
}

SCORE SubstitutionMatrix::getExtensionPenalty(size_t i) const
{
#ifdef DEBUG
    return ge_.at(i);
#else
    return ge_[i];
#endif
}

SCORE SubstitutionMatrix::getDifferencePenalty(size_t len) const
{
    if(len < threshold_[0])
	return ge_[0];
#ifdef DEBUG
    if(len < threshold_[1])
#else
    else
#endif
	return ge_[1];
#ifdef DEBUG
    static const int tsize = threshold_.size();
    for(int i = 2; i < tsize; ++i)
	if(len < threshold_[i])
	    return ge_[i];
    return limitScoreZero;
#endif
}

SCORE SubstitutionMatrix::getMinPenalty(size_t len) const
{
    if(len == 0)
	return limitScoreZero;

    SCORE t = go_[0] + len*ge_[0];

    static const int gosize = go_.size();
    for(int i = 1; i < gosize; ++i)
	if(t > go_[i] + len*ge_[i])
	    t = go_[i] + len*ge_[i];
    return t;
}

size_t SubstitutionMatrix::getMatSize() const
{
    return matrix_.getRowSize();
}

size_t SubstitutionMatrix::getGapSize() const
{
    return go_.size();
}

size_t SubstitutionMatrix::getNumGapThreshold() const
{
    return threshold_.size();
}

const vector<SCORE>& SubstitutionMatrix::getGapOpen() const
{
    return go_;
}

const vector<SCORE>& SubstitutionMatrix::getGapExtension() const
{
    return ge_;
}

size_t SubstitutionMatrix::getGapThreshold(size_t i) const
{
#ifdef DEBUG
    return threshold_.at(i);
#else
    return threshold_[i];
#endif
}

SubstitutionMatrix::SubstitutionMatrix(const string& mat_name, size_t nrestype, const string& gap_name)
    : nrt_(nrestype)
    , matrix_(nrt_, nrt_)
{
    setPenalty(gap_name);
    setMatrix(mat_name);
}

namespace
{
    void calcNonResidueScore(Matrix<SCORE>& matrix, RESCODE (*getRescode)(RESIDUE), size_t nrt)
    {
	// modification for 'B'
	RESCODE acid = getRescode('D');
	RESCODE nonacid = getRescode('N');
	RESCODE target = getRescode('B');
	for(RESCODE i = PenaltyGap+1; i < nrt; ++i)
	    matrix(target, i) = matrix(i, target)
		= roundDbl((matrix(i, acid) + matrix(i, nonacid)) * 0.5);

	// modification for 'Z'
	acid = getRescode('E');
	nonacid = getRescode('Q');
	target = getRescode('Z');
	for(RESCODE i = PenaltyGap+1; i < nrt; ++i)
	    matrix(target, i) = matrix(i, target)
		= roundDbl((matrix(i, acid) + matrix(i, nonacid)) * 0.5);
    }
}

void SubstitutionMatrix::setMatrix(const string& mat_name)
{
    string tmp;
    istringstream tline;

    ifstream mat_f(mat_name.c_str());
    assert(mat_f && "Substitution matrix can not open.");
    while(1)
    {
	getline(mat_f, tmp);
	if(tmp.size() != 0 && '#' != tmp[0])
	{
	    tline.str(tmp);
	    break;
	}
    }

    const bool isAAmode = (nrt_ == nAACode);
    RESCODE (*getRescode)(RESIDUE) = (isAAmode) ? &getAANum : &getNANum;

    vector<RESCODE> trc;
    trc.reserve(nrt_);
    RESIDUE t_res = '.';
    RESCODE nres = 0;
    while(tline >> t_res)
    {
	trc.push_back(getRescode(t_res));
	++nres;
    }
    assert(isAAmode || (trc.size() <= nNACode
		&& "Sustitution matrix for nucleic acid must be used."));

    for(RESCODE i = 0; i < nres; ++i)
    {
	getline(mat_f, tmp);
	tline.clear();
	tline.str(tmp);

	// first column must be a residue alphabet
	assert((tline >> t_res) && (getRescode(t_res) == trc[i])
		&& "Substitution matrix must be symmetric.");

	for(RESCODE j = 0; j < nres; ++j)
	    assert((tline >> matrix_(trc[i], trc[j]))
		    && "Substitution matrix must be symmetric.");
    }

    // need more refinement
    if(isAAmode && nres < nrt_-1)
	calcNonResidueScore(matrix_, getRescode, nrt_);

    // modification for gap extension penalty
    for(RESCODE i = PenaltyGap+1; i < nrt_; ++i)
    {
	matrix_(PenaltyGap, i) = - ge_[0];
	matrix_(i, PenaltyGap) = - ge_[0];
    }
    matrix_(PenaltyGap, PenaltyGap) = 0.0;
    assert(isSymmetric(matrix_)
	    && "Substitution matrix must be symmetric.");

#ifdef DEBUG
    const RESIDUE* rescode = getResCode(res_type);
    cout << "\nsubstitution matrix:\n";
    cout << "  ";
    for(size_t i = 0; i < nrt_; ++i)
	cout << setw(4) << rescode[i] << " ";
    cout << "\n";
    for(size_t i = 0; i < nrt_; ++i)
    {
	cout << rescode[i] << " ";
	for(size_t j = 0; j < nrt_; ++j)
	    cout << setw(4) << matrix_(i, j) << " ";
	cout << "\n";
    }
    cout << "\n";
#endif
}

void SubstitutionMatrix::setPenalty(const string& gap_name)
{
    ifstream fgap(gap_name.c_str());
    if(!fgap)
    {
	cerr << "\"" << gap_name << "\" cannot open:\n"
	    << "using default parameter: gap open = 9, gap extension = 1\n";
	go_.push_back(9);
	ge_.push_back(1);
	threshold_.push_back(limitSizeMax);
	return;
    }

    string tmp;
    istringstream gline;
    SCORE ge = limitScoreZero;

    // for affine gap cost
    while(getline(fgap, tmp))
    {
	if(tmp.size() == 0 || '#' == tmp[0])
	    continue;
	else
	{
	    SCORE go = limitScoreZero;

	    gline.str(tmp);
	    gline >> go >> ge;
	    go_.push_back(go);
	    ge_.push_back(ge);
	    gline.clear();
	    break;
	}
    }

    // for piecewise linear gap cost
    size_t clen = 0;
    while(getline(fgap, tmp))
    {
	if(tmp.size() == 0 || '#' == tmp[0])
	    continue;
	else
	{
	    gline.str(tmp);
	    gline >> clen >> ge;
	    threshold_.push_back(clen);
	    ge_.push_back(ge);
	    gline.clear();
	}
    }
    threshold_.push_back(limitSizeMax);

    const size_t gsize = ge_.size();
    for(size_t i = 0; i < gsize-1; ++i)
    {
	assert(ge_[i] > ge_[i+1]); 
	assert(threshold_[i] < threshold_[i+1]);
	go_.push_back((ge_[i]-ge_[i+1])*threshold_[i]+go_[i]);
	assert(go_[i] < go_[i+1]);
    }

#ifdef DEBUG
    cout << "\ngap penalty function:\ng(x) = - ";
    if(gsize > 1)
	cout << "min{";
    else
	cout << "(";
    for(size_t i = 0; i < gsize-1; ++i)
	cout << setw(3) << ge_[i] << "x + " << setw(3) << go_[i] << ", ";
    cout << setw(3) << ge_[gsize-1] << "x + " << setw(3) << go_[gsize-1];
    if(gsize > 1)
	cout << "}\n";
    else
	cout << ")\n";

    if(gsize > 1)
    {
	cout << "\nthreshold for piecewise linear gap cost:\n";
	for(size_t i = 0; i < gsize-1; ++i)
	    cout << setw(2) << threshold_[i] << ", ";
	cout << setw(2) << threshold_[gsize-1];
	cout << "\n";
    }
#endif
}
