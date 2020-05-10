// profile.cpp
//
// Last Modified: 3, Jun 2007
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

#include "profile.h"
#include "util.h"
#include <iostream>
#include <iomanip>
#include <cassert>

using namespace std;
using namespace prrn;

RESCODE Profile::begin_;
RESCODE Profile::end_;

Profile::Profile()
	: profile_(0)
	, frequency_(0)
	, nsqn_(0)
	, length_(0)
{}

Profile::Profile(size_t nseq, size_t len)
	: profile_(0)
	, frequency_(0)
	, nsqn_(nseq)
	, length_(len)
{}

Profile::Profile(const Profile& rp)
    : profile_(0)
    , frequency_(0)
    , nsqn_(rp.nsqn_)
    , length_(rp.length_)
{
    if(rp.profile_)
	profile_ = new Matrix<SCORE>(*rp.profile_);
    if(rp.frequency_)
	frequency_ = new Matrix<double>(*rp.frequency_);
}

Profile& Profile::operator=(const Profile& rp)
{
    if(this != &rp)
    {
	Matrix<SCORE>* pr = 0;
	Matrix<double>* pf = 0;
	if(rp.profile_)
	    pr = new Matrix<SCORE>(*rp.profile_);
	if(rp.frequency_)
	    pf = new Matrix<double>(*rp.frequency_);
	swap(pr, profile_);
	swap(pf, frequency_);
	delete pr;
	delete pf;
	nsqn_ = rp.nsqn_;
	length_ = rp.length_;
    }
    return *this;
}

Profile& Profile::operator+=(const Profile& rp)
{
#ifdef DEBUG
    assert(length_ == rp.length_);
#endif
    *profile_ += *rp.profile_;
    *frequency_ += *rp.frequency_;
    nsqn_ += rp.nsqn_;
    return *this;
}

size_t Profile::getLength() const
{
    return length_;
}

size_t Profile::getNumResidues() const
{
    return profile_->getRowSize();
}

size_t Profile::getNumSequences() const
{
    return nsqn_;
}

#ifdef DEBUG
void Profile::showProfile(const RESIDUE* rescode)
{
    cout << setw(5) << " ";
    for(int i = 0; i < end_; ++i)
	cout << setw(3) << rescode[i] << " ";
    cout << "\n";
    for(size_t i = 0; i < length_; ++i)
    {
	cout << setw(3) << i << ": ";
	for(int j = 0; j < end_; ++j)
	    cout << setw(3) << profile_->get(j, i) << " ";
	cout << "\n";
    }
    cout << "\n";
}

void Profile::showFrequency(const RESIDUE* rescode)
{
    cout << setw(5) << " ";
    for(int i = 0; i < end_; ++i)
	cout << setw(3) << rescode[i] << " ";
    cout << "\n";
    for(size_t i = 0; i < length_; ++i)
    {
	cout << setw(3) << i << ": ";
	for(int j = 0; j < end_; ++j)
	    cout << setw(3) << frequency_->get(j, i) << " ";
	cout << "\n";
    }
    cout << "\n";
}
#endif

void Profile::setResidueUseRange(RESCODE begin, RESCODE end)
{
#ifdef DEBUG
    assert(begin < end);
#endif
    begin_ = begin;
    end_ = end;
}

bool Profile::isEqualFrequency(Profile& pa, Profile& pb)
{
#ifdef DEBUG
    assert(pa.frequency_ && pb.frequency_);
    assert(pa.length_ == pb.length_);
#endif
    for(size_t i = 0; i < end_; ++i)
	for(size_t j = 0; j < pa.length_; ++j)
	    if(!isApproxEqual(pa.frequency_->get(i, j), pb.frequency_->get(i, j)))
		return false;
    return true;
}

bool Profile::isEqualProfile(Profile& pa, Profile& pb)
{
#ifdef DEBUG
    assert(pa.profile_ && pb.profile_);
    assert(pa.length_ == pb.length_);
#endif
    for(size_t i = 0; i < end_; ++i)
	for(size_t j = 0; j < pa.length_; ++j)
	    if(!isApproxEqual(pa.profile_->get(i, j), pb.profile_->get(i, j)))
		return false;
    return true;
}
