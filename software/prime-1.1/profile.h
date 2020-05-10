// profile.h
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

#ifndef __PROFILE_H__
#define __PROFILE_H__

#include "prrn.h"
#include "residue.h"
#include "alignment.h"
#include "matrix.hpp"

class Profile
{
    public:
	Profile();
	Profile(size_t, size_t);
	Profile(const Profile&);
	~Profile(){}

	Profile& operator=(const Profile&);
	Profile& operator+=(const Profile&);

	size_t getLength() const;

	size_t getNumResidues() const;
	size_t getNumSequences() const;

	prrn::SCORE getProfile(size_t row, size_t col) const
	{
	    return profile_->get(row, col);
	}

	double getFrequency(size_t row, size_t col) const
	{
	    return frequency_->get(row, col);
	}

#ifdef DEBUG
	void showProfile(const prrn::RESIDUE*);
	void showFrequency(const prrn::RESIDUE*);
#endif

	static void setResidueUseRange(prrn::RESCODE, prrn::RESCODE);
	static bool isEqualFrequency(Profile&, Profile&);
	static bool isEqualProfile(Profile&, Profile&);

    protected:
	Matrix<prrn::SCORE>* profile_;
	Matrix<double>* frequency_;
	size_t nsqn_;
	size_t length_;

	static prrn::RESCODE begin_, end_;
};

#endif
