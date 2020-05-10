// avg_profile.h
//
// Last Modefied: 4, Jun 2007
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

#ifndef __AVG_PROFILE_H__
#define __AVG_PROFILE_H__

#include "profile.h"
#include "gap_profile.h"
#include <vector>

class Sequence;
class Alignment;
class ProfileDPGlobal;
class ProfileDPLong;

class AverageProfile : public Profile
{
    private:
	friend class ProfileDPAlgorithm;
	friend class ProfileDPGlobal;
	friend class ProfileDPLong;

    public:
	AverageProfile()
	    : Profile()
	    , gvv_()
	    {}
	AverageProfile(const AverageProfile& rap)
	    : Profile(rap)
	    , gvv_(rap.gvv_)
	    {}
	AverageProfile(const Alignment&, const std::vector<double>&);
	AverageProfile(const Alignment&, size_t, double = 1.0);
	AverageProfile(const Alignment&, size_t, size_t, double = 1.0, double = 1.0);
	~AverageProfile();

	AverageProfile& operator=(const AverageProfile&);
	AverageProfile& operator+=(const AverageProfile&);

	const GapProfileVector& getGapProfileVector() const
	{
	    return gvv_;
	}

	void countFrequency(const Alignment&, const std::vector<double>&);

	void construct(const Sequence&, double);
	void makeProfile(const Alignment&, const std::vector<double>&);

	void multiplyWeight(double);
	void mergeProfile(const AverageProfile&, const AverageProfile&, double);

	static bool isEqualProfile(AverageProfile&, AverageProfile&);

    private:
	GapProfileVector gvv_;
};

#endif
