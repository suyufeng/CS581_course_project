// substitution_matrix.h
//
// Last Modified 12, Feb 2007
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

#ifndef __SUBSTITUTION_MATRIX_H__
#define __SUBSTITUTION_MATRIX_H__

#include "prrn.h"
#include "residue.h"
#include "matrix.hpp"
#include <string>


class SubstitutionMatrix
{
    public:
	~SubstitutionMatrix(){}

	static SubstitutionMatrix& getInstance(std::string = "pam250", size_t = prrn::nAACode, std::string = "gap");
	prrn::SCORE getOpenPenalty(size_t = 0) const;
	prrn::SCORE getExtensionPenalty(size_t = 0) const;
	prrn::SCORE getDifferencePenalty(size_t) const;
	prrn::SCORE getMinPenalty(size_t) const;
	prrn::SCORE getScore(prrn::RESCODE, prrn::RESCODE) const;
	size_t getMatSize() const;
	size_t getGapSize() const;
	size_t getNumGapThreshold() const;
	const std::vector<prrn::SCORE>& getGapOpen() const;
	const std::vector<prrn::SCORE>& getGapExtension() const;
	size_t getGapThreshold(size_t) const;

    private:
	SubstitutionMatrix();	// not defined
	SubstitutionMatrix(const std::string&, size_t, const std::string&);
	SubstitutionMatrix(const SubstitutionMatrix&);	// not defined
	SubstitutionMatrix operator=(const SubstitutionMatrix&);	// not defined
	void setMatrix(const std::string&);
	void setPenalty(const std::string&);

	size_t nrt_;
	Matrix<prrn::SCORE> matrix_;
	std::vector<prrn::SCORE> go_;
	std::vector<prrn::SCORE> ge_;
	std::vector<size_t> threshold_;
	static SubstitutionMatrix* pinstance_;
};

inline prrn::SCORE SubstitutionMatrix::getScore(prrn::RESCODE a, prrn::RESCODE b) const
{
#ifdef DEBUG
    assert(a < nrt_ && b < nrt_);
#endif
    return matrix_(a, b);
}

#endif
