// sequence.h
//
// Last Modified: 22, Mar 2007
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

#ifndef __SEQUENCE_H__
#define __SEQUENCE_H__

#include "residue.h"
#include <string>
#include <vector>

class Alignment;

class Sequence
{
    public:
	friend class Alignment;
	typedef std::vector<prrn::RESCODE>::iterator iterator;
	typedef std::vector<prrn::RESCODE>::const_iterator const_iterator;
	typedef std::vector<prrn::RESCODE>::reverse_iterator reverse_iterator;
	typedef std::vector<prrn::RESCODE>::const_reverse_iterator const_reverse_iterator;

	Sequence(){}
	Sequence(const std::string&, const std::string&);
	Sequence(const std::string& name)
	    : name_(name)
	    , sequence_()
	    {}
	Sequence(const std::string& name, const std::vector<prrn::RESCODE>& seq)
	    : name_(name)
	    , sequence_(seq)
	    {}
	~Sequence(){}

	void append(prrn::RESCODE r)
	{
	    sequence_.push_back(r);
	}

	void append(const_iterator begin, const_iterator end)
	{
	    sequence_.insert(sequence_.end(), begin, end);
	}

	void assign(const_iterator, const_iterator);

	void concatenate(const Sequence&);
	void replace(iterator, iterator, const Sequence&);

	void clear();

	void setName(const std::string&);
	const std::string& getName() const;

	std::string getSequence() const;
	std::string getSegment(size_t, size_t) const;

	prrn::RESCODE getResidueCode(size_t) const;

	size_t getLength() const;

	iterator begin();
	const_iterator begin() const;
	iterator end();
	const_iterator end() const;

	reverse_iterator rbegin();
	const_reverse_iterator rbegin() const;
	reverse_iterator rend();
	const_reverse_iterator rend() const;

    private:
	std::string name_;
	std::vector<prrn::RESCODE> sequence_;
};

inline prrn::RESCODE Sequence::getResidueCode(size_t pos) const
{
#ifdef DEBUG
    return sequence_.at(pos);
#endif
    return sequence_[pos];
}

void readSequence(const std::string&, std::vector<Sequence>&);

#endif
