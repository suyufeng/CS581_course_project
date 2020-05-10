// sequence.cpp
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

#include "sequence.h"
#include <iostream>
#ifdef DEBUG
#include <cassert>
#endif

using namespace std;
using namespace prrn;

Sequence::Sequence(const string& name, const string& seq)
    : name_(name)
{
    sequence_.reserve(seq.size());
    convert2rescode(seq.begin(), seq.end(), sequence_);
}

void Sequence::assign(const_iterator start, const_iterator stop)
{
    sequence_.assign(start, stop);
}

void Sequence::concatenate(const Sequence& s)
{
#ifdef DEBUG
    assert(name_ == s.name_
	    || "Sequence to be concatenated must be the same name.");
#endif
    sequence_.insert(sequence_.end(),
	    s.sequence_.begin(), s.sequence_.end());
}

void Sequence::replace(iterator start, iterator stop, const Sequence& target)
{
#ifdef DEBUG
    assert(name_ == target.name_
	    || "Sequence to be concatenated must be the same name.");
#endif
    sequence_.insert(sequence_.erase(start, stop),
	    target.begin(), target.end());
}

void Sequence::clear()
{
    name_.clear();
    sequence_.clear();
}

void Sequence::setName(const string& name)
{
    name_ = name;
}

const string& Sequence::getName() const
{
    return name_;
}

string Sequence::getSequence() const
{
    string s(sequence_.size(), '\0');
    convert2string(sequence_.begin(), sequence_.end(), s.begin());
    return s;
}

size_t Sequence::getLength() const
{
    return sequence_.size();
}

string Sequence::getSegment(size_t start, size_t length) const
{
#ifdef DEBUG
    assert(start+length-1 < sequence_.size());
#endif
    string s(length, '\0');
    vector<RESCODE>::const_iterator begin = sequence_.begin() + start;
    convert2string(begin, begin+length, s.begin());
    return s;
}

Sequence::iterator Sequence::begin()
{
    return sequence_.begin();
}

Sequence::const_iterator Sequence::begin() const
{
    return sequence_.begin();
}

Sequence::iterator Sequence::end()
{
    return sequence_.end();
}

Sequence::const_iterator Sequence::end() const
{
    return sequence_.end();
}

Sequence::reverse_iterator Sequence::rbegin()
{
    return sequence_.rbegin();
}

Sequence::const_reverse_iterator Sequence::rbegin() const
{
    return sequence_.rbegin();
}

Sequence::reverse_iterator Sequence::rend()
{
    return sequence_.rend();
}

Sequence::const_reverse_iterator Sequence::rend() const
{
    return sequence_.rend();
}
