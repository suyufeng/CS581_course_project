// path.h
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

#ifndef __PATH_H__
#define __PATH_H__

#include <vector>

class Alignment;

typedef struct _path
{
    size_t m;
    size_t n;

    _path(size_t p, size_t q)
	: m(p)
	, n(q)
    {}
} Path;

class PathList
{
    public:
	friend class Alignment;
	PathList(size_t);

	void setPath(size_t, size_t, size_t);
	void getPathIndexes(std::vector<size_t>&) const;

	size_t getSize() const;
#ifdef DEBUG
	void show() const;
	void showOptimalPath() const;
#endif
    private:
	std::vector<Path> pathlist_;
	std::vector<size_t> prevlist_;
};

inline void PathList::setPath(size_t m, size_t n, size_t p)
{
    pathlist_.push_back(Path(m, n));
    prevlist_.push_back(p);
}

#endif
