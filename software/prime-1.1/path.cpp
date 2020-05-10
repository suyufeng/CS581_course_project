// path.cpp
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

#include "path.h"
#include <iostream>
#include <iomanip>
#include <algorithm>

using namespace std;

PathList::PathList(size_t allocsize = 0)
{
    if(allocsize != 0)
    {
	pathlist_.reserve(allocsize);
	prevlist_.reserve(allocsize);
    }
}

void PathList::getPathIndexes(vector<size_t>& vidx) const
{
    size_t ptr = pathlist_.size() - 1;

    while(ptr != prevlist_[ptr])
    {
	vidx.push_back(ptr);
	ptr = prevlist_[ptr];
    }
    vidx.push_back(0);
    reverse(vidx.begin(), vidx.end());
}

size_t PathList::getSize() const
{
    return pathlist_.size();
}

#ifdef DEBUG
void PathList::show() const
{
    for(size_t i = 0; i < pathlist_.size(); ++i)
	cout << setw(3) << i
	    << ": (" << setw(3) << pathlist_[i].m
	    << ", " << setw(3) <<  pathlist_[i].n
	    << ") - " << setw(3) << prevlist_[i] << "\n";
}

void PathList::showOptimalPath() const
{
    vector<size_t> vs;
    getPathIndexes(vs);

    const size_t vssize = vs.size();
    for(size_t i = 0; i < vssize; ++i)
	cout << "(" << setw(3) << pathlist_[vs[i]].m << ", "
	    << setw(3) << pathlist_[vs[i]].n << ")\n";
}
#endif
