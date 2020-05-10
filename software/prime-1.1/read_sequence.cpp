// read_sequence.cpp
//
// Last Modified: 20, Apr 2007
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
#include <fstream>
#include <algorithm>

using namespace std;

void readSequence(const string& fName, vector<Sequence>& seqs)
{
    ifstream f(fName.c_str());
    if(!f)
    {
	cerr << fName << " can't open." << "\n";
	exit(1);
    }

    string tmp;
    string tName;
    string tSeq;

    if(!getline(f, tmp))
    {
	cerr << fName << " is empty.\n";
	exit(2);
    }
    if('>' != tmp[0])
    {
	cerr << "unexpected character: `" << tmp[0] << "'\n";
	exit(2);
    }
    else
	tName.assign(tmp.begin()+1, tmp.end());

    while(getline(f, tmp))
    {
	if('>' != tmp[0])
	    tSeq += tmp;
	else
	{
	    seqs.push_back(Sequence(tName, tSeq));
	    tSeq.clear();
	    tName.assign(tmp.begin()+1, tmp.end());
	}
    }
    seqs.push_back(Sequence(tName, tSeq));
}
