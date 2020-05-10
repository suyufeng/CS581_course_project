// alignment.cpp
//
// Last Modified: 4, Jun 2007
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

#include "alignment.h"
#include "prrn.h"
#include "path.h"
#include "pt_holder.h"
#include <iostream>
#include <iomanip>
#include <algorithm>
#include <cctype>
#include <cassert>

using namespace std;
using namespace prrn;

map<string, size_t> Alignment::midx_;

Alignment::Alignment(const Sequence& s)
{
    seqs_.push_back(s);
}

Alignment::Alignment(const vector<Sequence>& seqs)
{
#ifdef DEBUG
    const int nseq = seqs.size();
    const size_t len = (nseq) ? seqs[0].getLength() : 0;
    for(int i = 1; i < nseq; ++i)
	assert(seqs[i].getLength() == len
		&& "the vector of Sequence does not have alignment form");
#endif
    seqs_.reserve(seqs.size());
    copy(seqs.begin(), seqs.end(), back_inserter(seqs_));
}

Alignment::Alignment(const Sequence& a, const Sequence& b, const PathList& p)
{
    vector<RESCODE> aln_a;
    vector<RESCODE> aln_b;
    Sequence::const_iterator civr;

    vector<size_t> idxs;
    p.getPathIndexes(idxs);
    vector<size_t>::const_iterator i = idxs.begin();
    const vector<size_t>::const_iterator ie = idxs.end();

    size_t x = p.pathlist_[*i].m;
    size_t y = p.pathlist_[*i].n;
    ++i;
    size_t xx = p.pathlist_[*i].m;
    size_t yy = p.pathlist_[*i].n;

    for(;;)
    {
	size_t len_x = xx - x;
	size_t len_y = yy - y;
	if(len_x == 0)
	{
	    aln_a.insert(aln_a.end(), len_y, PenaltyGap);
	    civr = b.begin() + y;
	    copy(civr, civr+len_y, back_inserter(aln_b));
	}
	else if(len_y == 0)
	{
	    civr = a.begin() + x;
	    copy(civr, civr+len_x, back_inserter(aln_a));
	    aln_b.insert(aln_b.end(), len_x, PenaltyGap);
	}
	else
	{
	    size_t len = min(len_x, len_y);

	    civr = a.begin() + x;
	    copy(civr, civr+len, back_inserter(aln_a));
	    civr = b.begin() + y;
	    copy(civr, civr+len, back_inserter(aln_b));

	    x += len;
	    y += len;
	    continue;
	}

	x = xx;
	y = yy;

	if(++i == ie)
	    break;

	xx = p.pathlist_[*i].m;
	yy = p.pathlist_[*i].n;
    }

#ifdef DEBUG
    assert(aln_a.size() == aln_b.size()
	    && "alignment length error");
#endif
    // do not reserve operation
    seqs_.push_back(Sequence(a.getName(), aln_a));
    seqs_.push_back(Sequence(b.getName(), aln_b));
}

Alignment::Alignment(const Alignment& A, const Alignment& B, const PathList& p,
		const PTHolder* pth)
{
    const int M = static_cast<int>(A.seqs_.size());
    const int N = static_cast<int>(B.seqs_.size());
    const int tsize = M+N;

    vector<int> idx_a;
    idx_a.reserve(M);
    vector<int> idx_b;
    idx_b.reserve(N);

    const vector<Sequence>* pvs = (pth) ? &pth->getSequenceVector() : 0;
    const bool isSort(pvs && pvs->size() == static_cast<size_t>(tsize));
    if(isSort)
    {
	for(int m = 0; m < M; ++m)
	    idx_a.push_back(pth->getIndex(A.getName(m)));
	for(int n = 0; n < N; ++n)
	    idx_b.push_back(pth->getIndex(B.getName(n)));
    }
    else
    {
	for(int m = 0; m < M; ++m)
	    idx_a.push_back(m);
	for(int n = 0; n < N; ++n)
	    idx_b.push_back(M+n);
    }

    vector<vector<RESCODE> > aln(tsize);
    Sequence::const_iterator civr;

    vector<size_t> idxs;
    p.getPathIndexes(idxs);
    const size_t idxsize = idxs.size();

    size_t x = p.pathlist_[idxs[0]].m;
    size_t y = p.pathlist_[idxs[0]].n;

    size_t xx = p.pathlist_[idxs[1]].m;
    size_t yy = p.pathlist_[idxs[1]].n;
    size_t i = 1;

    for(;;)
    {
	size_t len_x = xx - x;
	size_t len_y = yy - y;

	if(len_x == 0)
	{
	    for(int m = 0; m < M; ++m)
		aln[idx_a[m]].insert(aln[idx_a[m]].end(), len_y, PenaltyGap);
	    for(int n = 0; n < N; ++n)
	    {
		civr = B.seqs_[n].begin() + y;
		copy(civr, civr+len_y, back_inserter(aln[idx_b[n]]));
	    }
	}
	else if(len_y == 0)
	{
	    for(int m = 0; m < M; ++m)
	    {
		civr = A.seqs_[m].begin() + x;
		copy(civr, civr+len_x, back_inserter(aln[idx_a[m]]));
	    }
	    for(int n = 0; n < N; ++n)
		aln[idx_b[n]].insert(aln[idx_b[n]].end(), len_x, PenaltyGap);
	}
	else
	{
	    size_t len = min(len_x, len_y);

	    for(int m = 0; m < M; ++m)
	    {
		civr = A.seqs_[m].begin() + x;
		copy(civr, civr+len, back_inserter(aln[idx_a[m]]));
	    }
	    for(int n = 0; n < N; ++n)
	    {
		civr = B.seqs_[n].begin() + y;
		copy(civr, civr+len, back_inserter(aln[idx_b[n]]));
	    }

	    x += len;
	    y += len;
	    continue;
	}

	x = xx;
	y = yy;

	if(++i == idxsize)
	    break;

	xx = p.pathlist_[idxs[i]].m;
	yy = p.pathlist_[idxs[i]].n;
    }

#ifdef DEBUG
    size_t length = aln[0].size();
    for(int i = 1; i < tsize; ++i)
	assert(aln[i].size() == length
		&& "alignment length error");
#endif
    seqs_.reserve(tsize);
    if(isSort)
    {
	for(int i = 0; i < tsize; ++i)
	    seqs_.push_back(Sequence((*pvs)[i].getName(), aln[i]));
    }
    else
    {
	for(int m = 0; m < M; ++m)
	    seqs_.push_back(Sequence(A.getName(m), aln[m]));
	for(int n = 0; n < N; ++n)
	    seqs_.push_back(Sequence(B.getName(n), aln[M+n]));
    }
}

Alignment::Alignment(const Alignment& a, size_t begin, size_t end)
{
#ifdef DEBUG
    assert(begin < end && end < a.getLength());
#endif
    seqs_.reserve(a.seqs_.size());
    for(vector<Sequence>::const_iterator i = a.seqs_.begin(), ie = a.seqs_.end(); i != ie; ++i)
    {
	seqs_.push_back(Sequence(i->getName())); 
	seqs_.back().append(i->begin()+begin, i->begin()+end+1);
    }
}

void Alignment::divide(const vector<size_t>& vs, Alignment& ra)
{
    const size_t vnseq = vs.size();
    ra.seqs_.clear();
    ra.seqs_.reserve(vnseq);

#ifdef DEBUG
    if(!vnseq)
	return;
    if(seqs_.size() == vnseq)
    {
	seqs_.swap(ra.seqs_);
	return;
    }
#endif

    vector<size_t>::const_iterator j = vs.begin();
    const vector<size_t>::const_iterator je = vs.end();
    const size_t tnseq = seqs_.size();
    vector<Sequence> temp;
    temp.reserve(tnseq-vnseq);
    for(size_t i = 0; i < tnseq; ++i)
    {
	if(j != je && i == *j)
	    ra.seqs_.push_back(seqs_[*(j++)]);
	else
	    temp.push_back(seqs_[i]);
    }
    seqs_.swap(temp);

    eraseNullColumn();
    ra.eraseNullColumn();
}

void Alignment::concatenate(const vector<Sequence>& vs)
{
    if(seqs_.empty())
    {
	seqs_.reserve(vs.size());
	copy(vs.begin(), vs.end(), back_inserter(seqs_));
	return;
    }
    const size_t nrow = seqs_.size();
#ifdef DEBUG
    assert(nrow == vs.size()
	    || "Alignment to be concatenated must have the same number of rows.");
#endif
    vector<size_t> idx;
    idx.reserve(nrow);
    for(size_t i = 0; i < nrow; ++i)
	idx.push_back(i);

    for(size_t i = 0; i < nrow; ++i)
    {
	for(vector<size_t>::iterator j = idx.begin(); j != idx.end(); ++j)
	{
	    if(seqs_[i].getName() == vs[*j].getName())
	    {
		seqs_[i].concatenate(vs[*j]);
		idx.erase(j);
		break;
	    }
	}
    }
}
void Alignment::concatenate(const Alignment& target)
{
    concatenate(target.seqs_);
}

void Alignment::replace(const Alignment& target, size_t start, size_t stop)
{
    const size_t nrow = seqs_.size();
#ifdef DEBUG
    assert(nrow == target.seqs_.size()
	    || "Alignment to be concatenated must have the same number of rows.");
    assert(start < stop);
    assert(stop < getLength());
#endif
    vector<size_t> idx;
    idx.reserve(nrow);
    for(size_t i = 0; i < nrow; ++i)
	idx.push_back(i);

    Sequence::iterator offset;
    for(size_t i = 0; i < nrow; ++i)
    {
	for(vector<size_t>::iterator j = idx.begin(); j != idx.end(); ++j)
	{
	    if(seqs_[i].getName() == target.seqs_[*j].getName())
	    {
		offset = seqs_[i].begin();
		seqs_[i].replace(offset+start, offset+stop+1, target.seqs_[*j]);
		idx.erase(j);
		break;
	    }
	}
    }
}

void Alignment::eraseNullColumn()
{
#ifdef DEBUG
    assert(!seqs_.empty());
#endif
    const size_t nseq = seqs_.size();

    if(nseq == 1)
    {
	seqs_[0].sequence_.erase(remove(seqs_[0].begin(),
		    seqs_[0].end(), PenaltyGap),
		seqs_[0].end());
	return;
    }

    const size_t length = seqs_[0].getLength();

    typedef char Bool;
    const Bool True = 1;
    const Bool False = 0;
    vector<Bool> isGapColumn(length, True);

    size_t ntrue = 0;
    const vector<Sequence>::iterator ib = seqs_.begin();
    const vector<Sequence>::iterator ie = seqs_.end();
    for(vector<Sequence>::iterator i = ib; i != ie; ++i)
    {
	for(size_t j = 0; j < length; ++j)
	{
	    if(isGapColumn[j] == True && isResidue(i->getResidueCode(j)))
	    {
		isGapColumn[j] = False;
		++ntrue;
	    }
	}
    }
    if(length == ntrue)
	return;

    for(vector<Sequence>::iterator i = ib; i != ie; ++i)
    {
	for(size_t j = 0; j < length; ++j)
	    if(isGapColumn[j] == True)
		i->sequence_[j] = limitRescodeMax;
	i->sequence_.erase(remove(i->begin(), i->end(),
		    limitRescodeMax), i->end());
    }
}

void Alignment::changeTerminalGap()
{
    const vector<Sequence>::iterator send = seqs_.end();
    for(vector<Sequence>::iterator i = seqs_.begin(); i != send; ++i)
    {
	fill(i->begin(),
		find_if(i->begin(), i->end(),
		    bind2nd(greater<RESCODE>(), PenaltyGap)),
		nonPenaltyGap);
	fill(i->rbegin(),
		find_if(i->rbegin(), i->rend(),
		    bind2nd(greater<RESCODE>(), PenaltyGap)),
		nonPenaltyGap);
    }
}

void Alignment::resetTerminalGap()
{
    using std::replace;
    const vector<Sequence>::iterator send = seqs_.end();
    for(vector<Sequence>::iterator i = seqs_.begin(); i != send; ++i)
    {
	replace(i->begin(),
		find_if(i->begin(), i->end(),
		    bind2nd(greater<RESCODE>(), PenaltyGap)),
		nonPenaltyGap, PenaltyGap);
	replace(i->rbegin(),
		find_if(i->rbegin(), i->rend(),
		    bind2nd(greater<RESCODE>(), PenaltyGap)),
		nonPenaltyGap, PenaltyGap);
    }
}

size_t Alignment::getNumSequences() const
{
    return seqs_.size();
}

size_t Alignment::getLength() const
{
#ifdef DEBUG
    return seqs_.at(0).getLength();
#else
    return seqs_[0].getLength();
#endif
}

const string& Alignment::getName(size_t n) const
{
#ifdef DEBUG
    return seqs_.at(n).getName();
#else
    return seqs_[n].getName();
#endif
}

const Sequence& Alignment::getSequence(size_t n) const
{
#ifdef DEBUG
    return seqs_.at(n);
#else
    return seqs_[n];
#endif
}

namespace
{
    void setIndeces(const Alignment& a, const vector<Sequence>* vs,
	    map<string, size_t>& idx_map, vector<size_t>& ref_idx)
    {
	if(!vs)
	{
	    const size_t rsize = ref_idx.size();
	    for(size_t i = 0; i < rsize; ++i)
		ref_idx[i] = i;
	}
	else
	{
	    idx_map.clear();
	    const size_t vsize = vs->size();
	    for(size_t i = 0; i < vsize; ++i)
		idx_map[(*vs)[i].getName()] = i;

	    for(size_t i = 0; i < vsize; ++i)
		ref_idx[i] = idx_map[a.getName(i)];
	}
    }
}

void Alignment::showMSF(const Alignment& a, ostream& os,
	const vector<Sequence>* vs, size_t row_size)
{
#ifdef DEBUG
    assert(!a.seqs_.empty());
    if(vs)
	assert(vs->size() == static_cast<size_t>(a.seqs_.size()));
#endif

    const size_t ssize = a.seqs_.size();
    vector<size_t> ref_index(ssize);
    setIndeces(a, vs, midx_, ref_index);

    size_t name_length = a.getName(ref_index[0]).size();
    for(size_t i = 1; i < ssize; ++i)
	name_length = max(name_length, a.getName(ref_index[i]).size());

    // show preamble
    const size_t length = a.seqs_[0].getLength();
    os << "PileUp\n\n\n\n"
	<< "   MSF: ";
    os.width(4);
    os << length << "  Type: P    Check:  0000   ..\n\n";
    for(size_t i = 0; i < ssize; ++i)
    {
	os << " Name: ";
	os.width(name_length);
	os << left << a.getName(ref_index[i]) << " oo  Len: ";
	os.width(4);
	os << right << length << "  Check:  0000  Weight:  1.000\n";
    }
    os << "\n//\n\n\n";

    // show alignment
    const size_t window = 10;
    for(size_t i = 0; i < length; i += row_size)
    {
	os << "\n";
	const size_t trs = min(row_size, length-i);
	for(size_t j = 0; j < ssize; ++j)
	{
	    os.width(name_length);
	    os << left << a.getName(ref_index[j]);
	    os.width(5);
	    os << " ";
	    for(size_t k = 0; k < trs; k += window)
	    {
		os << " " << a.seqs_[ref_index[j]].getSegment(i+k, min(window, length-i-k));
	    }
	    os << "\n";
	}
	os << "\n";
    }
}

void Alignment::showFASTA(const Alignment& a, ostream& os,
	const vector<Sequence>* vs, size_t row_size)
{
#ifdef DEBUG
    assert(!a.seqs_.empty());
    if(vs)
	assert(vs->size() == static_cast<size_t>(a.seqs_.size()));
#endif

    const size_t ssize = a.seqs_.size();
    vector<size_t> ref_index(ssize);
    setIndeces(a, vs, midx_, ref_index);

    const size_t length = a.seqs_[0].getLength();
    for(size_t i = 0; i < ssize; ++i)
    {
	os << ">" << a.getName(ref_index[i]) << "\n";
	for(size_t j = 0; j < length; j += row_size)
	    os << a.seqs_[ref_index[i]].getSegment(j, min(row_size, length-j)) << "\n";
    }
}

void Alignment::showPHYLIP(const Alignment& a, ostream& os,
	const vector<Sequence>* vs, size_t row_size)
{
#ifdef DEBUG
    assert(!a.seqs_.empty());
    if(vs)
	assert(vs->size() == static_cast<size_t>(a.seqs_.size()));
#endif

    const size_t ssize = a.seqs_.size();
    vector<size_t> ref_index(ssize);
    setIndeces(a, vs, midx_, ref_index);

    size_t name_length = a.getName(ref_index[0]).size();
    for(size_t i = 1; i < ssize; ++i)
	name_length = max(name_length, a.getName(ref_index[i]).size());

    // show preamble
    os << setw(6) << ssize << " " << setw(6) << a.getLength();

    // show alignment
    const size_t window = 10;
    const size_t length = a.seqs_[0].getLength();
    for(size_t i = 0; i < length; i += row_size)
    {
	os << "\n";
	const size_t trs = min(row_size, length-i);
	for(size_t j = 0; j < ssize; ++j)
	{
	    os.width(name_length);
	    if(i == 0)
		os << left << a.getName(ref_index[j]);
	    else
		os << " ";
	    for(size_t k = 0; k < trs; k += window)
		os << " " << a.seqs_[ref_index[j]].getSegment(i+k, min(window, length-i-k));
	    os << "\n";
	}
    }
}

void Alignment::showGDE(const Alignment& a, ostream& os,
	const vector<Sequence>* vs, size_t row_size)
{
#ifdef DEBUG
    assert(!a.seqs_.empty());
    if(vs)
	assert(vs->size() == static_cast<size_t>(a.seqs_.size()));
#endif

    const size_t ssize = a.seqs_.size();
    vector<size_t> ref_index(ssize);
    setIndeces(a, vs, midx_, ref_index);

    const size_t length = a.seqs_[0].getLength();
    string tmpstr;
    for(size_t i = 0; i < ssize; ++i)
    {
	os << "%" << a.getName(ref_index[i]) << "\n";
	for(size_t j = 0; j < length; j += row_size)
	{
	    tmpstr = a.seqs_[ref_index[i]].getSegment(j, min(row_size, length-j));
	    const string::iterator ie = tmpstr.end();
	    for(string::iterator i = tmpstr.begin(); i != ie; ++i)
		*i = tolower(*i);
	    os << tmpstr << "\n";
	}
    }
}

double getSequenceIdentity(const Alignment& a, size_t p, size_t q, bool isNeglectTerminal)
{
    const size_t size = a.getLength();
#ifdef DEBUG
    assert(p < size && q < size
	    && "out of range error");
#endif
    size_t start = 0;
    if(isNeglectTerminal)
	for(size_t i = 0; i != size; ++i)
	    if(isResidue(a.getResidueCode(p, i)) && isResidue(a.getResidueCode(q, i)))
	    {
		start = i;
		break;
	    }
    size_t end = size;
    if(isNeglectTerminal)
	for(size_t i = size-1; i != 0; --i)
	    if(isResidue(a.getResidueCode(p, i)) && isResidue(a.getResidueCode(q, i)))
	    {
		end = i;
		break;
	    }

    size_t length = 0;
    size_t nIdentity = 0;


    RESCODE arc = 0;
    RESCODE brc = 0;
    for(size_t i = start; i < end; ++i)
    {
	arc = a.getResidueCode(p, i);
	brc = a.getResidueCode(q, i);
	if(isResidue(arc) || isResidue(brc))
	{
	    ++length;
	    if(arc == brc)
		++nIdentity;
	}
    }
    return (double)nIdentity / (double)length;
}

double getAverageIdentity(const Alignment& a, size_t idx, bool isNeglectTerminal)
{
    double seqid = 0.0;
    const size_t nrow = a.getNumSequences();
    for(size_t i = 0; i < nrow; ++i)
	seqid += (i != idx) ? getSequenceIdentity(a, idx, i, isNeglectTerminal) : 0.0;
    return seqid / (double)(nrow-1);
}

double getAverageIdentity(const Alignment& a, const vector<size_t>& idxs, bool isNeglectTerminal)
{
    double seqid = 0.0;
    const size_t nrow = a.getNumSequences();
    const vector<size_t>::const_iterator jb = idxs.begin();
    const vector<size_t>::const_iterator je = idxs.end();
    for(size_t i = 0; i < nrow; ++i)
	if(je == find(jb, je, i))
	    for(vector<size_t>::const_iterator j = jb; j != je; ++j)
		seqid += getSequenceIdentity(a, i, *j, isNeglectTerminal);
    const size_t isize = idxs.size();
    return seqid / (double)(isize * (nrow - isize));
}
