// alignment.h
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

#ifndef __ALIGNMENT_H__
#define __ALIGNMENT_H__

#include "residue.h"
#include "sequence.h"
#include <vector>
#include <map>

class PathList;
class PTHolder;

class Alignment
{
    public:
	Alignment(){}

	Alignment(const Sequence&);

	Alignment(const std::vector<Sequence>&);

	Alignment(const Sequence&, const Sequence&, const PathList&);
	Alignment(const Alignment&, const Alignment&, const PathList&,
		const PTHolder* = 0);

	Alignment(const Alignment&, size_t, size_t);

	~Alignment(){}

	void divide(const std::vector<size_t>&, Alignment&);
	void concatenate(const Alignment&);
	void concatenate(const std::vector<Sequence>&);

	void replace(const Alignment&, size_t, size_t);

	void eraseNullColumn();
	void changeTerminalGap();
	void resetTerminalGap();

	size_t getNumSequences() const;

	size_t getLength() const;

	const std::string& getName(size_t) const;
	const Sequence& getSequence(size_t) const;

	prrn::RESCODE getResidueCode(size_t, size_t) const;

	static void showMSF(const Alignment&, std::ostream&,
		const std::vector<Sequence>* = 0, size_t = 50);
	static void showFASTA(const Alignment&, std::ostream&,
		const std::vector<Sequence>* = 0, size_t = 60);
	static void showPHYLIP(const Alignment&, std::ostream&,
		const std::vector<Sequence>* = 0, size_t = 50);
	static void showGDE(const Alignment&, std::ostream&,
		const std::vector<Sequence>* = 0, size_t = 60);

    private:
	std::vector<Sequence> seqs_;
	static std::map<std::string, size_t> midx_;
};

enum opfmt {FASTA, MSF, PHYLIP, GDE};

double getSequenceIdentity(const Alignment&, size_t, size_t, bool = false);
double getAverageIdentity(const Alignment&, size_t , bool = false);
double getAverageIdentity(const Alignment&, const std::vector<size_t>&, bool = false);

inline prrn::RESCODE Alignment::getResidueCode(size_t seqno, size_t colno) const
{
#ifdef DEBUG
    return seqs_.at(seqno).getResidueCode(colno);
#else
    return seqs_[seqno].getResidueCode(colno);
#endif
}

#endif
