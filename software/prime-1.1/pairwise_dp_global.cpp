// pairwise_dp_global.cpp
//
// Last Modified: 10, May 2007
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

#include "pairwise_dp_algorithm.h"
#include "substitution_matrix.h"
#include "sequence.h"
#include "alignment.h"
#include "path.h"
#include <iostream>
#ifdef MSA_MATRIX
#include <iomanip>
#endif

using namespace std;
using namespace prrn;

PairwiseDPGlobal::PairwiseDPGlobal()
    : PairwiseDPAlgorithm()
    , ge_(sub_.getExtensionPenalty())
{}

void PairwiseDPGlobal::showAlignmentMode() const
{
    cerr << "pairwise alignment mode is global (terminal gaps are penalized)\n";
    cerr << "Gap penalty function g(x) = " << ge_ << "x + " << go_ << "\n\n";
}

namespace
{
    typedef char pathdir;
    const pathdir fromDiagonal = 0;
    const pathdir fromVertical = 1;
    const pathdir fromHorizontal = 2;
}

Alignment PairwiseDPGlobal::getAlignment(const Sequence &a, const Sequence &b) const
{
    // sequence lengths
    const size_t I = a.getLength();
    const size_t J = b.getLength();

    // scores
    vector<SCORE> scr_max(J+1, limitScoreZero);		// maximum score
    vector<SCORE> scr_vert(J+1, - limitScoreMax);	// vertical score
    SCORE tmp_scr_diag = limitScoreZero;	// temporary diagonal score
    SCORE tmp_scr_hori = limitScoreZero;	// temporary horizontal score

    // for optimal path
    PathList path_list(I+J);
    vector<size_t> path_max(J+1, 0);
    vector<size_t> path_vert(J+1, 0);
    size_t tmp_path_diag = 0;
    size_t tmp_path_hori = 0;

    // for path direction
    vector<pathdir> path_direction(J+1, fromHorizontal);
    pathdir diag_path_direction = fromDiagonal;

    // initialize
    scr_max[0] -= go_;
    for(size_t j = 1; j <= J; ++j)
	scr_max[j] = scr_max[j-1] - ge_;
    path_list.setPath(0, 0, 0);

#ifdef MATRIX
    const RESIDUE* rescode = getResCode(res_type);
    cout << setw(6) << " ";
    for(size_t j = 0; j < J; ++j)
	cout << setw(3) << rescode[b.getResidueCode(j)] << " ";
    cout << "\n";	
    cout << setw(2) << " ";
    for(size_t j = 0; j <= J; ++j)
	cout << setw(3) << scr_max[j] << " ";
    cout << "\n";
#endif

    // alignment
    for(size_t i = 1; i <= I; ++i)
    {
	// initialize
	scr_max[0] -= ge_;
	tmp_scr_hori = - limitScoreMax;

	tmp_path_diag = tmp_path_hori = 0;

	// calculation
	for(size_t j = 1; j <= J; ++j)
	{
	    // diagonal
	    SCORE t_max = tmp_scr_diag + sub_.getScore(a.getResidueCode(i-1), b.getResidueCode(j-1));
	    size_t t_path = tmp_path_diag;
	    pathdir t_path_direction = fromDiagonal;

	    // vertical
	    if(scr_max[j] - go_ > scr_vert[j])
	    {
		scr_vert[j] = scr_max[j] - go_;
		path_vert[j] = path_max[j];
	    }
	    scr_vert[j] -= ge_;
	    if(t_max < scr_vert[j])
	    {
		t_max = scr_vert[j];
		t_path = path_vert[j];
		t_path_direction = fromVertical;
	    }

	    // horizontal
	    if(scr_max[j-1] - go_ > tmp_scr_hori)
	    {
		tmp_scr_hori = scr_max[j-1] - go_;
		tmp_path_hori = path_max[j-1];
	    }
	    tmp_scr_hori -= ge_;
	    if(t_max < tmp_scr_hori)
	    {
		t_max = tmp_scr_hori;
		t_path = tmp_path_hori;
		t_path_direction = fromHorizontal;
	    }

	    // set path
	    tmp_path_diag = path_max[j];
	    if(t_path_direction == fromDiagonal && diag_path_direction != fromDiagonal)
	    {
		path_list.setPath(i-1, j-1, t_path);
		path_max[j] = path_list.getSize() - 1;
	    }
	    else
		path_max[j] = t_path;

	    // set path direction
	    diag_path_direction = path_direction[j];
	    path_direction[j] = t_path_direction;

	    // set score
	    tmp_scr_diag = scr_max[j];
	    scr_max[j] = t_max;
	}
#ifdef MATRIX
	cout << rescode[a.getResidueCode(i-1)] << " ";
	for(size_t j = 0; j <= J; ++j)
	    cout << setw(3) << scr_max[j] << " ";
	cout << "\n";
#endif

	// initialize
	tmp_scr_diag = scr_max[0];
	diag_path_direction = fromVertical;
    }
    path_list.setPath(I, J, path_max[J]);

#ifdef PSA_TEST
    cout << setw(4) << scr_max[J] << " : ";
#endif

    return Alignment(a, b, path_list);
}
