// pairwise_dp_long.cpp
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
#include <algorithm>
#ifdef MSA_MATRIX
#include <iomanip>
#endif

using namespace std;
using namespace prrn;

PairwiseDPLong::PairwiseDPLong()
    : PairwiseDPAlgorithm()
    , L(sub_.getGapSize())
    , go_(sub_.getGapOpen())
    , ge_(sub_.getGapExtension())
{}

void PairwiseDPLong::showAlignmentMode() const
{
    cerr << "Pairwise alignment mode is long (terminal gaps are penalized)\n";
    cerr << "Gap penalty function g(x) = min{";
    for(size_t i = 0; i < L; ++i)
    {
	cerr << ge_[i] << "x + ";
	cerr << go_[i];
	if(i != L-1)
	    cerr << ", ";
    }
    cerr <<"}\n\n";
}

namespace
{
    typedef char pathdir;
    const pathdir fromDiagonal = 0;
    const pathdir fromVertical = 1;
    const pathdir fromHorizontal = 2;
}

Alignment PairwiseDPLong::getAlignment(const Sequence &a, const Sequence &b) const
{
    // sequence lengths
    const size_t I = a.getLength();
    const size_t J = b.getLength();

    // scores
    // long gap assumed
    vector<SCORE> score_max(J+1, limitScoreZero);	// maximum score
    vector<size_t> path_max(J+1, 0);	// path for maximum score

    Matrix<SCORE> scores_vert(L, J+1);	// vertical scores
    Matrix<size_t> paths_vert(L, J+1);	// paths for vertical score

    vector<SCORE> tmp_scores_hori(L, limitScoreZero);
    vector<size_t> tmp_paths_hori(L, 0);

    SCORE tmp_score_diag = limitScoreZero;	// temporary diagonal scores
    size_t tmp_path_diag = 0;	// temporary path for diagonal score

    // for optimal path
    PathList path_list(I+J);
    path_list.setPath(0, 0, 0);

    // for path direction
    vector<pathdir> path_direction(J+1, fromHorizontal);
    pathdir diag_path_direction = fromDiagonal;

    // for maximum score and its path
    struct _max{
	SCORE score;
	size_t path;
	pair<size_t, size_t> pos;
    } maximum;

    // initialize
    // for maximum
    for(size_t j = 1; j <= J; ++j)
	score_max[j] = - sub_.getMinPenalty(j);

    // for vertical
    for(size_t i = 0; i < L; ++i)
	for(size_t j = 0; j <= J; ++j)
	    scores_vert(i, j) = - limitScoreMax;

#ifdef MATRIX
    const RESIDUE* rescode = getResCode(res_type);
    cout << setw(7) << " ";
    for(size_t j = 0; j < J; ++j)
	cout << setw(4) << rescode[b.getResidueCode(j)] << " ";
    cout << "\n";	
    cout << setw(2) << " ";
    for(size_t j = 0; j <= J; ++j)
	cout << setw(4) << score_max[j] << " ";
    cout << "\n";
#endif

    // alignment
    for(size_t i = 1; i <= I; ++i)
    {
	// initialize
	tmp_score_diag = score_max[0];

	score_max[0] = - sub_.getMinPenalty(i);
	fill(tmp_scores_hori.begin(), tmp_scores_hori.end(), - limitScoreMax);

	tmp_path_diag = 0;
	fill(tmp_paths_hori.begin(), tmp_paths_hori.end(), 0);

	// calculation
	for(size_t j = 1; j <= J; ++j)
	{
	    // diagonal
	    SCORE t_max = tmp_score_diag + sub_.getScore(a.getResidueCode(i-1), b.getResidueCode(j-1));
	    size_t t_path = tmp_path_diag;
	    pathdir t_path_direction = fromDiagonal;

	    // vertical
	    for(size_t l = 0; l < L; ++l)
	    {
		if(score_max[j] - go_[l] > scores_vert(l, j))
		{
		    scores_vert(l, j) = score_max[j] - go_[l];
		    paths_vert(l, j) = path_max[j];
		}
		scores_vert(l, j) -= ge_[l];
		if(t_max < scores_vert(l, j))
		{
		    t_max = scores_vert(l, j);
		    t_path = paths_vert(l, j);
		    t_path_direction = fromVertical;
		}
	    }

	    // horizontal
	    for(size_t l = 0; l < L; ++l)
	    {
		if(score_max[j-1] - go_[l] > tmp_scores_hori[l])
		{
		    tmp_scores_hori[l] = score_max[j-1] - go_[l];
		    tmp_paths_hori[l] = path_max[j-1];
		}
		tmp_scores_hori[l] -= ge_[l];
		if(t_max < tmp_scores_hori[l])
		{
		    t_max = tmp_scores_hori[l];
		    t_path = tmp_paths_hori[l];
		    t_path_direction = fromHorizontal;
		}
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
	    tmp_score_diag = score_max[j];
	    score_max[j] = t_max;
	}
#ifdef MATRIX
	cout << rescode[a.getResidueCode(i-1)] << " ";
	for(size_t j = 0; j <= J; ++j)
	    cout << setw(4) << score_max[j] << " ";
	cout << "\n";
#endif
	// initialize
	diag_path_direction = fromVertical;
    }
    path_list.setPath(I, J, path_max[J]);

#ifdef PSA_TEST
    cout << setw(4) << score_max[J] << " : ";
#endif
    return Alignment(a, b, path_list);
}
