// prrn.cpp
//
// Last Modified: 11, Jul 2007
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

#include "params.h"
#include "pt_holder.h"
#include "pairwise_dp_algorithm.h"
#include "distance_matrix.h"
#include "profile.h"
#include "profile_dp_global.h"
#include "profile_dp_long.h"
#include "res_conserv.h"
#include "anchor.h"
#include <iostream>
#include <iomanip>
#include "string.h"
using namespace std;
using namespace prrn;

template<typename C, typename T>
void showAlignment(Alignment& a, opfmt opformat, basic_ostream<C, T>& bos)
{
    switch(opformat)
    {
	case FASTA:
	    Alignment::showFASTA(a, bos, 0, 60);
	    break;
	case MSF:
	    Alignment::showMSF(a, bos, 0, 60);
	    break;
	case PHYLIP:
	    Alignment::showPHYLIP(a, bos, 0, 50);
	    break;
	case GDE:
	    Alignment::showGDE(a, bos, 0, 60);
	    break;
    }
}

int main(int argc, char **argv)
{
    cerr << "PRIME version " << version << " (coded by " << author << ")\n";

    if(argc < 2)
    {
	showUsage(argv);
	cerr << "\nFor help, type `" << argv[0] << " --help'\n";
	return -1;
    }

    if(strcmp("--help", argv[1]) == 0)
    {
	showHelp(argv);
	showLicense();
	return 0;
    }
    setParameters(argc, argv);

    vector<Sequence> s;
    readSequence(ip_name, s);
    if(s.size() == 1)
    {
	cerr << ip_name << " contains only one sequence.\n";
	return 1;
    }
    showParameters(s.size());

    static_cast<void>(SubstitutionMatrix::getInstance(sm_name, nResType, gp_name));

    PairwiseDPAlgorithm* pda = 0;
    switch(psaa)
    {
	case SeqLong:
	    pda = new PairwiseDPLong;
	    break;
	case SeqAffine:
	    pda = new PairwiseDPGlobal;
	    break;
	default:
	    cerr << "Unsupported mode has been specified.\n";
	    return 1;
    }
    pda->showAlignmentMode();

    if(s.size() == 2)
    {
	cerr << "Calculating pairwise alignment" << flush;
	Alignment a = pda->getAlignment(s[0], s[1]);
	cerr << ": done.\n";
	if(op_name.empty())
	    showAlignment(a, opformat, cout);
	else
	{
	    out.open(op_name.c_str(), ios::out);
	    showAlignment(a, opformat, out);
	}
	delete pda;
	return 0;
    }

    MultipleDPAlgorithm* mda = 0;
    switch(msaa)
    {
	case ProfLong:
	    Profile::setResidueUseRange(PenaltyGap+1, nResType);
	    mda = new ProfileDPLong;
	    break;
	case ProfAffine:
	    Profile::setResidueUseRange(0, nResType);
	    mda = new ProfileDPGlobal;
	    break;
	default:
	    cerr << "Unsupported mode has been specified.\n";
	    exit(1);
    }
    mda->showAlignmentMode();

    Matrix<size_t> dist_mat(s.size(), s.size());
    switch(idm)
    {
	case Psa:
	    calculateDistanceMatrix(s, *pda, dist_mat, isVerbose);
	    break;
	case Oligo:
	    calculateDistanceMatrix(s, dist_mat, isVerbose);
	    break;
    }

    if(isVerbose)
	cerr << "\n===== constructing guide tree =====\n";
    else
	cerr << "constructing guide tree " << flush;
    PTHolder pth(s, dist_mat, eqfactor);
    if(isVerbose)
    {
	pth.show();
	pth.showNode();
	if(pth.isWeightMode())
	{
	    cerr << "\nEdge Weights\n";
	    pth.showWeight();
	    cerr << "\nPair Weights\n";
	    pth.showPairWeight();
	}
    }
    else
	cerr << "... done.\n";

    Alignment a;
    getProgressiveAlignment(s, *pda, *mda, pth, a, isVerbose);

    for(int i = 0; i < nreconst; ++i)
    {
	calculateDistanceMatrix(a, dist_mat, isVerbose);
	cerr << "reconstructing guide tree " << flush;
	pth.construct(dist_mat);
	cerr << " ... done.\n";
	getProgressiveAlignment(s, *pda, *mda, pth, a);
    }

    vector<int> branch_order;
    branch_order.reserve(2*s.size() - 3);
    getBranchOrder(2*s.size()-3, branch_order);
    eraseConservedBranch(a, pth.getPhylogeneticTree(), branch_order, sfbundl_thr, isVerbose);

    vector<pair<size_t, size_t> > anchorPoint;
    vector<pair<size_t, size_t> > nonConsRegion;
    vector<Alignment> partialMSA;
    if(ucfix == UCR || ucfix == UCRS)
    {
	// for anchors based on progressive alignment
	//calcBackground(s, nResType);
	//calcAnchorPoint(a, anchorPoint, anchoring);

	calcNonConservedRegion(a.getLength(), anchorPoint, nonConsRegion);
	decomposeAlignment(a, nonConsRegion, partialMSA);
#ifdef DEBUG
	assert(nonConsRegion.size() == partialMSA.size());
#endif
    }
    // for unchanged region
    Alignment prevA(a);
    for(int i = 0; i < outitr; ++i)
    {
	calculateDistanceMatrix(a, dist_mat);
	cerr << "reconstructing guide tree " << flush;
	pth.construct(dist_mat);
	cerr << " ... done.\n";
	if(isVerbose)
	{
	    pth.showNode();
	    cerr << "\nEdge Weights\n";
	    pth.showWeight();
	}

	random_shuffle(branch_order.begin(), branch_order.end());

	cerr << "\n===== refine the alignment iteratively =====\n";
	const SCORE startscore = mda->getSPScore(a, &pth);
	cerr << "start score = " << setprecision(8) << startscore << endl;

	if(ucfix == NoUC || ucfix == UCS)
	    refineAlignment(*mda, inritr, pth, branch_order, a);
	else
	{
	    vector<pair<size_t, size_t> >::reverse_iterator k = nonConsRegion.rbegin();
	    const vector<Alignment>::reverse_iterator je = partialMSA.rend();
	    for(vector<Alignment>::reverse_iterator j = partialMSA.rbegin(); j != je; ++j, ++k)
	    {
		const SCORE before = mda->getSPScore(*j, &pth);
		refineAlignment(*mda, inritr, pth, branch_order, *j);
		if(mda->getSPScore(*j, &pth) - before > 0.0)
		    a.replace(*j, k->first, k->second);
	    }
	}

	const SCORE refratio = fabs((mda->getSPScore(a, &pth) - startscore) / startscore) * 100.0;
	cerr << "refinement ratio: " << setprecision(3) << refratio << " %\n";
	if(refratio < cutoff)
	    break;

	if(ucfix != NoUC)
	{
	    if(ucfix != UCR)
		eraseUnchangedBranch(prevA, a, pth.getPhylogeneticTree(), branch_order);

	    if(ucfix != UCS)
	    {
		anchorPoint.clear();
		detectUnchangedRegion(prevA, a, anchorPoint);

		nonConsRegion.clear();
		calcNonConservedRegion(a.getLength(), anchorPoint, nonConsRegion);

		partialMSA.clear();
		decomposeAlignment(a, nonConsRegion, partialMSA);
#ifdef DEBUG
		assert(nonConsRegion.size() == partialMSA.size());
#endif
	    }
	    prevA = a;
	}
    }

#ifdef DEBUG
    cerr << "\nPhlogenetic Tree\n";
    pth.show();
    if(pth.isWeightMode())
    {
	cerr << "\nEqualization Factor = " << eqfactor << "\n";
	cerr << "\nEdge Weights\n";
	pth.showWeight();
	cerr << "\nPair Weights\n";
	pth.showPairWeight();
    }
#endif

    cerr << "\nFinal alignment score = " << mda->getSPScore(a, &pth) << "\n";

    if(op_name.empty())
	showAlignment(a, opformat, cout);
    else
    {
	out.open(op_name.c_str(), ios::out);
	showAlignment(a, opformat, out);
    }

    delete pda;
    delete mda;

    return 0;
}
