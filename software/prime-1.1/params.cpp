// params.cpp
//
// Last Modified: 12, Jul 2007
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
#include <iostream>
#include <iomanip>
#include "string.h"

using namespace std;
using namespace prrn;

namespace
{
    int getHashValue(const char* const str)
    {
	const size_t len = strlen(str);
	return len * (str[0] + str[len-1]);
	/*
	const int Prmno = 211;
	unsigned int hval = 0, g = 0;
	for(char* p = str; *p; ++p)
	{
	    hval = (hval << 4) + (*p);
	    if((g = hval & 0xf0000000) != 0)
	    {
		hval ^= g>>24;
		hval ^= g;
	    }
	}
	return hval % Prmno;
	*/
    }

    const int hip = 1105;	//getHashValue("input");
    const int hop = 1362;	//getHashValue("output");
    const int hfmt = 1308;	//getHashValue("format");
    const int hrt = 2580;	//getHashValue("residue-type");
    const int hsm = 4465;	//getHashValue("substitution-matrix");
    const int hgp = 1752;	//getHashValue("gap-cost");
    const int halg = 1854;	//getHashValue("algorithm");
    const int hdst = 1608;	//getHashValue("distance");
    const int hnr = 3740;	//getHashValue("num-recalculation");
    const int hanc = 1266;	//getHashValue("anchor");
    const int hsfb = 3456;	//getHashValue("subfamily-bundle");
    const int hucf = 3081;	//getHashValue("unchanged-fix");
    const int hpt = 3621;	//getHashValue("phylogenetic-tree");
    const int hni = 4180;	//getHashValue("num-outer-iteration");
    const int hcut = 1206;	//getHashValue("cutoff");
    const int heqf = 4085;	//getHashValue("equalization-factor");

    const char cip[] = "input";
    const char cop[] = "output";
    const char cfmt[] = "format";
    const char crt[] = "residue-type";
    const char csm[] = "substitution-matrix";
    const char cgp[] = "gap-cost";
    const char calg[] = "algorithm";
    const char cdst[] = "distance";
    const char cnr[] = "num-recalculation";
    const char canc[] = "anchor";
    const char csfb[] = "subfamily-bundle";
    const char cucf[] = "unchanged-fix";
    const char cpt[] = "phylogenetic-tree";
    const char cni[] = "num-outer-iteration";
    const char ccut[] = "cutoff";
    const char ceqf[] = "equalization-factor";

    void showHash()
    {
	cerr << "\nHash values:\n"
	    << setw(20) << cip << "\t" << getHashValue(cip) << "\n"
	    << setw(20) << cop << "\t" << getHashValue(cop) << "\n"
	    << setw(20) << cfmt << "\t" << getHashValue(cfmt) << "\n"
	    << setw(20) << crt << "\t" << getHashValue(crt) << "\n"
	    << setw(20) << csm << "\t" << getHashValue(csm) << "\n"
	    << setw(20) << cgp << "\t" << getHashValue(cgp) << "\n"
	    << setw(20) << calg << "\t" << getHashValue(calg) << "\n"
	    << setw(20) << cdst << "\t" << getHashValue(cdst) << "\n"
	    << setw(20) << cnr << "\t" << getHashValue(cnr) << "\n"
	    << setw(20) << canc << "\t" << getHashValue(canc) << "\n"
	    << setw(20) << csfb << "\t" << getHashValue(csfb) << "\n"
	    << setw(20) << cucf << "\t" << getHashValue(cucf) << "\n"
	    << setw(20) << cpt << "\t" << getHashValue(cpt) << "\n"
	    << setw(20) << cni << "\t" << getHashValue(cni) << "\n"
	    << setw(20) << ccut << "\t" << getHashValue(ccut) << "\n"
	    << setw(20) << ceqf << "\t" << getHashValue(ceqf) << "\n";
    }

    void setFormat(const char* fmt)
    {
	if(strcmp("msf", fmt) == 0)
	    opformat = MSF;
	else if(strcmp("phylip", fmt) == 0)
	    opformat = PHYLIP;
	else if(strcmp("gde", fmt) == 0)
	    opformat = GDE;
	else if(strcmp("fasta", fmt) != 0)
	{
	    cerr << "\nArgument of `-f (--format)' option must be `fasta', `msf', `phylip', or `gde'.\n";
	    cerr << "default option `-f fasta' is used.\n\n";
	}
    }

    void setResidueType(const char* rt)
    {
	if(strcmp("amino", rt) == 0)
	{
	    res_type = AA;
	    nResType = nAACode;
	}
	else if(strcmp("dna", rt) == 0)
	{
	    res_type = DNA;
	    nResType = nNACode;
	}
	else if(strcmp("rna", rt) == 0)
	{
	    res_type = RNA;
	    nResType = nNACode;
	}
	else
	{
	    cerr << "\nArgument of `-t (--residue-type)' option must be `amino', `dna', or `rna'.\n";
	    cerr << "default option `-t amino' is used.\n\n";
	}
    }

    void setDistanceMethod(const char* mthd)
    {
	if(strcmp("oligo", mthd) == 0)
	    idm = Oligo;
	else if(strcmp("psa", mthd) != 0)
	{
	    cerr << "\nArgument of `-d (--distance)' option must be either `psa' or `oligo'.\n";
	    cerr << "default option `-d psa' is used.\n\n";
	}
    }

    void setGapCost(const char* gcname)
    {
	if(strcmp("affine", gcname) == 0)
	{
	    psaa = SeqAffine;
	    msaa = ProfAffine;
	}
	else if(strcmp("piecewise", gcname) != 0)
	{
	    cerr << "\nArgument of `-a (--algo)' option must be either `piecewise' or `affine'.\n";
	    cerr << "default option `-a piecewise' is used.\n\n";
	}
    }

    void setAnchorMethod(const char* mthd)
    {
	if(strcmp("nongap", mthd) == 0)
	    anchoring = NonGap;
	else if(strcmp("cons", mthd) == 0)
	    anchoring = Conservation;
	else if(strcmp("none", mthd) != 0)
	{
	    cerr << "\nArgument of `-n (--anchor)' option must be either `none', `nongap' or `cons'.\n";
	    cerr << "default option `-n none' is used.\n\n";
	}
    }

    const char* showAnchorMethod()
    {
	switch(anchoring)
	{
	    case NoAnchor:
		return "None";
	    case Conservation:
		return "Conservation";
	    case NonGap:
	    default:
		return "Non gap";
	}
    }

    void setUnchangedStrategy(const char* strategy)
    {
	if(strcmp("ucr", strategy) == 0)
	    ucfix = UCR;
	else if(strcmp("ucs", strategy) == 0)
	    ucfix = UCS;
	else if(strcmp("ucrs", strategy) == 0)
	    ucfix = UCRS;
	else if(strcmp("nouc", strategy) == 0)
	    ucfix = NoUC;
	else
	{
	    cerr << "\nArgument of `-u (--unchanged-fix)' option must be either `ucr', `ucs', `ucrs' or `none'.\n";
	    cerr << "default option `-u ucrs' is used.\n\n";
	}
    }

    const char* showUnchangedStrategy()
    {
	switch(ucfix)
	{
	    case UCR:
		return "unchanged region only";
	    case UCS:
		return "unchanged subfamily only";
	    case UCRS:
	    default:
		return "unchanged region and subfamily";
	    case NoUC:
		return "disabled";
	}
    }

    void setLongNameOptions(const char* optname, const char* param)
    {
	switch(getHashValue(optname))
	{
	    case hip:
		if(strcmp(cip, optname) == 0)
		{
		    ip_name = param;
		    return;
		}
		break;
	    case hop:
		if(strcmp(cop, optname) == 0)
		{
		    op_name = param;
		    return;
		}
		break;
	    case hfmt:
		if(strcmp(cfmt, optname) == 0)
		{
		    setFormat(param);
		    return;
		}
		break;
	    case hrt:
		if(strcmp(crt, optname) == 0)
		{
		    setResidueType(param);
		    return;
		}
		break;
	    case hsm:
		if(strcmp(csm, optname) == 0)
		{
		    sm_name = param;
		    return;
		}
		break;
	    case hgp:
		if(strcmp(cgp, optname) == 0)
		{
		    gp_name = param;
		    return;
		}
		break;
	    case halg:
		if(strcmp(calg, optname) == 0)
		{
		    setGapCost(param);
		    return;
		}
		break;
	    case hdst:
		if(strcmp(cdst, optname) == 0)
		{
		    setDistanceMethod(param);
		    return;
		}
		break;
	    case hnr:
		if(strcmp(cnr, optname) == 0)
		{
		    nreconst = atoi(param);
		    return;
		}
		break;
	    case hanc:
		if(strcmp(canc, optname) == 0)
		{
		    setAnchorMethod(param);
		    return;
		}
		break;
	    case hsfb:
		if(strcmp(csfb, optname) == 0)
		{
		    sfbundl_thr = atof(param);
		    return;
		}
		break;
	    case hucf:
		if(strcmp(cucf, optname) == 0)
		{
		    setUnchangedStrategy(param);
		    return;
		}
		break;
	    case hpt:
		if(strcmp(cpt, optname) == 0)
		{
		    cerr << "\nCurrent version does not support `--"
			<< cpt << "' option; it is ignored.\n\n";
		    return;
		}
		break;
	    case hni:
		if(strcmp(cni, optname) == 0)
		{
		    outitr = atoi(param);
		    return;
		}
		break;
	    case hcut:
		if(strcmp(ccut, optname) == 0)
		{
		    cutoff = atof(param);
		    return;
		}
		break;
	    case heqf:
		if(strcmp(ceqf, optname) == 0)
		{
		    eqfactor = atof(param);
		    return;
		}
		break;
	    default:
		cerr << "\n" << optname << " is unsupported option.\n";
		exit(2);
	}
	cerr << "\noption `" << optname << "' may be typo.\n";
	exit(2);
    }
}

namespace prrn
{
    // basic options
    //
    size_t nResType = nAACode;
    rtype res_type = AA;
    string sm_name = "BLOSUM62";
    string gp_name = "gap";

    psa_alg psaa = SeqLong;
    gsa_alg msaa = ProfLong;

    dist_calc idm = Psa;

    string ip_name;
    string op_name;
    fstream out;

    opfmt opformat = FASTA;

    bool isVerbose = false;

    anchor_mtd anchoring = NoAnchor;
    // 2.0 achieves most similar results when no heuristics is employed
    // while 2.5 does slightly better results
    // if sfbundl_thr <= 0.0, subfamily bundling is disabled
    double sfbundl_thr = 2.5;
    ucf_strategy ucfix = UCRS;

    // if nreconst = 0, it does not reconstruct progressive alignment
    int nreconst = 0;


    // detailed options
    // 
    // if outitr = 0, it execute only progressive alignment phase
    int outitr = 5;
    int inritr = 8;

    SCORE cutoff = 1.0e-2;

    double eqfactor = 1.0;

    void showUsage(char** argv)
    {
	cerr << "\nUsage: " << argv[0] << " [options] [-i] input_file [-o output_file]\n";
    }

    void showHelp(char** argv)
    {
	if(argv)
	    showUsage(argv);

	cerr << "\n\nBasic options:\n\n"
	    << "  -i, --" << cip << " FILE              \tread FILE as input (FILE must be FASTA format)\n"
	    << "  -o, --" << cop << " FILE             \toutput alignment to FILE\n"
	    << "                                \tif this option is not specified, output to `stdout'\n"
	    << "  -f, --" << cfmt << "  FMT              \tuse FMT format for output\n"
	    << "                                \tFMT is `fasta' (default), `msf', `phylip', or `gde'\n"
	    << "\n"
	    << "  -t, --" << crt << " TYPE            \tinput sequences are treated as that of TYPE\n"
	    << "                                \tTYPE is `amino' (default), `dna', or `rna'\n"
	    << "  -s, --" << csm << " FILE\tuse substitution matrix FILE\n"
	    << "  -g, --" << cgp << " FILE           \tset gap cost parameters specified by FILE\n"
	    << "  -a, --" << calg << " GCOST         \tuse GCOST for group-to-group sequence alignment\n"
	    << "                                \tGCOST is `piecewise' (default) or `affine'\n"
	    << "\n"
	    << "  -d, --" << cdst << " MTHD           \tuse MTHD for initial distance calculation\n"
	    << "                                \tMTHD is `psa' (default) or `oligo'\n"
	    << "  -r, --" << cnr << " NUM   \tset no. of initial distance recalculation to NUM\n"
	    << "\n"
	    //<< "  -n, --" << canc << " MTHD             \tuse MTHD for anchor point calculation\n"
	    //<< "                                \tMTHD is `none' (calc no anchor) (default), `nongap', or `cons'\n"
	    << "  -b, --" << csfb << " VAL    \tuse VAL as a threshold of bundling subfamilies;\n"
	    << "                                \tif VAL <= 0.0, bundling is disabled\n"
	    << "                                \t(default) VAL = " << sfbundl_thr << "\n"
	    << "  -u, --" << cucf << " MTHD             \tuse MTHD for unchanged fixing\n"
	    << "                                \tMTHD is `ucr'(unchanged region) (default), `ucs'(unchanged subfamily),\n"
	    << "                                \t`ucrs'(unchanged region and subfamily), or `nouc'(disabled)\n"
	    << "\n"
	    << "  -p, --" << cpt << " ALG   \tuse ALG to calculate phylogenetic tree\n"
	    << "                                \t(currently unsupported)\n"
	    << "\n"
	    << "      --help                    \tprint this message (do not specify the other options)\n";

	cerr << "\n\nDetailed options:\n\n"
	    << "      --" << cni << " NUM \tset no. of outer iteration to NUM\n"
	    << "                                \t(default) NUM = " << outitr << "\n"
	    << "      --" << ccut << " VAL              \tset cutoff value of outer iteration to VAL\n"
	    << "                                \t(default) VAL = " << cutoff << "\n"
	    << "      --" << ceqf << " VAL \tset equalization factor to VAL\n"
	    << "                                \t(default) VAL = " << eqfactor << "\n";

	cerr << "\nIf you have any problem, please send an e-mail to <" << mail_address << ">.\n";
    }

    void showLicense()
    {
	cerr << "\n\nPRIME is free software; you can redistribute it and/or modify it\n"
	    << "under the terms of the GNU General Public License as published by\n"
	    << "the Free Software Foundation; either version 2 of the License, or\n"
	    << "(at your option) any later version.\n\n"
	    << "PRIME is distributed in the hope that it will be useful,\n"
	    << "but WITHOUT ANY WARRANTY; without even the implied warranty of\n"
	    << "MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the\n"
	    << "GNU General Public License for more details.\n\n"
	    << "You should have received a copy of the GNU General Public License\n"
	    << "along with PRIME; if not, write to the Free Software\n"
	    << "Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA\n\n";
    }

    void showArgs(int argc, const char** argv)
    {
	cerr << "\nargc:\t" << argc << "\n";
	cerr << "argv:\n";
	for(int i = 0; i < argc; ++i)
	    cerr << argv[i] << "\n";
	cerr << "\n\n";
    }

    void setParameters(int argc, char** argv)
    {
	char* p = 0;
	// skip execute command
	--argc, ++argv;
	while(argc > 0 && (p = argv[0])[0] == '-' && *++p != '\0')
	{
	    do
	    {
		if(argc == 1)
		{
		    cerr << "All options must have an argument\n";
		    cerr << "\nFor help, type `" << argv[0] << " --help'\n";
		    exit(1);
		}
		--argc, ++argv;
		switch(*p)
		{
		    case 'i':
			ip_name = argv[0];
			break;
		    case 'o':
			op_name = argv[0];
			break;
		    case 'f':
			setFormat(argv[0]);
			break;
		    case 't':
			setResidueType(argv[0]);
			break;
		    case 's':
			sm_name = argv[0];
			break;
		    case 'g':
			gp_name = argv[0];
			break;
		    case 'a':
			setGapCost(argv[0]);
			break;
		    case 'd':
			setDistanceMethod(argv[0]);
			break;
		    case 'r':
			nreconst = atoi(argv[0]);
			break;
		    case 'n':
			cerr << "\nCurrent version does not support '-n' option; it takes no effect.\n\n";
			setAnchorMethod(argv[0]);
			break;
		    case 'b':
			sfbundl_thr = atof(argv[0]);
			break;
		    case 'u':
			setUnchangedStrategy(argv[0]);
			break;
		    case 'p':
			cerr << "\nCurrent version does not support '-p' option; it is ignored.\n\n";
			break;
		    case '-':
			if(*++p != '\0')
			{
			    setLongNameOptions(p, argv[0]);
			    goto nextopts;
			}
			else
			{
			    cerr << "\n'--' is unsupported option.\n";
			    exit(2);
			}
			break;
		    default:
			cerr << "\n" << *p << " is unsupported option.\n";
			exit(2);
		}
	    }
	    while(*++p);
nextopts:
	    --argc, ++argv;
	}
	if(argc > 0 && argv != 0 && ip_name.empty())
	    ip_name = argv[0];
	if(nResType == nNACode && idm == Oligo)
	{
	    idm = Psa;
	    cerr << "Oligomer counting for nucleic acid sequences is still unsupported;\n"
		<< "Distances are calculated by PSA.\n"
		<< "Sorry for inconvenience.\n";
	}
#ifdef MSA_TEST
	isVerbose = true;
#endif
    }

    void showParameters(size_t nseq)
    {
	cerr << "\nParameters:\n"
	    << "                   input file:\t" << ip_name << " (#seq = " << nseq << ")\n"
	    << "                              \t(treated as ";
	switch(res_type)
	{
	    case DNA:
		cerr << "dna";
		break;
	    case RNA:
		cerr << "rna";
		break;
	    case AA:
	    default:
		cerr << "amino";
	}
	cerr << " sequences)\n";
	if(!op_name.empty())
	    cerr << "                  output file:\t" << op_name << "\n";
	cerr << "                output format:\t" << ((opformat == FASTA) ? "FASTA" : "MSF") << "\n"
	    << "          substitution matrix:\t" << sm_name << "\n"
	    << "           gap cost parameter:\t" << gp_name << "\n"
	    << "                     gap cost:\t" << ((psaa == SeqLong) ? "piecewise linear" : "affine") << "\n"
	    << "         distance calculation:\t" << ((idm == Psa) ? "PSA" : "Oligomer counting") << "\n"
	    << "no. of distance recalculation:\t" << nreconst << "\n"
	    //<< "     anchor point calculation:\t" << showAnchorMethod() << "\n"
	    << "           subfamily bundling:\t" << sfbundl_thr << "\n"
	    << "             unchanged fixing:\t" << showUnchangedStrategy() << "\n"
	    << "       no. of outer iteration:\t" << outitr << "\n"
	    << "       outer iteration cutoff:\t" << cutoff << "\n"
	    << "          equalization factor:\t" << eqfactor << "\n"
	    << "\n";
    }
}
