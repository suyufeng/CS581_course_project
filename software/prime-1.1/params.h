// params.h
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

#ifndef __PARAMS_H__
#define __PARAMS_H__

#include "prrn.h"
#include "residue.h"
#include "alignment.h"
#include "anchor.h"
#include <string>
#include <fstream>

namespace prrn
{
    // administrative information
    //
    const std::string author = "Shinsuke Yamada";
    const std::string mail_address = "prime@yama.info.waseda.ac.jp";
    const std::string version = "1.1";

    // basic options
    //
    extern size_t nResType;
    extern rtype res_type;
    extern std::string sm_name;
    extern std::string gp_name;

    extern psa_alg psaa;
    extern gsa_alg msaa;

    extern dist_calc idm;

    extern std::string ip_name;
    extern std::string op_name;
    extern std::fstream out;

    extern opfmt opformat;

    extern bool isVerbose;

    extern anchor_mtd anchoring;
    // if sfbundl_thr <= 0.0, subfamily bundling is disabled
    extern double sfbundl_thr;
    extern ucf_strategy ucfix;


    // if nreconst = 0, it does not reconstruct progressive alignment
    extern int nreconst;


    // detailed options
    // 
    // if outitr = 0, it execute only progressive alignment phase
    extern int outitr;
    extern int inritr;

    extern SCORE cutoff;

    extern double eqfactor;

    void showUsage(char**);
    void showLicense();
    void showHelp(char** = 0);
    void showArgs(int, const char**);

    void setParameters(int, char**);

    void showParameters(size_t);
}

#endif
