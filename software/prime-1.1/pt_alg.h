// pt_alg.h
//
// Last Modified: 11, May 2007
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

#ifndef __PT_ALG_H__
#define __PT_ALG_H__

#include "matrix.hpp"

class Sequence;

namespace nj
{
    class PTAlg
    {
	private:
	    typedef std::pair<size_t, double> elm_type;
	    typedef std::vector<elm_type>::iterator iterator;
	    typedef std::vector<elm_type>::const_iterator const_iterator;
	    iterator get_iter(size_t, size_t) ;
	    const_iterator get_iter(size_t, size_t) const;

	public:
	    PTAlg(const Matrix<size_t>&);
	    ~PTAlg();

	    size_t size() const;
	    void selectOTU(double&, size_t&, size_t&, const std::vector<size_t>&, size_t) const;
	    double distance(size_t, size_t) const;
	    double length(size_t, double) const;

	    void update(size_t, size_t, const std::vector<size_t>&);

	    void show() const;

	private:
	    size_t nelm_;
	    size_t insidx_;
	    std::vector<std::vector<elm_type>*> dm_;
	    std::vector<double> totdist_;
    };
}


#endif
