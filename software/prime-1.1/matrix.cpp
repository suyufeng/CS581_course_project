// matrix.cpp
//
// Last Modified: 13, Mar 2007
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


#include "matrix.hpp"

// NOTE: If a compiler does not support a keyword `export', comment out the keyword
//export
template<>
bool isSymmetric<double>(const Matrix<double>& m)
{
    const size_t size = m.getRowSize();
    if(size != m.getColumnSize())
	return false;
    for(size_t i = 0; i < size-1; ++i)
	for(size_t j = i+1; j < size; ++j)
	    if(!isApproxEqual(m(i, j), m(j, i)))
		return false;
    return true;
}
