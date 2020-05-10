// matrix.hpp
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

#ifndef __MATRIX_H__
#define __MATRIX_H__

#include "util.h"
#include <iosfwd>
#include <new>
#include <vector>
#include <algorithm>
#include <functional>
#ifdef DEBUG
#include <cassert>
#endif

template<typename T>
class Matrix
{
public:
	Matrix(){}

	Matrix(size_t row, size_t col)
	    : elm_(std::vector<std::vector<T> >(row, std::vector<T>(col)))
	{}

	Matrix& operator+=(const Matrix&);
	Matrix& operator*=(const T&);
	~Matrix(){}

	T& operator()(size_t, size_t);
	const T& operator()(size_t, size_t) const;

	void resize(size_t, size_t);

	void set(size_t, size_t, const T&);
	T get(size_t, size_t) const;

	size_t getRowSize() const;
	size_t getColumnSize() const;
private:
	std::vector<std::vector<T> > elm_;
};

template<typename T>
Matrix<T> &Matrix<T>::operator+=(const Matrix& sm)
{
    const size_t row = elm_.size();
#ifdef DEBUG
    assert(row == sm.elm_.size());
#endif
    const size_t col = elm_.begin()->size();
#ifdef DEBUG
    assert(col == sm.elm_.begin()->size());
#endif

    for(size_t i = 0; i < row; ++i)
	for(size_t j = 0; j < col; ++j)
	    elm_[i][j] += sm.elm_[i][j];

    return *this;
}

template<typename T>
Matrix<T> &Matrix<T>::operator*=(const T& val)
{
    const size_t row = elm_.size();
    const size_t col = elm_.begin()->size();

    for(size_t i = 0; i < row; ++i)
	for(size_t j = 0; j < col; ++j)
	    elm_[i][j] *= val;

    return *this;
}

template<typename T>
T& Matrix<T>::operator()(size_t i, size_t j)
{
#ifdef DEBUG
    return elm_.at(i).at(j);
#else
    return elm_[i][j];
#endif
}

template<typename T>
const T& Matrix<T>::operator()(size_t i, size_t j) const
{
#ifdef DEBUG
    return elm_.at(i).at(j);
#else
    return elm_[i][j];
#endif
}

template<typename T>
void Matrix<T>::resize(size_t row, size_t col)
{
    elm_.resize(row);
    std::for_each(elm_.begin(), elm_.end(),
	    std::bind2nd(std::mem_fun_ref(&std::vector<T>::resize), col));
}

template<typename T>
void Matrix<T>::set(size_t row, size_t col, const T& val)
{
#ifdef DEBUG
    elm_.at(row).at(col) = val;
#else
    elm_[row][col] = val;
#endif
}

template<typename T>
inline T Matrix<T>::get(size_t row, size_t col) const
{
#ifdef DEBUG
    return elm_.at(row).at(col);
#else
    return elm_[row][col];
#endif
}

template<typename T>
size_t Matrix<T>::getRowSize() const
{
    return elm_.size();
}

template<typename T>
size_t Matrix<T>::getColumnSize() const
{
#ifdef DEBUG
    assert(!elm_.empty());
#endif
    return elm_.begin()->size();
}

template<typename T>
bool isSymmetric(const Matrix<T>& m)
{
    const size_t size = m.getRowSize();
    if(size != m.getColumnSize())
	return false;
    for(size_t i = 0; i < size-1; ++i)
	for(size_t j = i+1; j < size; ++j)
	    if(m(i, j) != m(j, i))
		return false;
    return true;
}

template<>
bool isSymmetric<double>(const Matrix<double>&);

#endif
