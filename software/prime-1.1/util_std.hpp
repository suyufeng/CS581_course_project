// util_std.hpp
//
// Last Modefied 8, Mar 2007
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

#ifndef __UTIL_STD_HPP__
#define __UTIL_STD_HPP__

#include <utility>
#include <functional>

struct DeleteObject
{
    template<class T>
	void operator()(const T* ptr) const
	{
	    delete ptr;
	}
};

template<typename T, typename U = T>
    struct add1st
	: public std::binary_function<std::pair<T, U>, T, std::pair<T, U> >
{
    std::pair<T, U> operator()(std::pair<T, U>& a,
	    const T& b) const
    {
	a.first += b;
	return a;
    }
};

template<typename T, typename U = T>
    struct add2nd
	: public std::binary_function<std::pair<T, U>, U, std::pair<T, U> >
{
    std::pair<T, U> operator()(std::pair<T, U>& a,
	    const U& b) const
    {
	a.second += b;
	return a;
    }
};

template<typename T, typename U = T>
    struct less_1st
	: public std::binary_function<std::pair<T, U>, T, bool>
{
    bool operator()(const std::pair<T, U>& a,
	    const T& b) const
    {
	return (a.first < b);
    }
};

template<typename T, typename U = T>
    struct less_2nd
	: public std::binary_function<std::pair<T, U>, U, bool>
{
    bool operator()(const std::pair<T, U>& a,
	    const U& b) const
    {
	return (a.second < b);
    }
};

template<typename T, typename U = T>
    struct equal_to_1st
	: public std::binary_function<std::pair<T, U>, T, bool>
{
    bool operator()(const std::pair<T, U>& a,
	    const T& b) const
    {
	return (a.first == b);
    }
};

template<typename T, typename U = T>
    struct equal_to_2nd
	: public std::binary_function<std::pair<T, U>, U, bool>
{
    bool operator()(const std::pair<T, U>& a,
	    const U& b) const
    {
	return (a.second == b);
    }
};

template<typename T>
    struct equal_to_2nd<T, double>
	: public std::binary_function<std::pair<T, double>, double, bool>
{
    bool operator()(const std::pair<T, double>& a,
	    const double& b) const
    {
	return isApproxEqual(a.second, b);
    }
};

template<typename T, typename U = T>
    struct less1st
	: public std::binary_function<std::pair<T, U>, std::pair<T, U>, bool>
{
    bool operator()(const std::pair<T, U>& a,
	    const std::pair<T, U>& b) const
    {
	return (a.first < b.first);
    }
};

template<typename T, typename U = T>
    struct less2nd
	: public std::binary_function<std::pair<T, U>, std::pair<T, U>, bool>
{
    bool operator()(const std::pair<T, U>& a,
	    const std::pair<T, U>& b) const
    {
	return (a.second < b.second);
    }
};

#endif
