// Copyright (C) 2016-2018 by Keith Pedersen (Keith.David.Pedersen@gmail.com)

// Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
// The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

#include <vector>

#ifndef KDP_STD_VECTOR_MATH
#define KDP_STD_VECTOR_MATH

/*! @file kdpStdVectorMath.hpp
 *  @brief Defines operators that allow std::vector to be added, subtracted, 
 *  and multiplied (very much like numpy or FORTRAN arrays).
 * 
 *  They are declared outside the kdp namespace to make them easier to call.
 * 
 *  @author Copyright (C) 2018 Keith Pedersen (Keith.David.Pedersen@gmail.com)
 *  @date 2018
*/

template<typename T>
std::vector<T>& operator+=(std::vector<T>& lhs, std::vector<T> const& rhs)
{
	if(rhs.size() not_eq lhs.size())
		throw std::range_error("operator+= (std::vector<T>): vectors have different lengths");
	
	for(size_t i = 0; i < lhs.size(); ++i)
		lhs[i] += rhs[i];
	
	return lhs;
}

template<typename T>
std::vector<T> operator+(std::vector<T> const& lhs, std::vector<T> const& rhs)
{
	std::vector<T> copy(lhs);
	return copy += rhs;
}

template<typename T>
std::vector<T>& operator-=(std::vector<T>& lhs, std::vector<T> const& rhs)
{
	if(rhs.size() not_eq lhs.size())
		throw std::range_error("kdp::operator-= (std::vector<T>): vectors have different lengths");
	
	for(size_t i = 0; i < lhs.size(); ++i)
		lhs[i] -= rhs[i];
	
	return lhs;
}

template<typename T>
std::vector<T> operator-(std::vector<T> const& lhs, std::vector<T> const& rhs)
{
	std::vector<T> copy(lhs);
	return copy -= rhs;
}

template<typename T>
std::vector<T>& operator*=(std::vector<T>& lhs, T const& rhs)
{
	for(size_t i = 0; i < lhs.size(); ++i)
		lhs[i] *= rhs;
	
	return lhs;
}

template<typename T>
std::vector<T> operator*(std::vector<T> const& lhs, T const& rhs)
{
	std::vector<T> copy(lhs);
	return copy *= rhs;
}

template<typename T>
std::vector<T>& operator*=(std::vector<T>& lhs, std::vector<T> const& rhs)
{
	if(rhs.size() not_eq lhs.size())
		throw std::range_error("kdp::operator*= (std::vector<T>): vectors have different lengths");
	
	for(size_t i = 0; i < lhs.size(); ++i)
		lhs[i] *= rhs[i];
	
	return lhs;
}

template<typename T>
std::vector<T> operator*(std::vector<T> const& lhs, std::vector<T> const& rhs)
{
	std::vector<T> copy(lhs);
	return copy *= rhs;
}

#endif
