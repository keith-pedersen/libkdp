// Copyright (C) 2014-2018 by Keith Pedersen (Keith.David.Pedersen@gmail.com)

// Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
// The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

#ifndef KDP_TOOLS
#define KDP_TOOLS

/*! @file kdpTools.hpp
 *  @brief Some short utility functions and classes.
 * 
 *  \note As a reminder, a constexpr function \em could be evaluated at compile time
 *  (although sometimes it cannot).
 * 
 *  @author Copyright (C) 2014-2018 Keith Pedersen (Keith.David.Pedersen@gmail.com)
*/

#include <assert.h>
#include <string>
#include <vector>
#include <array>
#include <sstream>
#include <cmath>
#include <limits>
#include <mutex>
#include <exception>
#include <fstream>

// Limit memory
#include <sys/time.h>
#include <sys/resource.h>

// Testing
//~ #include <iostream>

////////////////////////////////////////////////////////////////////////

// If we're using GCC, we can tell Tell GCC to not issue warnings for
// expressions that have been reviewed and approved by the author.
#ifdef __GNUC__

/* Use case: 
 * We pass -Wfloat-equal so that GCC warns for exact floating point comparison, 
 * which is generally unsafe given floating point rounding error.
 * 
 * For example, the following is not a safe test that a vector is normalized:
 * 	
 * 	bool IsNormalized_unsafe(vec3_t const& vec)
 * 	{
 * 		return (vec.Mag2() == 1.);
 * 	}
 * 
 * because rounding error might case the square magnitude to equal 1 +/- MachineEpsilon. 
 * A safe normalization test looks like:
 * 
 * 	bool IsNormalized_safe(vec3_t const& vec)
 * 	{
 * 		return (std::fabs(vec.Mag2() - 1.) < 1e-15); // assuming vec3_t internally uses double
 * 	}
 * 
 * But, if we're writing a function that normalizes a vector,
 * we need to explicitly check for a null vector (which will not be caused by rounding error).
 * Hence, to suppress warnings, we can GCC_IGNORE_PUSH the warning we wish to ignore. 
 * To reactivate the warning, we simply call GCC_IGNORE_POP.
 * 
 * 	vec3_t Normalize(vec3_t const& vec)
 * 	{
 * 		real_t const length = vec.Mag();
 * 
 * 								GCC_IGNORE_PUSH(-Wfloat-equal)
 * 		if(length == 0.)
 * 			return vec;
 * 		else 
 * 			return vec / length;
 * 								GCC_IGNORE_POP
 *		}
 * 
 * See RelDiff() for another example.
*/ 

// This requires some pre-processor magic that is beyond my short attention span
// https://stackoverflow.com/questions/8724644/how-do-i-implement-a-macro-that-creates-a-quoted-string-for-pragma
#define HELPER0(x) #x
#define HELPER1(x) HELPER0(GCC diagnostic ignored x)
#define HELPER2(y) HELPER1(#y)

// Push back existing state, ignore warning
#define GCC_IGNORE_PUSH(x) _Pragma("GCC diagnostic push") _Pragma(HELPER2(x))
// Restore state when GCC_IGNORE_PUSH was called
#define GCC_IGNORE_POP _Pragma("GCC diagnostic pop")

#else

// If we're not using GCC, ignore these macros
#define GCC_IGNORE_PUSH(x) 
#define GCC_IGNORE_POP

#endif

namespace kdp
{
	
// Declare the @mainpage inside the namespace, so we can use unqualified names
	
	
/*! @mainpage 
 * 
 *  @brief libkdp;
 *  a collection of tools which I frequently use for many projects.
 *  
 *  @author Copyright (C) 2014-2018 Keith Pedersen (Keith.David.Pedersen@gmail.com)
 * 
 *  One of the focuses of this library is numerical stability, 
 *  as primarily exemplified by my classes for 3-vectors (with stable rotations)
 *  and 4-vectors (with stable boosts) found in kdpVectors.hpp. 
 *  Additionally, BinaryAccumulate is useful for mitigating rounding error
 *  when accumulating a large array.
 * 
 *  As of July 2018, the Python library includes all the classes in 
 *  kdpVectors.hpp, but none of the other classes or functions.
 *  
 *  \warning If compiling to a python2 library, 
 *  only "true" division is supported for vector classes.
 *  To activate true division, one can launch "python -Qnew".
*/ 

// Reminder: If there is a non-templated function or class we wish to include here, 
// we can simply declare it inline.	This will alter function lookup semantics.
// declared inline per https://stackoverflow.com/questions/3973218/header-only-libraries-and-multiple-definition-errors

////////////////////////////////////////////////////////////////////////

//! @brief Return \p val squared (better/faster than std::pow for power = 2).
template<typename real_t>
constexpr real_t Squared(real_t const val) {return val*val;}

////////////////////////////////////////////////////////////////////////

//! @brief Return a^2 - b^2 precisely (better than Squared(a) - Squared(b) when b is close to a).
template<typename real_t>
constexpr real_t Diff2(real_t const a, real_t const b)
{
	return (a - b)*(a + b);
}
 
////////////////////////////////////////////////////////////////////////

/*! @brief Biased relative error: (calculated - correct)/correct
 *  
 *  This function is safe if either input is infinity;
 *  it will only return NAN when either/both input is NAN, or both are zero.
*/ 
template<typename real_t>
real_t RelError(real_t const calculated, real_t const correct)
{
											GCC_IGNORE_PUSH(-Wfloat-equal)
	if(calculated == correct)
	{
		// Formally, lim 0/x as x->0 = 0. So even though the relative error 
		// makes no sense when correct = 0, the usable answer is zero.
		return real_t(0);
	}
											GCC_IGNORE_POP
	
	real_t relError = (calculated - correct) / correct;
	
	// We'll get nan if either argument is infinity
	if(std::isnan(relError))
	{
		// This expression will work if one argument is +/- infinity, 
		relError = calculated/correct - real_t(1);
		
		// This expression will work if both arguments are +/- infinity.
		if(std::isnan(relError))
			relError = -real_t(2)*real_t(calculated < correct);
	}

	return relError;
}

////////////////////////////////////////////////////////////////////////

//! @brief The absolute value of the biased relative error (from RelError)
template<typename real_t>
inline real_t AbsRelError(real_t const calculated, real_t const correct)
{
	return std::fabs(RelError(calculated, correct));
}

////////////////////////////////////////////////////////////////////////

/*! @brief Unbiased relative difference: (me - you)/(me + you)
 * 
 *  This version is safe if either input is infinity;
 *  it will only return nan when either/both input is nan or zero.
*/
template<typename real_t>
real_t RelDiff(real_t const me, real_t const you)
{
											GCC_IGNORE_PUSH(-Wfloat-equal)
	if(me == you)
	{
		// Formally, lim 0/(2x) as x->0 = 0. So even though the relative difference 
		// makes no sense when both are zero, the usable answer is zero.
		return real_t(0);
	}
											GCC_IGNORE_POP
	
	// This version is safe if either input is infinity;
	// it will only return nan when either/both input is nan
	const real_t alpha = me/you;
	
	// If both inputs are infinity, we'll get nan
	if(std::isnan(alpha))
	{
		// We can salvage every situation but where one/both are nan
		return real_t(1)-real_t(2)*(me < you);
	}
	else
	{
		const real_t beta = you/me;

		// When a ~= b, errors arise primarily in the numerator
		// On average, though, this form is numerically better than the original for small differences
		return (alpha - beta)/(real_t(2) + alpha + beta);
	}
}

////////////////////////////////////////////////////////////////////////

//! @brief The absolute value of the unbiased relative error (from RelDiff)
template<typename real_t>
inline real_t AbsRelDiff(real_t const me, real_t const you)
{
	return std::fabs(RelDiff(me, you));
}

////////////////////////////////////////////////////////////////////////

//! @brief Convert degrees to radians
template <typename real_t>
inline real_t ToRadians(real_t const degrees)
{
	return real_t(M_PI/180.)*degrees;
}

////////////////////////////////////////////////////////////////////////

//! @brief Convert radians to degrees
template <typename real_t>
inline real_t ToDegrees(real_t const radians)
{
	return real_t(180./M_PI)*radians;
}

////////////////////////////////////////////////////////////////////////

/*! @brief Parse a string (e.g. from an INI file) to read the angle,
 *  allowing input in both radians and degrees.
 *  
 *  Assume radians if the angle is purely numeric. Accept suffixes containing 
 *  two known substrings:
 *   - "deg" means degrees (so {degs, degree, or degrees} will also work).
 *   - "rad" means radians (so {rads, radian, or radians} will also work).
 * 
 *  \returns the parsed angle, in radians
 *  
 *  \note There must be a space between the numbers and the suffix.
*/
template<typename real_t>
real_t ReadAngle(std::string const& toParse)
{
	// static means we need only one instance of the variable for every call to the function,
	// which in this case is true and useful
	// constexpr means that the value is known at compile time,
	// so the compiler can potentially optimize the variable away.
	static constexpr char const* degreeKey = "deg"; // a std::string cannot be a constexpr
	static constexpr char const* radianKey = "rad";
		
	real_t angle;
	std::stringstream parser(toParse);
	parser >> angle;
		
	if(not parser)
		throw std::runtime_error(("ReadAngle: cannot parse <" + toParse + ">. Is there a space between the number and its units?").c_str());
	
	// Angles with no units are assumed radians; attempt to parse the units
	std::string units;
	parser >> units;
		
	if(units.length() > 0) // A suffix existed
	{
		// Look for degreeKey in the suffix (starting at position 0). 
		if(units.find(degreeKey, 0) not_eq units.npos)
		{
			//~ printf("Converted degrees to radians.\n");
			angle = kdp::ToRadians(angle);
		}
		// Else if we find radians, nothing to do (no conversion needed)
		else if(units.find(radianKey, 0) not_eq units.npos) {}
		// A non-empty string which doesn't contain 'rad' or 'deg' has an unsupported unit
		else
		{
			throw std::runtime_error(("ReadAngle: unsupported angular units ... '" + units + "'."
				 + " Only a string containing 'deg' or 'rad' work (e.g. 'deg', 'degs', 'degree' or 'degrees'"
				 + " and the same variations on 'rad')."
				 + " Radians are assumed in the absence of a qualified degree keyword,"
				  + " but we don't know what to do with this.").c_str());
		}
		// Otherwise we found radians 
	}
	// else no units to extract, return angle (in radians) "as is"
	
	return angle;
}

////////////////////////////////////////////////////////////////////////

//! @brief Partition \p n things into partitions of size \p partSize. How many do you need?
template<typename uint_t>
uint_t MinPartitions(uint_t const n, uint_t const partSize)
{
	return (n - 1)/partSize + 1;
}

////////////////////////////////////////////////////////////////////////

//! @brief Round a value to the nearest pitch (e.g. if pitch is 3, round value = 10 to 9). 
template<typename real_t>
real_t RoundToNearestPitch(real_t const value, real_t const pitch)
{
	return pitch * std::round(value / pitch);
}

////////////////////////////////////////////////////////////////////////

//! @brief sum of {1, 2, ..., n}
template<typename uint_t>
uint_t GaussSum(uint_t const n)
{
	return (n*(n + 1))/2;
}

////////////////////////////////////////////////////////////////////////

//! @brief Return the value of the smallest set bit (2**index, not index)
template<typename uint_t>
constexpr uint_t SmallestSetBit(uint_t const x)
{
	static_assert(not std::numeric_limits<uint_t>::is_signed, "SmallestBit: type must be UN-signed");
	
	// Developed from 
	// http://www.exploringbinary.com/ten-ways-to-check-if-an-integer-is-a-power-of-two-in-c/
	return (x bitand (compl(x) + 1));
}

////////////////////////////////////////////////////////////////////////

/*! @brief Return the value of the largest set bit (2**index, not index)
 * 
 *  Does a stupid search (not x86 BSR instruction), 
 *  b/c x86 tutorials are difficult to understand (and not portable to other arch).
 *  A binary search *might* be faster, but uses much more code.
*/ 
template<typename uint_t>
uint_t LargestBit(uint_t x)
{
	static_assert(not std::numeric_limits<uint_t>::is_signed, "LargestBit: type must be UN-signed");
	
	if(x == 0) return 0;
	
	// We can start with the smallest set bit and move on from there
	uint_t largestBit = SmallestSetBit(x);
	x /= (2 * largestBit);
	
	while(x)
	{
		x /= 2;
		largestBit *= 2;
	}
	return largestBit;
}

////////////////////////////////////////////////////////////////////////

/*! @brief Determines if the argument is an exact power of 2
*/ 
template<typename uint_t>
constexpr bool IsPowerOfTwo(uint_t const x)
{
	static_assert(not std::numeric_limits<uint_t>::is_signed, "IsPowerOfTwo: type must be UN-signed");
	
	return ((x not_eq 0) and (x == SmallestSetBit(x)));
}

////////////////////////////////////////////////////////////////////////

/*! @brief Binary accumulate a std::array with \p size = (exact power of 2).
 * 
 *  @warning The array (passed by ref) is accumulated in place (i.e. destructively altered).
*/
template<typename T, size_t arraySize>
T BinaryAccumulate_Destructive(std::array<T, arraySize>& vec, size_t size = arraySize)
{
	// We can use a static assert because std::array has static size
	static_assert(IsPowerOfTwo(arraySize), "BinaryAccumulate(std::array): array size must be an exact power of 2");
		
	if(not kdp::IsPowerOfTwo(size))
		throw std::runtime_error("BinaryAccumulate: starting size must be a power of two!");
	
	while(size > 1)
	{
		size /= 2; // Equivalent to >>=, but more readible
		for(size_t i = 0; i < size; ++i)
			vec[i] += vec[size + i];
	}
	
	return vec.front();
}

////////////////////////////////////////////////////////////////////////

/*! @brief Binary accumulate a std::array with size = (exact power of 2).
 * 
 *  The accumulation is non-destructive (array is passed by value, 
 *  then passed to BinaryAccumulate_Destructive()).
*/
template<typename T, size_t arraySize>
T BinaryAccumulate(std::array<T, arraySize> vec, size_t const size = arraySize)
{
	return BinaryAccumulate_Destructive(vec, size);
}

////////////////////////////////////////////////////////////////////////

/*! @brief Binary accumulate a std::vector of any size.
 * 
 *  @warning The array (passed by ref) is accumulated in place (i.e. destructively altered).
 *  For efficiency, the first reduction does not span 100% of elements, 
 *  but makes the vector's new size an exact power of two. This speeds up all 
 *  subsequent reductions, which do not require an oddness check
 *  and will have good alignment in memory.
*/  
template<typename T>
T BinaryAccumulate_Destructive(std::vector<T>& vec)
{
	if(vec.size())
	{
		size_t size = vec.size();
		
		if(size > 1) // If size == 1, skip the reduction
		{
			// 1. Find the smallest power-of-two less than vec.size()
			size = LargestBit(vec.size());
			
			{
				// 2. Determine how many values are past the power of two
				size_t const overflow = vec.size() - size;
				assert(overflow < size); // unsigned arithmetic, so also asserts overflow >= 0
			
				// 3. Add the overflow to the front of the list (assert guarantees enough room)
				for(size_t i = 0; i < overflow; ++i)
					vec[i] += vec[size + i];
			}
			
			// Now we can do a fast, power-of-two accumulate (no oddness check)
			vec.resize(size); // <== Is this helpfull?
			
			while(size > 1)
			{
				size /= 2; // Equivalent to >>=, but more readible
				for(size_t i = 0; i < size; ++i)
					vec[i] += vec[size + i];
			}
		}
		
		return vec.front();
	}
	return T(); // Is this correct when nothing is summed? Should we instead throw an exception?
}

////////////////////////////////////////////////////////////////////////

/*! @brief Binary accumulate a std::vector of any size.
 * 
 *  The accumulation is non-destructive (array is passed by value, 
 *  then passed to BinaryAccumulate_Destructive()).
*/
template<typename T>
T BinaryAccumulate(std::vector<T> vec)
{
	return BinaryAccumulate_Destructive(vec);
}

////////////////////////////////////////////////////////////////////////

/*! @brief Use a linux system call to set a soft limit on this processes' memory use.
 * 
 *  \throws throws \c std::bad::alloc when the limit is reached
 *  
 *  \note Only the total memory use can be constrained, including the instruction stack, 
 *  so expect at least 1 MB to be unavailable for use as data.
 * 
 *  This function was developed after numerous cases where a buggy program
 *  (or a perfectly fine one which was given an insane task)
 *  chewed through the physical memory and started using swap.
 *  This totally hangs my system, even if the program has nice +20, 
 *  because the CPU is spending all it's time swapping, 
 *  and swapd is nice'd at -20. 
*/ 
inline void LimitMemory(double const num_GB = 4.)
{
	// Declare inline so it can be defined in a header only without leading to multiple copies
	
	static_assert(sizeof(long) >= 8, "Must have 64-bit addressing");
	
	struct rlimit current;
	int const resource = RLIMIT_AS; // Virtual address space. 
	// I tried RLIMIT_DATA, but it offered no protection for std::vector
	
	int ret = getrlimit(resource, &current);
	if(ret)
		throw std::runtime_error("LimitMemory: Cannot get current limits.");
		
	current.rlim_cur = long(num_GB * double(long(1) << 30));
	if(current.rlim_cur >= current.rlim_max)
		throw std::runtime_error("LimitMemory: Supplied limit too large!");
		
	ret = setrlimit(resource, &current);
	if(ret)
		throw std::runtime_error("LimitMemory: Cannot set new soft limit.");
}

////////////////////////////////////////////////////////////////////////

/*! @brief Use Welford's algorithm to calculate a running mean and variance
 * 
 *  This provides a numerically stable estimate of the variance, 
 *  but without having to store the entire vector.
 * 
 *  \note There is no function for the standard deviation because some T may not
 *  fit into std::sqrt (e.g. a std::vector<double>).
*/ 
template <typename T>
class WelfordEstimate
{
	private:
		T mean; 
		T sumOfSquareDeviations;
		
		size_t n;
		
	public:
		WelfordEstimate():
			mean(), sumOfSquareDeviations(), n(0) {}
		WelfordEstimate(WelfordEstimate const&) = default;
		WelfordEstimate& operator=(WelfordEstimate const&) = default;
		
		// Enable move semantics, in case T is a vector
		WelfordEstimate(WelfordEstimate&&) = default;
		WelfordEstimate& operator=(WelfordEstimate&&) = default;
		
		WelfordEstimate& operator+=(T const& val)
		{
			T deltaOld = (val - mean);
			
			// Use double for denominator (even if T = float) because n could be very large
			mean += deltaOld/double(++n);
			sumOfSquareDeviations += deltaOld*(val - mean);
					
			if(n == 0lu)
				throw std::runtime_error("WelfordEstimate: sample size exceeded capacity of size_t; I'm impressed.");
				
			return *this;
		}
		
		T const& Mean() const {return mean;}
		T Variance() const {return sumOfSquareDeviations/double(n);}
		T Variance_Unbiased() const {return sumOfSquareDeviations/double(n - size_t(1));}
		
		size_t Size() {return n;}
};

template <typename T>
T PythagoreanSum(std::vector<T> vec)
{
	T sum(0); // zero-initializing constructor
	
	for(auto const& val : vec)
		sum += kdp::Squared(val);
		
	return std::sqrt(sum);	
}

template <typename T>
T PythagoreanSum(std::initializer_list<T> vec)
{
	return PythagoreanSum(std::vector<T>(vec));
}

/*! @brief Test that a file can be opened and read.
 * 
 *  Using ifstream is most portable across OS, but not the fastest.
 *  However, as long as this is not used often, its speed is irrelevant.
*/ 
// Declare inline to allow compiler to consolidate multiple copies of function in 
// various object files gathered into a single program/library.
inline bool FileIsReadable(std::string const& filePath)
{
	return std::ifstream(filePath).is_open();
}

/*! @brief Quickly check if a file is visible.
 * 
 *  Mainly useful for control logic before sending file to another function for internal use.
*/
// Based on: http://stackoverflow.com/questions/12774207/fastest-way-to-check-if-a-file-exist-using-standard-c-c11-c
//~ inline bool FileIsVisible(const std::string& name) 
//~ {
	//~ // per: http://stackoverflow.com/questions/23329382/function-and-struct-having-the-same-name-in-c,
	//~ // A struct and a function/variable/etc can have the same name. 
	//~ // The function/variable/etc hides the struct, 
	//~ // which is identified by the struct prefix
	//~ struct stat buffer; // POSIX uses this name hiding paradigm to associate structs with functions
	//~ // stat::Return(): Upon successful completion, 0 shall be returned. 
	//~ // Otherwise, -1 shall be returned and errno set to indicate the error.
	//~ return (stat(name.c_str(), &buffer) == 0); // Any error means we can't read the file.
//~ }

}// end namespace

#endif
