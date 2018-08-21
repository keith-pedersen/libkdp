#ifndef KDP_SETTINGS
#define KDP_SETTINGS

// Copyright (C) 2018 by Keith Pedersen (Keith.David.Pedersen@gmail.com)
//
// Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
// The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

#include <QtCore/QSettings>

namespace kdp
{
	
//! @brief A Settings class providing a generic parameter class
struct Settings_Base
{
	/*! @brief A parameter class for use with QSettings
	 *  
	 *  Read() will search a QSettings object for the parameter, 
	 *  specified by its \ref key, and use the \ref value_default if
	 *  they key is not found or no value is specified.
	 *  The purpose of this class is to define the key and default value in one line, 
	 *  for readability
	 * 
	 *  	struct Settings : Settings_Base
	 * 	{
	 *  		Param<double> setting1 = Param<double>("setting1", 3.);
	 *  		Param<size_t> setting2 = Param<size_t>("setting2", 14);
	 * 	};
	 * 
	 *  Ideally \ref key and \ref value_default would be static members of the class, 
	 *  or even better template parameters. Unfortunately, C++ doesn't support string literals
	 *  as template parameters, and defining static members defeats 
	 *  the whole purpose of this class.
	 * 
	 *  A hack to emulate static members is to make \ref key and \ref value_default constant,
	 *  which prevents them from being changed outside the derived Settings class declaration.
	 *  The only side effect to this hack is that it deletes operator= for Settings, 
	 *  which will delete operator= for any class containing Settings.
	 *  To fix this, we overload operator= to only change the value.
	 * 
	 *  Implicit conversion to T is allowed through operators, 
	 *  but some functions may not be able to deduce the correct type, 
	 *  and one may be forced to explicitly access \ref value.
	 */ 
	template<typename T>
	struct Param
	{
		std::string const key;
		std::string const value_default;
		T value;
		
		//! @brief Define the key and default value; default-initialize the value
		Param(std::string const& key_in, std::string const*& value_default_in):
			key(key_in), value_default(value_default_in), value() {}
		
		//! @brief Define the key and default value; default-initialize the value	
		Param(std::string const& key_in, char const* const value_default_in):
			key(key_in), value_default(value_default_in), value() {}		
		
		//! @brief Define the key and default value; default-initialize the value
		template<typename U>
		Param(std::string const& key_in, U const& value_default_in):
			key(key_in), value_default(std::to_string(value_default_in)), value() {}
		
		operator T&() {return value;}
		operator T const& () const {return value;}
			
		//! @brief Set the value of the parameter
		Param& operator=(T const& val)
		{
			value = val; 
			return *this;
		}
		
		//! @brief Set one parameter equal to another (transferring only its value).
		Param& operator=(Param const& that)
		{
			this->value = that.value; 
			return *this;
		}
		
		/*! @brief Read the value from the INI file (looking for \ref key inside the given \p section).
		 * 
		 *  For many types, we can rely on template function QVariant::value<T>.
		 *  However, I provide special instances for the following types:
		 *  
		 *   - \c int, int64_t, uint64_t: First convert to double, then round to size_t.
		 *     This allows the specification of \em large integers via
		 *     floating point specification (e.g., 1e4).
		 *     Additionally, Qt will natively round 13.1 to 0 for integer types
		 *     \warning Any integers larger than (2^53-1) will lose precision
		 *     during the conversion to double, and any larger than (2^64-1)
		 *     will round to zero.
		 * 
		 *   - std::string: Qt does not define QVariant::value<std::string>, 
		 *     only QVariant::value<QString>.
		 * 
		 *  \note Using std::stringstream for the conversion solves the 
		 *  rouding 13.1 -> 0 for integers problem, but does not provide
		 *  natural conversion of "true" -> true for bool, 
		 *  and also interprets "1e3" -> 1 for integer types, 
		 *  so it's an even worse solution.
		 * 
		 *  \param section
		 *  The section of the INI file where the key should be found, 
		 *  creating the full key "section/key"
		*/  
		T& Read(QSettings const& parsedINI, std::string const& section = "");
		
		/*! @brief Read the QVariant without setting the value
		 * 
		 *  \param section
		 *  The section of the INI file where the key should be found, 
		 *  creating the full key "section/key"
		*/  
		QVariant ReadVariant(QSettings const& parsedINI, std::string const& section) const
		{
			 return parsedINI.value((section + "/" + key).c_str(), value_default.c_str());
		}
	};
};

// These must be declared inline so we don't get multiple definitions

template<typename T>
inline T& Settings_Base::Param<T>::Read(QSettings const& parsedINI, std::string const& section)
{
	// I do not understand why I need to add this static_cast; but otherwise the compiler can't interpret
	return (value = static_cast<QVariant>(ReadVariant(parsedINI, section)).value<T>());
}
			
template<>
inline std::string& Settings_Base::Param<std::string>::Read(QSettings const& parsedINI, 
	std::string const& section)
{
	return (value = ReadVariant(parsedINI, section).toString().toStdString());
}

template<>
inline uint64_t& Settings_Base::Param<uint64_t>::Read(QSettings const& parsedINI, 
	std::string const& section)
{
	return (value = size_t(ReadVariant(parsedINI, section).toDouble()));
}

template<>
inline int& Settings_Base::Param<int>::Read(QSettings const& parsedINI, 
	std::string const& section)
{
	return (value = int(ReadVariant(parsedINI, section).toDouble()));
}

template<>
inline int64_t& Settings_Base::Param<int64_t>::Read(QSettings const& parsedINI, 
	std::string const& section)
{
	return (value = int64_t(ReadVariant(parsedINI, section).toDouble()));
}

// End namespace
}

#endif
