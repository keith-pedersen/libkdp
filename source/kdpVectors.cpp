// Copyright (C) 2014-2018 by Keith Pedersen (Keith.David.Pedersen@gmail.com)

// Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
// The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

#include "kdpVectors.hpp"
#include <limits>
#include <assert.h>
#include <array>
#include <stdexcept>

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	 __     __        _            ____  
	 \ \   / /__  ___| |_ ___  _ _|___ \ 
	  \ \ / / _ \/ __| __/ _ \| '__|__) |
		\ V /  __/ (__| || (_) | |  / __/ 
		 \_/ \___|\___|\__\___/|_| |_____|
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

template<typename real_t>  //+++
kdp::Vector2<real_t>::Vector2(): x1(real_t(0)), x2(real_t(0)) {}

template<typename real_t>  //+++
kdp::Vector2<real_t>::Vector2(bool) {} // Don't initilize, you know what you're doing

template<typename real_t>  //+++
kdp::Vector2<real_t>::Vector2(real_t const x1_init, real_t const x2_init):
	x1(x1_init), x2(x2_init) {}
  
template<typename real_t>  //+++
kdp::Vector2<real_t>::Vector2(real_t const w1, real_t const w2, Vec2from const argFormat):
	Vector2<real_t>(false)
{
	switch(argFormat)
	{
		// w1 = length
		// w2 = phi	
		case Vec2from::LengthPhi:
			x1 = w1 * std::cos(w2);
			x2 = w1 * std::sin(w2);
		break;
	}
}
	
template<typename real_t>  //+++
kdp::Vector2<real_t>& kdp::Vector2<real_t>::operator+=(Vector2 const& that)
{
	this->x1 += that.x1;
	this->x2 += that.x2;
	return *this;
}

template<typename real_t>  //+++
kdp::Vector2<real_t>& kdp::Vector2<real_t>::operator-=(Vector2 const& that)
{
	this->x1 -= that.x1;
	this->x2 -= that.x2;
	return *this;
}

template<typename real_t>  //+++
kdp::Vector2<real_t>& kdp::Vector2<real_t>::operator*=(real_t const scale)
{
	this->x1 *= scale;
	this->x2 *= scale;
	return *this;
}

template<typename real_t>  //+++
kdp::Vector2<real_t>& kdp::Vector2<real_t>::operator/=(real_t const scale)
{
	this->x1 /= scale;
	this->x2 /= scale;
	return *this;
}

template<typename real_t>  //+++
kdp::Vector2<real_t> kdp::Vector2<real_t>::operator-() const
{
	return Vector2(-x1, -x2);
}

template<typename real_t>  //+++
kdp::Vector2<real_t>& kdp::Vector2<real_t>::Normalize()
{
	real_t const length = Mag();
	if(length > real_t(0))
		(*this) /= length;
	return *this;
}

template<typename real_t>  //+++
real_t kdp::Vector2<real_t>::Mag() const
{
	return std::sqrt(Mag2());
}

template<typename real_t>  //++
real_t kdp::Vector2<real_t>::Mag2() const
{
	return std::fma(x2, x2, x1*x1);
}

template<typename real_t> //+
real_t kdp::Vector2<real_t>::Phi() const
{
	return std::atan2(x2, x1);
}

template<typename real_t> //+
real_t kdp::Vector2<real_t>::DeltaPhi(Vector2<real_t> const& that) const
{
	return this->Vector2::InteriorAngle(that);
}

template<typename real_t>  //+++
real_t kdp::Vector2<real_t>::Dot(Vector2 const& that) const
{
	return std::fma(this->x2, that.x2, this->x1 * that.x1);
}

template<typename real_t>  //+++
real_t kdp::Vector2<real_t>::Cross(Vector2 const& that) const
{
	return std::fma(this->x1, that.x2, -this->x2*that.x1);
}

template<typename real_t>  //+++
real_t kdp::Vector2<real_t>::InteriorAngle(Vector2 const& that) const
{
	// This version has better precision than acos of the normalized dot product.
	// Plus, in rare cases, the normalized dot product can actually have a
	// absolute magnitude greater than one, causing acos() to return non-a-number.  
	return std::atan2(std::fabs(this->Cross(that)), this->Dot(that));
}

template<typename real_t>
template<typename convertTo>  //+++
kdp::Vector2<real_t>::operator kdp::Vector2<convertTo>() const
{
	return kdp::Vector2<convertTo>(convertTo(x1), convertTo(x2));
}

// Instantiate common types
template class kdp::Vector2<double>;
template class kdp::Vector2<float>;
//template class kdp::Vector2<long double>;

// Instantiate conversion functions between common types
// (nested template, not instantiated by struct instantiation. 
template kdp::Vector2<float>::operator kdp::Vector2<double>() const;
template kdp::Vector2<double>::operator kdp::Vector2<float>() const;


/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	 __     __        _            _____ 
	 \ \   / /__  ___| |_ ___  _ _|___ / 
	  \ \ / / _ \/ __| __/ _ \| '__||_ \ 
		\ V /  __/ (__| || (_) | |  ___) |
		 \_/ \___|\___|\__\___/|_| |____/ 
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

template<typename real_t> //+++
kdp::Vector3<real_t>::Vector3(): Vector2<real_t>(), x3(real_t(0)) {}

template<typename real_t>  //+++
kdp::Vector3<real_t>::Vector3(bool): Vector2<real_t>(false) {} // Don't initilize, you know what you're doing

template<typename real_t> //+++
kdp::Vector3<real_t>::Vector3(real_t const x1_init, real_t const x2_init, real_t const x3_init):
	Vector2<real_t>(x1_init, x2_init), x3(x3_init) {}
	
template<typename real_t> //+++
kdp::Vector3<real_t>::Vector3(real_t const w1, real_t const w2, real_t const w3, 
	Vec3from const argFormat):
Vector3<real_t>(false)
{
	switch(argFormat)
	{
		case Vec3from::LengthEtaPhi:
		{
			// w1 -> length
			// w2 -> pseudorapidity eta
			// w3 -> azimuthal angle phi
			real_t const coshEta = std::cosh(w2);
			
			x1 = (w1 * std::cos(w3))/coshEta;
			x2 = (w1 * std::sin(w3))/coshEta;
			x3 = w1 * std::tanh(w2);
		}
		break;
		
		case kdp::Vec3from::LengthThetaPhi:
		{
			// w1 -> length
			// w2 -> polar angle theta
			// w3 -> azimuthal angle phi
			real_t const sinTheta = std::sin(w2);
			
			x1 = w1 * sinTheta * std::cos(w3);
			x2 = w1 * sinTheta * std::sin(w3);
			x3 = w1 * std::cos(w2);
		}
		break;
	}
}

template<typename real_t> //+++
kdp::Vector3<real_t>& kdp::Vector3<real_t>::operator+=(Vector3 const& that)
{
	this->Vector2<real_t>::operator += (that);
	this->x3 += that.x3;
	return *this;
}

template<typename real_t> //+++
kdp::Vector3<real_t>& kdp::Vector3<real_t>::operator-=(Vector3 const& that)
{
	this->Vector2<real_t>::operator -= (that);
	this->x3 -= that.x3;
	return *this;
}

template<typename real_t> //+++
kdp::Vector3<real_t>& kdp::Vector3<real_t>::operator*=(real_t const scale)
{
	this->Vector2<real_t>::operator *= (scale);
	this->x3 *= scale;
	return *this;
}

template<typename real_t> //+++
kdp::Vector3<real_t>& kdp::Vector3<real_t>::operator/=(real_t const scale)
{
	this->Vector2<real_t>::operator /= (scale);
	this->x3 /= scale;
	return *this;
}

template<typename real_t> //++
kdp::Vector3<real_t> kdp::Vector3<real_t>::operator-() const
{
	return Vector3(-x1, -x2, -x3);
}

template<typename real_t>  //++
kdp::Vector3<real_t>& kdp::Vector3<real_t>::Normalize()
{
	real_t const length = this->Mag();
	if(length > real_t(0))
		(*this) /= length;
	return *this;
}

template<typename real_t> //++
real_t kdp::Vector3<real_t>::Mag() const
{
	return std::sqrt(Mag2());
}

template<typename real_t> //++
real_t kdp::Vector3<real_t>::Mag2() const
{
	return std::fma(x3, x3, Vector2<real_t>::Mag2());
}

template<typename real_t> //++
real_t kdp::Vector3<real_t>::Eta() const
{
	/* pseudorapidity = 0.5*log((|x| + x3)/(|x| - x3)) = 0.5*log1p(2*x3/(|x|-x3))
	 *
	 * Unlike rapidity, the largeness of the argument to log1p is
	 * caused by a CATASTROPHIC cancellation. Thus, we will do
	 * some floating point gymnastics and branch on sign(x3)
	 * (since sign(eta) = sign(x3)).*/

	// 2*x3/(sqrt(x^2) - x3) = 2*x3*(sqrt(x^2) + x3)/(x^2 - x3^2)
	//                       = 2*x3*(sqrt(x^2) + x3)/(x1^2 + x2^2)
	const real_t xT2 = T().Mag2();
	const real_t absZ = std::fabs(x3);
	return std::copysign(real_t(0.5)*std::log1p(
		(real_t(2)*(absZ*(std::sqrt(std::fma(absZ, absZ, xT2)) + absZ)))/xT2), x3);
}

template<typename real_t> //++
real_t kdp::Vector3<real_t>::DeltaEta(Vector3 const& that) const
{
	return this->Eta() - that.Eta();
}

template<typename real_t> //++
real_t kdp::Vector3<real_t>::Theta() const
{
	return std::atan2(T().Mag(), x3);
}

template<typename real_t> //++
real_t kdp::Vector3<real_t>::Dot(Vector3 const& that) const
{
	return std::fma(this->x3, that.x3, this->Vector2<real_t>::Dot(that));
}

template<typename real_t> //+++
kdp::Vector3<real_t> kdp::Vector3<real_t>::Cross(Vector3 const& that) const
{
	return Vector3(
	 std::fma(this->x2, that.x3, -this->x3*that.x2),
	 std::fma(this->x3, that.x1, -this->x1*that.x3),
	 this->Vector2<real_t>::Cross(that));
}

template<typename real_t> //++
real_t kdp::Vector3<real_t>::InteriorAngle(Vector3 const& that) const
{
	// This version has better precision than acos of the normalized dot product.
	// Plus, in rare cases, the the normalized dot product can actually have a
	// absolute magnitude greater than one, leading to nan.  
	return std::atan2(this->Cross(that).Mag(), this->Dot(that));
}

template<typename real_t>
template<typename convertTo>  //++
kdp::Vector3<real_t>::operator kdp::Vector3<convertTo>() const
{
	return kdp::Vector3<convertTo>(convertTo(x1), convertTo(x2), convertTo(x3));
}

// Insantiate common types
template class kdp::Vector3<double>;
template class kdp::Vector3<float>;
//template class kdp::Vector2<long double>;

// Insantiate conversion functions between common types
// (nested template, not instantiated by struct instantiation 
template kdp::Vector3<float>::operator kdp::Vector3<double>() const;
template kdp::Vector3<double>::operator kdp::Vector3<float>() const;

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	 __     __        _             _  _   
	 \ \   / /__  ___| |_ ___  _ __| || |  
	  \ \ / / _ \/ __| __/ _ \| '__| || |_ 
		\ V /  __/ (__| || (_) | |  |__   _|
		 \_/ \___|\___|\__\___/|_|     |_|  
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/  

template<typename real_t>  //++
kdp::Vector4<real_t>::Vector4(): Vector3<real_t>(), x0(real_t(0)) {}

template<typename real_t>  //++
kdp::Vector4<real_t>::Vector4(bool): Vector3<real_t>(false) {} // Don't initilize, you know what you're doing

template<typename real_t>  //++
kdp::Vector4<real_t>::Vector4(real_t const x0_init,
real_t const x1_init, real_t const x2_init, real_t const x3_init):
Vector3<real_t>(x1_init, x2_init, x3_init), x0(x0_init) {}

template<typename real_t>  //+
kdp::Vector4<real_t>::Vector4(real_t const w0, real_t const w1, real_t const w2, real_t const w3,
	Vec4from4 const argFormat): Vector3<real_t>(false)
{
								GCC_IGNORE_PUSH(-Wfloat-equal)
	switch(argFormat)
	{
		case kdp::Vec4from4::EnergyEtaPhiMass:
		{
			{
				// w0 -> E,
				// w1 -> eta
				// w2 -> phi
				// w3 -> mass
				real_t const momentum = (w3 == real_t(0)) ? w0 : std::sqrt(kdp::Diff2(w0, w3));
				x() = Vector3<real_t>(momentum, w1, w2, Vec3from::LengthEtaPhi);
			}
			x0 = w0;
		}
		break;
	}
								GCC_IGNORE_POP
}

template<typename real_t>
kdp::Vector4<real_t>::Vector4(Vector3<real_t> const& x_in): 
	Vector3<real_t>(x_in), x0(this->x().Mag()) {}

template<typename real_t>
kdp::Vector4<real_t>::Vector4(real_t const w0, Vector3<real_t> const& x_in, 
	Vec4from2 const w0type):
	Vector3<real_t>(x_in)
{
	switch(w0type)
	{
		case Vec4from2::Energy:
		case Vec4from2::Time:
			x0 = w0;
		break;
		
		case Vec4from2::Length:
		case Vec4from2::Mass:
			if(w0 < real_t(0))
				throw std::domain_error("Vector4: negative mass");
				
			x0 = std::sqrt(std::fma(w0, w0, x().Mag2()));
		break;	
		
		default:
			throw std::runtime_error("Vector4::Vecto4. Unrecognized Vec4from2 enum");
		break;
	}
}

//~ case Vec4from2::Boost_preserve_p3:
	//~ if(w0 < real_t(1))
		//~ throw std::domain_error("Vector4: boost factor cannot be less than 1.");
	
	//~ x0 = p().Mag() / BetaFrom_Gamma(w0);
//~ break;

//~ case Vec4from2::BoostMinusOne_preserve_p3:
	//~ if(w0 < real_t(0))
		//~ throw std::domain_error("Vector4: boost factor cannot be less than 1.");
	
	//~ x0 = p().Mag() / BetaFrom_GammaMinusOne(w0);
//~ break;

//~ case Vec4from2::Boost_preserve_E:
	//~ if(w0 < real_t(1))
		//~ throw std::domain_error("Vector4: boost factor cannot be less than 1.");
		
	//~ x0 = p().Mag();
	//~ p() *= BetaFrom_Gamma(w0);
//~ break;

//~ case Vec4from2::BoostMinusOne_preserve_E:
	//~ if(w0 < real_t(0))
		//~ throw std::domain_error("Vector4: boost factor cannot be less than 1.");
		
	//~ x0 = p().Mag();
	//~ p() *= BetaFrom_GammaMinusOne(w0);
//~ break;

template<typename real_t>  //++
real_t kdp::Vector4<real_t>::relDiffThreshold = std::sqrt(real_t(16)*std::numeric_limits<real_t>::epsilon());

template<typename real_t>  //++
void kdp::Vector4<real_t>::SetLengthRelDiffThreshold(real_t const newThreshold)
{
	relDiffThreshold = std::fabs(newThreshold);
}

template<typename real_t>  //++
kdp::Vector4<real_t>& kdp::Vector4<real_t>::operator+=(Vector4 const& that)
{
	this->Vector3<real_t>::operator += (that);
	this->x0 += that.x0;
	return *this;
}

template<typename real_t>  //++
kdp::Vector4<real_t>& kdp::Vector4<real_t>::operator-=(Vector4 const& that)
{
	this->Vector3<real_t>::operator -= (that);
	this->x0 -= that.x0;
	return *this;
}

template<typename real_t>  //++
kdp::Vector4<real_t>& kdp::Vector4<real_t>::operator*=(real_t const scale)
{
	this->Vector3<real_t>::operator *= (scale);
	this->x0 *= scale;
	return *this;
}

template<typename real_t>  //++
kdp::Vector4<real_t>& kdp::Vector4<real_t>::operator/=(real_t const scale)
{
	this->Vector3<real_t>::operator /= (scale);
	this->x0 /= scale;
	return *this;
}

template<typename real_t> //++
kdp::Vector4<real_t> kdp::Vector4<real_t>::operator-() const
{
	return Vector4(-x0, -x1, -x2, -x3);
}

template<typename real_t>  //++
kdp::Vector4<real_t>::operator std::vector<real_t>() const
{
	return {x0, x1, x2, x3};
}

// No cancellation protection here, since there's no expectation
// that (this) and (that) are the same vector
template<typename real_t>  //++
real_t kdp::Vector4<real_t>::Contract(Vector4 const& that) const
{
	return std::fma(this->x0, that.x0, -this->Vector3<real_t>::Dot(that));
}

// For the vector magnitude, we use cancellation protection
// (for massless particles encountering rounding error)
// Such protection can be deactivated with:
//    Vector4::SetLengthRelDiffThreshold(0.);

#include <stdio.h>
template<typename real_t>  //++
real_t kdp::Vector4<real_t>::Length2() const
{
	const real_t vec3mag2 = Vector3<real_t>::Mag2();
	
	if((relDiffThreshold > real_t(0)) and
		(std::fabs(RelDiff(Squared(x0), vec3mag2)) < relDiffThreshold))
	{
		return real_t(0);
	}
	else return std::fma(x0, x0, -vec3mag2);
}

// Return a negative distance instead of an imaginary distance
template<typename real_t>  //++
real_t kdp::Vector4<real_t>::Length() const
{
	const real_t length2 = Length2();
	return std::copysign(std::sqrt(std::fabs(length2)), length2);
}

template<typename real_t>  //+
real_t kdp::Vector4<real_t>::Length2_T() const
{
	return kdp::Diff2(x0, x3);
}

// Return a negative distance instead of an imaginary distance
template<typename real_t>  //+
real_t kdp::Vector4<real_t>::Length_T() const
{
	const real_t length2_T = Length2_T();
	return std::copysign(std::sqrt(std::fabs(length2_T)), length2_T);
}

template<typename real_t>
real_t kdp::Vector4<real_t>::Rapidity() const
{
	/* rapidity = 0.5*log((x0 + x3)/(x0 - x3)) = 0.5*log1p(2*x3/(x0-x3)
	 * 
	 * Normally, we would avoid using log1p(y) for y < 0, 
	 * because y cannot get arbitrarily close to -1,
	 * but CAN get ``arbitrarily'' close to infinity.
	 * However, 2*x3/(x0-x3) cannot actually get arbitrarily close to inf
	 * because its largeness is caused by a benign cancellation (when x3 ~= x0).
	 * Thus, we don't need to branch on the sign of x3. 
	*/

	// This is actually the definition of std::atanh(x3/x0) in glibc,
	// but we don't know what library we are linking against
	return real_t(0.5)*std::log1p((real_t(2)*x3)/(x0 - x3));
}

template<typename real_t>
real_t kdp::Vector4<real_t>::DeltaR2_rap(Vector4 const& that) const
{
	return kdp::Squared(this->Rapidity() - that.Rapidity()) + kdp::Squared(this->DeltaPhi(that));
}

template<typename real_t>
real_t kdp::Vector4<real_t>::DeltaR_rap(Vector4 const& that) const
{
	return std::sqrt(this->DeltaR2_rap(that));
}

template<typename real_t>
real_t kdp::Vector4<real_t>::DeltaR2_pseudo(Vector4 const& that) const
{
	return kdp::Squared(this->DeltaEta(that)) + kdp::Squared(this->DeltaPhi(that));
}

template<typename real_t>
real_t kdp::Vector4<real_t>::DeltaR_pseudo(Vector4 const& that) const
{
	return std::sqrt(this->DeltaR2_pseudo(that));
}

template<typename real_t>
real_t kdp::Vector4<real_t>::Beta() const
{
	return p().Mag() / x0;
}

template<typename real_t>
real_t kdp::Vector4<real_t>::Beta(kdp::Vector3<real_t> const& p3, real_t const mass)
{
									GCC_IGNORE_PUSH(-Wfloat-equal)
	if(mass < real_t(0))
		throw std::domain_error("Vector4::Beta(p3, mass): space-like particle.");
	else if (mass == real_t(0))
		return real_t(1);
	else 
	{
		real_t const pSquared = p3.Mag2();
		return std::sqrt(pSquared/(pSquared + kdp::Squared(mass)));
	}
									GCC_IGNORE_POP
}

template<typename real_t>
kdp::Vector3<real_t> kdp::Vector4<real_t>::BetaVec() const
{
	return p() / x0;
}

template<typename real_t>
template<typename convertTo>  //++
kdp::Vector4<real_t>::operator kdp::Vector4<convertTo>() const
{
	return kdp::Vector4<convertTo>(convertTo(x0), convertTo(x1), convertTo(x2), convertTo(x3));
}

//~ template<typename real_t>
//~ real_t kdp::Vector4<real_t>::BetaFrom_Gamma(real_t const gamma)
//~ {
	//~ // This is the cutoff where beta > 0.5
	//~ static constexpr real_t gamma_cutoff = std::sqrt(4./3.);
	
	//~ // beta = sqrt(1 - gamma**-2)	
	//~ if(gamma < gamma_cutoff)
	//~ {
		//~ // To get the most accurate answer from gamma, we must use Diff2
		//~ // However, the most accurate answer will come from |p| / E
		//~ return std::sqrt(kdp::Diff2(real_t(1), real_t(1)/gamma));
	//~ }
	//~ else // beta is close to one, we can get a more accurate answer 
	//~ {
		//~ // 1 - (1 - beta) is actually the most accurate beta
		//~ // 1 - beta = 1 - sqrt(1 - gamma*-2) = 1/(gamma (gamma + sqrt(gamma**2 - 1)))
		//~ return real_t(1) - real_t(1)/(gamma * 
			//~ (gamma + std::sqrt(kdp::Diff2(gamma, real_t(1)))));
	//~ }	
//~ }

//~ template<typename real_t>
//~ real_t kdp::Vector4<real_t>::BetaFrom_GammaMinusOne(real_t const gm1)
//~ {
	//~ // This is the cutoff where beta > 0.5. gamma_cutoff = sqrt(4/3)
	//~ // gm1_cutoff = sqrt(4/3) - 1 = (4/3 - 1)/(sqrt(4/3)+1) = 1/(3*sqrt(4/3) + 3)
	//~ static constexpr real_t gm1_cutoff = real_t(1)/(std::sqrt(real_t(12)) + real_t(3));
	
	//~ // beta = sqrt(1 - gamma**-2)
	//~ if(gm1 < gm1_cutoff)
		//~ return std::sqrt(gm1*(gm1 + real_t(2)))/(gm1 + real_t(1));
	//~ else // beta is close to one, we can get a more accurate answer by subtracting
	//~ {
		//~ // 1 - beta = 1 - sqrt(1 - gamma*-2) = 1/(gamma2 (1 + sqrt(1 - gamma2**-2)))
		//~ // 1 - (1 - beta) is actually the most accurate beta
		//~ return real_t(1) - real_t(1)/((gm1 + real_t(1))*
			//~ (real_t(1) + (gm1 + std::sqrt(gm1*(gm1 + real_t(2))))));
	//~ }
//~ }

// Insantiate common types
template class kdp::Vector4<double>;
template class kdp::Vector4<float>;
//template class kdp::Vector2<long double>;

// Instantiate conversion functions between common types 
// (nested template, not instantiated by struct instantiation 
template kdp::Vector4<float>::operator kdp::Vector4<double>() const;
template kdp::Vector4<double>::operator kdp::Vector4<float>() const;

////////////////////////////////////////////////////////////////////////

template<typename real_t>		
kdp::Rotate3<real_t>::Rotate3(vec3_t const& axis_in, real_t const psi):
	axis(axis_in),
	cos_psi(std::cos(psi)),
	oneMcos_psi(real_t(2)*kdp::Squared(std::sin(real_t(0.5)*psi))), // Better precision for 1 - cos(x)
	sin_psi(std::sin(psi))
{
	axis.Normalize();
}

template<typename real_t>		
kdp::Rotate3<real_t>::Rotate3(vec3_t const& u, vec3_t const& v, real_t const omega):
	axis(false) // skip initialize
{
	real_t const u_Mag2 = u.Mag2();
	real_t const v_Mag2 = v.Mag2();
													GCC_IGNORE_PUSH(-Wfloat-equal)
	if((u_Mag2 * v_Mag2) == real_t(0))
		throw std::invalid_argument("Rotate3: no rotation can be defined ... at least one null vector");
													GCC_IGNORE_POP
	
	vec3_t const u_hat = u / std::sqrt(u_Mag2);
	vec3_t const v_hat = v / std::sqrt(v_Mag2);
	
	if(u_hat.Dot(v_hat) + 1. < real_t(3) * std::numeric_limits<real_t>::epsilon())
		throw std::invalid_argument(std::string("Rotate3: u and v are antiparallel, the rotation is ambiguous.")
			 + "This ambiguity can only be broken by defining the orthogonal axis of rotation, "
			 + "in which case one should define the rotation via axis/angle=pi.");
	
	// Do interior angle manually because we need Mag2 for other things
	real_t const halfTheta = real_t(0.5)*std::atan2(std::sqrt(u.Cross(v).Mag2() / (u_Mag2 * v_Mag2)), 
		u.Dot(v) / std::sqrt(u_Mag2 * v_Mag2));
	real_t const halfOmega = real_t(0.5)*omega;
	
	vec3_t const x1 = u.Cross(v).Normalize();
	vec3_t const x2 = (u_hat + v_hat).Normalize();
	
	real_t const sin_halfOmega = std::sin(halfOmega);
	real_t const sin_halfTheta = std::sin(halfTheta);
	real_t const c_Theta = std::sqrt(
		real_t(2)*kdp::Squared(sin_halfTheta)
		+ real_t(2)*kdp::Squared(sin_halfOmega)
		+ kdp::Squared(std::sin(real_t(2)*(halfTheta + halfOmega)))
		+ kdp::Squared(std::sin(real_t(2)*(halfTheta - halfOmega))));
	//~ real_t const c_Theta = std::sqrt(real_t(3) -
		//~ std::cos(omega) - std::cos(real_t(2)*halfTheta)*(real_t(1) + std::cos(omega)));
		
	if(c_Theta < std::numeric_limits<real_t>::epsilon()) // u || v and omega = 0 -> identify matrix
	{
		axis = x2; // 
		cos_psi = real_t(1);
		sin_psi = oneMcos_psi = real_t(0);
	}
	else
	{	
		real_t const a = real_t(2) * std::cos(halfOmega) * sin_halfTheta / c_Theta;
		real_t const b = real_t(2) * sin_halfOmega / c_Theta;
			
		axis = (x1 * a + x2 * b).Normalize();
		
		vec3_t const u_perp = u_hat - axis * (u_hat.Dot(axis));	
		vec3_t const v_perp = v_hat - axis * (v_hat.Dot(axis));
		
		real_t const psi = std::copysign(u_perp.InteriorAngle(v_perp), a);
		
		cos_psi = std::cos(psi);
		oneMcos_psi = real_t(2)*kdp::Squared(std::sin(real_t(0.5)*psi));
		sin_psi = std::sin(psi);
	}
}

template<typename real_t>
typename kdp::Rotate3<real_t>::vec3_t& kdp::Rotate3<real_t>::operator()(vec3_t& vict) const
{
	/* Length correction implemented to correct length-mangling for 
	 * near-pi rotations. The length correction step only takes ~15% longer,
	 * so it doesn't hurt enough to justify removing 
	 * (especially when rotating unit vectors).
	 * 
	 * A possibly faster solution for large rotations is to compose the 
	 * rotation from two parity flips (reverse sign of two of the axis), 
	 * with a small rotation to finish the job. This only works if we do 
	 * an EVEN number of parity flips --- otherwise it is an improper rotation.
	 * 
	 * Note: when cross is small, it's from cancellation, when dot is small, it's from cancellation
	 * I don't think there is any way to get a more accurate rotation from two vectors.
	 * When the b is in the orthogonal plane, the dot with the axis
	 * has a decent amount of relative error, but that's because the dot is small.
	*/

	real_t const orig_mag2 = vict.Mag2();
	// This form used here is better than
	// 	parallel + (b - parallel)*cos_phi + ...
	// because there is one less vector addition.
	vec3_t const parallel = axis * axis.Dot(vict);
	vict = parallel * oneMcos_psi + (vict * cos_psi + axis.Cross(vict) * sin_psi);
	//~ vict = parallel  + ((vict - parallel) * cos_psi + axis.Cross(vict) * sin_psi);
	vict *= std::sqrt(orig_mag2 / vict.Mag2());
	
	return vict;
}
		
		
template<typename real_t>
typename kdp::Rotate3<real_t>::vec3_t kdp::Rotate3<real_t>::operator()(vec3_t const& b) const
{
	vec3_t copy = b;
	return (*this)(copy);
}

template<typename real_t>
real_t kdp::Rotate3<real_t>::Angle() const
{
	return std::atan2(sin_psi, cos_psi);
}

template class kdp::Rotate3<float>;
template class kdp::Rotate3<double>;

////////////////////////////////////////////////////////////////////////

template<typename real_t>
kdp::LorentzBoost<real_t>::LorentzBoost(Vector3<real_t> const& beta):
	axis(beta)
{
	axis.Normalize();
	
											GCC_IGNORE_PUSH(-Wfloat-equal)
	if(axis.Mag2() == real_t(0))
		gamma_m_1 = betaGamma = real_t(0);
	else
	{	
		// (09.02.2018 @ 11:22) Double-checked math
		real_t const beta2 = beta.Mag2();
		
		if(beta2 >= real_t(1))
			throw std::invalid_argument("kdp::LorentzBoost: cannot boost to beta = 1");
		
		real_t const oneOverGamma = std::sqrt(real_t(1) - beta2);
		gamma_m_1 = beta2 / (oneOverGamma*(real_t(1) + oneOverGamma));
		betaGamma = std::sqrt(beta2 / (real_t(1) - beta2));
	}
											GCC_IGNORE_POP
}

// (09.02.2018 @ 11:22) Double-checked math
template<typename real_t>
kdp::LorentzBoost<real_t>::LorentzBoost(Vector3<real_t> const& axis_in, 
	real_t const gamma):
	axis(axis_in),
	gamma_m_1(gamma - real_t(1)),
	betaGamma(std::sqrt(kdp::Diff2(gamma, real_t(1))))
{
	axis.Normalize();
											GCC_IGNORE_PUSH(-Wfloat-equal)
	if(gamma_m_1 < real_t(0))
		throw std::domain_error("kdp::LorentzBoost: cannot boost to gamma < 1");
	if(std::isinf(gamma_m_1))
		throw std::domain_error("kdp::LorentzBoost: cannot boost to beta = 1");
		
	if(gamma_m_1 == real_t(0))
		axis = Vector3<real_t>(); // If the boost is null, destroy the meaningless axis
	else if(axis.Mag2() == real_t(0))
		throw std::invalid_argument("kdp::LorentzBoost: null axis supplied ... which way should we boost?");
	
											GCC_IGNORE_POP
}

template<typename real_t>
kdp::Vector4<real_t>& kdp::LorentzBoost<real_t>::Boost_sign(Vector4<real_t>& victim, 
	real_t const sign) const
{
	// sign = -1 indicates a backward boost, which alters the sign of beta
	real_t const betaGamma_thisBoost = std::copysign(betaGamma, sign);
	
	// (29.07.2018 @ 21:36) Tripe-checked math and corrected sign error
	real_t const p_L = victim.p().Dot(axis);
	// CAREFUL: must alter momentum first, because x0 is used for new-p calculation 
	// (whereas p_L is previously calculated for new-x0 calculation).
	victim.p() += axis * (betaGamma_thisBoost * victim.x0 + gamma_m_1 * p_L);
	victim.x0 += gamma_m_1 * victim.x0 + betaGamma_thisBoost * p_L;
	return victim;
}

template<typename real_t>
kdp::Vector4<real_t>& kdp::LorentzBoost<real_t>::Forward(Vector4<real_t>& victim) const
{
	return this->Boost_sign(victim, real_t(1));
}

template<typename real_t>
kdp::Vector4<real_t>& kdp::LorentzBoost<real_t>::Backward(Vector4<real_t>& victim) const
{
	return this->Boost_sign(victim, real_t(-1));
}

template<typename real_t>
kdp::Vector4<real_t> kdp::LorentzBoost<real_t>::Forward(Vector4<real_t> const& orig) const
{
	auto victim = orig;
	return this->Forward(victim);
}

template<typename real_t>
kdp::Vector4<real_t> kdp::LorentzBoost<real_t>::Backward(Vector4<real_t> const& orig) const
{
	auto victim = orig;
	return this->Backward(victim);
}

template class kdp::LorentzBoost<float>;
template class kdp::LorentzBoost<double>;
