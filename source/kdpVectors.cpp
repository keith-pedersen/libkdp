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
kdp::Vector2<real_t>::Vector2(): x1(0.), x2(0.) {}

template<typename real_t>  //+++
kdp::Vector2<real_t>::Vector2(bool) {} // Don't initilize, you know what you're doing

template<typename real_t>  //+++
kdp::Vector2<real_t>::Vector2(real_t const x1_init, real_t const x2_init):
   x1(x1_init), x2(x2_init) {}
  
template<typename real_t>  //+++
kdp::Vector2<real_t>::Vector2(real_t const w1, real_t const w2, Vec2from const argFormat)
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

template<typename real_t>  //++
kdp::Vector2<real_t>& kdp::Vector2<real_t>::operator*=(real_t const scale)
{
   this->x1 *= scale;
   this->x2 *= scale;
   return *this;
}

template<typename real_t>  //--
kdp::Vector2<real_t>& kdp::Vector2<real_t>::operator/=(real_t const scale)
{
   this->x1 /= scale;
   this->x2 /= scale;
   return *this;
}

template<typename real_t>  //++
kdp::Vector2<real_t> kdp::Vector2<real_t>::operator-() const
{
   return Vector2(-x1, -x2);
}

template<typename real_t>  //+++
kdp::Vector2<real_t>& kdp::Vector2<real_t>::Normalize()
{
	real_t const length = Mag();
	if(length > 0.)
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

template<typename real_t>  //+++
real_t kdp::Vector2<real_t>::Dot(Vector2 const& that) const
{
	return std::fma(this->x2, that.x2, this->x1 * that.x1);
}

template<typename real_t>  //+++
real_t kdp::Vector2<real_t>::Cross(Vector2 const& that) const
{
   return std::fma(this->x1, that.x2,  -this->x2*that.x1);
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
template<typename convertTo>  //++
kdp::Vector2<real_t>::operator kdp::Vector2<convertTo>() const
{
   return kdp::Vector2<convertTo>(convertTo(x1), convertTo(x2));
}

// Insantiate common types
template class kdp::Vector2<double>;
template class kdp::Vector2<float>;
//template class kdp::Vector2<long double>;

// Insantiate conversion functions between common types
// (nested template, not instantiated by struct instantiation 
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
kdp::Vector3<real_t>::Vector3(): Vector2<real_t>(), x3(0.) {}

template<typename real_t>  //+++
kdp::Vector3<real_t>::Vector3(bool): Vector2<real_t>(false) {} // Don't initilize, you know what you're doing

template<typename real_t> //+++
kdp::Vector3<real_t>::Vector3(real_t const x1_init, real_t const x2_init, real_t const x3_init):
   Vector2<real_t>(x1_init, x2_init), x3(x3_init) {}
   
template<typename real_t> //+++
kdp::Vector3<real_t>::Vector3(real_t const w1, real_t const w2, real_t const w3, 
	Vec3from const argFormat)
{
   switch(argFormat)
   {
      case Vec3from::LengthEtaPhi:
      {
			// w1 -> length
			// w2 -> eta
			// w3 -> phi
			real_t const coshEta = std::cosh(w2);
			
			x1 = (w1 * std::cos(w3))/coshEta;
			x2 = (w1 * std::sin(w3))/coshEta;
			x3 = w1 * std::tanh(w2);
      }
      break;
      
      case kdp::Vec3from::LengthThetaPhi:
      {
			// w1 -> length
			// w2 -> theta
			// w3 -> phi
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

template<typename real_t> //++
kdp::Vector3<real_t>& kdp::Vector3<real_t>::operator-=(Vector3 const& that)
{
   this->Vector2<real_t>::operator -= (that);
   this->x3 -= that.x3;
   return *this;
}

template<typename real_t> //++
kdp::Vector3<real_t>& kdp::Vector3<real_t>::operator*=(real_t const scale)
{
   this->Vector2<real_t>::operator *= (scale);
   this->x3 *= scale;
   return *this;
}

template<typename real_t> //++
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
	if(length > real_t(0.))
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
   // Return +/- infinity for a particle parallel to z-axis
   // Return nan for stationary particle (but what else is more accurate)

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
      (real_t(2.)*(absZ*(std::sqrt(std::fma(absZ, absZ, xT2)) + absZ)))/xT2), x3);
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

//~ template<typename real_t> //++
//~ real_t kdp::Vector3<real_t>::Tprojection() const
//~ {
   //~ return std::sqrt(Tprojection2());
//~ }

//~ template<typename real_t> //++
//~ real_t kdp::Vector3<real_t>::Tprojection2() const
//~ {
   //~ const real_t T2 = Vector2<real_t>::Mag2();
   //~ return T2 / std::fma(x3, x3, T2);
//~ }

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
kdp::Vector4<real_t>::Vector4(): Vector3<real_t>(), x0(0.) {}

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
   switch(argFormat)
   {
      case kdp::Vec4from4::EnergyEtaPhiM:
      {
         {
            // w0 -> E,
            // w1 -> eta
            // w2 -> phi
            // w3 -> mass
            real_t const momentum = (w3 == 0.) ? w0 : std::sqrt((w0 - w3)*(w0 + w3));
            x() = Vector3<real_t>(momentum, w1, w2, Vec3from::LengthEtaPhi);
         }
         x0 = w0;
      }
      break;
   }
}

template<typename real_t>
kdp::Vector4<real_t>::Vector4(Vector3<real_t> const& x_in): 
	Vector3<real_t>(x_in),	x0(this->x().Mag()) {}

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
		
		case Vec4from2::Boost_preserve_p3:
			if(w0 < real_t(1))
				throw std::domain_error("Vector4: boost factor cannot be less than 1.");
			
			x0 = p().Mag() / BetaFrom_Gamma(w0);
		break;
		
		case Vec4from2::BoostMinusOne_preserve_p3:
			if(w0 < real_t(1))
				throw std::domain_error("Vector4: boost factor cannot be less than 1.");
			
			x0 = p().Mag() / BetaFrom_GammaMinusOne(w0);
		break;
		
		case Vec4from2::Boost_preserve_E:
			if(w0 < real_t(1))
				throw std::domain_error("Vector4: boost factor cannot be less than 1.");
				
			x0 = p().Mag();
			p() *= BetaFrom_Gamma(w0);
		break;
		
		case Vec4from2::BoostMinusOne_preserve_E:
			if(w0 < real_t(1))
				throw std::domain_error("Vector4: boost factor cannot be less than 1.");
				
			x0 = p().Mag();
			p() *= BetaFrom_GammaMinusOne(w0);
		break;
		
		default:
			throw std::runtime_error("Vector4::Vecto4. Unrecognized Vec4from2 enum");
		break;
	}
}

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
   
   if((relDiffThreshold > 0.) and
		(std::fabs(RelDiff(Squared(x0), vec3mag2)) < relDiffThreshold))
   {
      return 0.;
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

template<typename real_t>
real_t kdp::Vector4<real_t>::Rapidity() const
{
   // Return +/- infinity for massless particle parallel to z-axis
   // Return nan for totally null or space-like particle (but what else is more accurate)
   
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
   return real_t(0.5)*std::log1p((real_t(2.)*x3)/(x0 - x3));
}

template<typename real_t>
template<typename convertTo>  //++
kdp::Vector4<real_t>::operator kdp::Vector4<convertTo>() const
{
   return kdp::Vector4<convertTo>(convertTo(x0), convertTo(x1), convertTo(x2), convertTo(x3));
}

template<typename real_t>
real_t kdp::Vector4<real_t>::BetaFrom_Gamma(real_t const gamma)
{
	// This is the cutoff where beta > 0.5
	static constexpr real_t gamma_cutoff = std::sqrt(4./3.);
	
	// beta = sqrt(1 - gamma**-2)	
	if(gamma < gamma_cutoff)
	{
		// To get the most accurate answer from gamma, we must use Diff2
		// However, the most accurate answer will come from |p| / E
		return std::sqrt(kdp::Diff2(real_t(1), real_t(1)/gamma));
	}
	else // beta is close to one, we can get a more accurate answer 
	{
		// 1 - (1 - beta) is actually the most accurate beta
		// 1 - beta = 1 - sqrt(1 - gamma*-2) = 1/(gamma (gamma + sqrt(gamma**2 - 1)))
		return real_t(1) - real_t(1)/(gamma * 
			(gamma + std::sqrt(kdp::Diff2(gamma, real_t(1)))));
	}	
}

template<typename real_t>
real_t kdp::Vector4<real_t>::BetaFrom_GammaMinusOne(real_t const gm1)
{
	// This is the cutoff where beta > 0.5. gamma_cutoff = sqrt(4/3)
	// gm1_cutoff = sqrt(4/3) - 1 = (4/3 - 1)/(sqrt(4/3)+1) = 1/(3*sqrt(4/3) + 3)
	static constexpr real_t gm1_cutoff = real_t(1)/(std::sqrt(real_t(12)) + real_t(3));
	
	// beta = sqrt(1 - gamma**-2)
	if(gm1 < gm1_cutoff)
		return std::sqrt(gm1*(gm1 + real_t(2)))/(gm1 + real_t(1));
	else // beta is close to one, we can get a more accurate answer by subtracting
	{
		// 1 - beta = 1 - sqrt(1 - gamma*-2) = 1/(gamma2 (1 + sqrt(1 - gamma2**-2)))
		// 1 - (1 - beta) is actually the most accurate beta
		return real_t(1) - real_t(1)/((gm1 + real_t(1))*
			(real_t(1) + (gm1 + std::sqrt(gm1*(gm1 + real_t(2))))));
	}
}

template<typename real_t>
real_t kdp::Vector4<real_t>::BetaFrom_Mass_pSquared(real_t const mass, real_t const pSquared)
{
	return std::sqrt(pSquared/(pSquared + kdp::Squared(mass)));
}

//~ template <typename real_t>
//~ kdp::Vector4<real_t> kdp::MasslessVec4_EnergyEtaPhi(real_t const E, real_t const eta, real_t const phi)
//~ {
   //~ const real_t coshEta = std::cosh(eta);
   
   //~ return Vec4(E, (E*std::cos(phi))/coshEta, (E*std::sin(phi))/coshEta, E*std::tanh(eta));
//~ }

// Insantiate common types
template class kdp::Vector4<double>;
template class kdp::Vector4<float>;
//template class kdp::Vector2<long double>;

// Insantiate conversion functions between common types
// (nested template, not instantiated by struct instantiation 
template kdp::Vector4<float>::operator kdp::Vector4<double>() const;
template kdp::Vector4<double>::operator kdp::Vector4<float>() const;

//~ template kdp::Vector4<double> kdp::MasslessVec4_EnergyEtaPhi(double const E, double const eta, double const phi);
//~ template kdp::Vector4<float> kdp::MasslessVec4_EnergyEtaPhi(float const E, float const eta, float const phi);


//! @brief Construct the object that takes u> to v> 
		
template<typename real_t>		
kdp::Rotate3<real_t>::Rotate3(vec3_t const& u, vec3_t const& v):
	axis_NN(u.Cross(v)), 
	axis_mag2(axis_NN.Mag2()),
	uv_mag2(u.Mag2() * v.Mag2()),
	cos_phi(u.Dot(v) / std::sqrt(uv_mag2))
{
	if(uv_mag2 == real_t(0))
		throw std::invalid_argument("Rotate2: no rotation can be defined ... at least one null vector");
}

template<typename real_t>
typename kdp::Rotate3<real_t>::vec3_t& kdp::Rotate3<real_t>::operator()(vec3_t& b) const
{
	// Only possible when (u == +/- v), so either identity or parity
	if((axis_mag2 == real_t(0)) and (cos_phi < real_t(0)))
		b = -b;
	else
	{
		vec3_t const parallel = axis_NN * (axis_NN.Dot(b)/axis_mag2);
		double const orig_mag2 = b.Mag2();
		// This form is better because than
		// 	parallel + (b - parallel)*cos_phi + ...
		// because there is one less vector addition.
		// Note: when cross is small, it's from cancellation, when dot is small, it's from cancellation
		// I don't think there is any way to get a more accurate rotation from two vectors.
		// When the b is in the orthogonal plane, the dot with the axis
		// has a decent amount of relative error, but that's because the dot is small.					
		b = parallel*(real_t(1) - cos_phi) + (b*cos_phi + axis_NN.Cross(b)/std::sqrt(uv_mag2));
		b *= std::sqrt(orig_mag2 / b.Mag2());
	}
	
	return b;
}
		
template<typename real_t>
typename kdp::Rotate3<real_t>::vec3_t kdp::Rotate3<real_t>::operator()(vec3_t const& b) const
{
	vec3_t copy = b;
	return (*this)(copy);
}
		
template<typename real_t>
typename kdp::Rotate3<real_t>::vec3_t kdp::Rotate3<real_t>::Axis() const
{
	return axis_NN/std::sqrt(axis_mag2);
}

template class kdp::Rotate3<float>;
template class kdp::Rotate3<double>;

typedef kdp::Rotate3<double> Rot3;
typedef kdp::Rotate3<float> Rot3_f;

////////////////////////////////////////////////////////////////////////

template<typename real_t>
kdp::LorentzBoost<real_t>::LorentzBoost(Vector4<real_t> const& target, 
	bool const CMtoTarget):
// Don't init the 4-vectors
lambda_0(false), lambda_1(false), lambda_2(false), lambda_3(false)
{
	// See ArbitraryBoost.pdf for an explanation of the math
	real_t const p2 = target.p().Mag2();
	real_t const mass2 = target.Length2();
	
	if(mass2 <= 0.)
		throw std::runtime_error("LorentzBoost: boost matrix does not exist for massless of space-like Vector4");
	
	real_t const denom = (mass2 + std::sqrt(mass2 * target.x0 * target.x0));
	
	lambda_0.x0 = p2 / denom;
			
	lambda_0.p() = (CMtoTarget ? target.p() : -target.p());
	lambda_0.p() /= std::sqrt(mass2);
	
	lambda_1.x0 = lambda_0.x1;
	lambda_2.x0 = lambda_0.x2;
	lambda_3.x0 = lambda_0.x3;
	
	lambda_1.x1 = (target.p().x1 * target.p().x1) / denom;
	lambda_2.x1 = lambda_1.x2 = (target.p().x1 * target.p().x2) / denom;
	lambda_3.x1 = lambda_1.x3 = (target.p().x1 * target.p().x3) / denom;
	
	lambda_2.x2 = (target.p().x2 * target.p().x2) / denom;
	lambda_3.x2 = lambda_2.x3 = (target.p().x2 * target.p().x3) / denom;
	
	lambda_3.x3 = (target.p().x3 * target.p().x3) / denom;
}

template<typename real_t>
kdp::Vector4<real_t> kdp::LorentzBoost<real_t>::operator()(Vector4<real_t> const& victim) const
{
	// First construct the 4-vector shift
	kdp::Vector4<real_t> boosted = lambda_0 * victim.x0;
	boosted += lambda_1 * victim.x1;
	boosted += lambda_2 * victim.x2;
	boosted += lambda_3 * victim.x3;
	
	// Then we add the shift to the original vector
	boosted += victim;
	
	return boosted;
}

template class kdp::LorentzBoost<float>;
template class kdp::LorentzBoost<double>;
