// Copyright (C) 2014-2018 by Keith Pedersen (Keith.David.Pedersen@gmail.com)

// Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
// The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

#ifndef KDP_VECTORS
#define KDP_VECTORS

// includes must be outside of namespace
#include "kdpTools.hpp"

namespace kdp
{
	
GCC_IGNORE_PUSH(-Wpadded)

/*! @file kdpVectors.hpp
 *  @brief Defines 2-vectors and 3-vectors (Euclidean metric) and 
 *  4-vectors (Minkowski metric, (+,-,-,-)), along with various operations.
 * 
 *  These classes take unusual care to provide the highest level of floating-point stability.
 *  Additionally, vectors are defined as template classes, for easy specification
 *  of the floating-point types storing the components.
 *  
 *  \note There is no equality operation (==) for these vector classes
 *  because exact equality is not safe for floating-point calculations.
 *  A good test for equality is (aVec - bVec).Mag() / (aVec + bVec).Mag() < smallValue.
 *  
 *  \note (Design decisions): Vec4 inherits Vec3, which inherits Vec2.
 *  However, the inheritance is protected (i.e. hidden from the public), 
 *  which abuses the  Vec3 "is a" Vec2 test for inheritance.
 *  In fact, "composition" (Vector3 "has a" Vector2) makes the most sense for vectors,
 *  but creates cumbersome access to the vector components, complicating math expressions.
 *  We choose protected inheritance because it simplifies expressions and
 *  protects against unwanted conversion (Vec2 t2 = Vec3(1, 2, 3)).
 *  Additionally, pretty much all functions need to be redefined
 *  for each vector class. If they have the same name,
 *  the lower versions are simply hidden (still accessible from within the class, 
 *  but not from outside). Vec4 explicitly renames its Mag and Dot functions, 
 *  to remind the user of the different metric.
 *  With public inheritance, we would have to delete Mag and Dot,
 *  but his would violate the C++11 standard:
 *  	A function with a deleted definition shall not
 *  	override a function that does not have a deleted definition.
 *
 *  @author Copyright (C) 2014-2018 Keith Pedersen (Keith.David.Pedersen@gmail.com)
*/

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    __     __        _            ____  
    \ \   / /__  ___| |_ ___  _ _|___ \ 
     \ \ / / _ \/ __| __/ _ \| '__|__) |
      \ V /  __/ (__| || (_) | |  / __/ 
       \_/ \___|\___|\__\___/|_| |_____|
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

//! @brief When a Vector2 is constructed from two generic d.o.f., this enum identifies those d.o.f.
enum class Vec2from {LengthPhi}; // Add more options as needed
// We declare the enums outside of the class which uses them
// so they do not require template specialization
                                          
//! @brief Cartesian, Euclidean 2-vector
template<typename real_t>
struct Vector2
{
   real_t x1;
   real_t x2;

   Vector2(); //!< @brief Construct a null vector (0, 0).
   explicit Vector2(bool); //!< @brief Don't initialize the vector (regardless of bool's value).
   explicit Vector2(real_t const x1_init, real_t const x2_init);
   //! @brief Construct from two generic d.o.f., identified via \p argFormat.
   explicit Vector2(real_t const w1_init, real_t const w2_init, 
		Vec2from const argFormat);
	virtual ~Vector2() {}

   //! @brief Convert the internal real_t floating-point type.
   template<typename convertTo>
   operator Vector2<convertTo>() const;

   Vector2& operator+=(Vector2 const&);
   Vector2& operator-=(Vector2 const&);
   Vector2& operator*=(real_t const);
   Vector2& operator/=(real_t const);

   inline Vector2 operator+(Vector2 const& that) const {return Vector2(*this) += that;}
   inline Vector2 operator-(Vector2 const& that) const {return Vector2(*this) -= that;}
   inline Vector2 operator*(real_t const scale) const {return Vector2(*this) *= scale;}
   inline Vector2 operator/(real_t const scale) const {return Vector2(*this) /= scale;}

   Vector2 operator-() const; //!< @brief Reverse all coordinates (180 degree rotation).
   
   Vector2& Normalize(); //!< @brief Scale to unit length.

   real_t Mag()  const; //!< @brief The vector's Euclidean magnitude.
   real_t Mag2() const; //!< @brief The vector's squared Euclidean magnitude.
   
   real_t Phi() const; //!< @brief Vec2 = [cos(phi), sin(phi)]
   real_t DeltaPhi(Vector2 const&) const;//<! @brief The azimuthal angle between two vectors (in radians).
   
   real_t Dot          (Vector2 const&) const; //!< @brief The dot/inner product.
   real_t Cross        (Vector2 const&) const; //!< @brief The cross/outer product.
   real_t InteriorAngle(Vector2 const&) const; //!< @brief The angle between two vectors (in radians).
};

// Create nicknames for common Vector2 types
using Vec2 = Vector2<double>;
using Vec2_f = Vector2<float>;

//! @brief The interior angle between two azimuthal angles from [-Pi, Pi]
template<typename real_t>
real_t DeltaPhi(real_t const phi_1, real_t const phi_2)
{
	return real_t(M_PI - std::fabs(M_PI - std::fabs(phi_1 - phi_2)));
}

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    __     __        _            _____ 
    \ \   / /__  ___| |_ ___  _ _|___ / 
     \ \ / / _ \/ __| __/ _ \| '__||_ \ 
      \ V /  __/ (__| || (_) | |  ___) |
       \_/ \___|\___|\__\___/|_| |____/ 
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/                                      

//! @brief When a Vector3 is constructed from three generic d.o.f., this enum identifies those d.o.f.
enum class Vec3from {LengthEtaPhi, LengthThetaPhi}; // Add more options as needed

//! @brief Cartesian, Euclidean 3-vector.
template<typename real_t>
struct Vector3 : protected Vector2<real_t>
{
   // Give the user access to the (T)ransverse (x1, x2) degrees of freedom.
   using Vector2<real_t>::x1;
   using Vector2<real_t>::x2;
   using Vector2<real_t>::Phi;
   real_t x3; //!< @brief The longitudenal component (x1 and x2 being the transverse).

   Vector3(); //!< @brief Construct a null vector (0, 0, 0).
   explicit Vector3(bool); //!< @brief Don't initialize the vector (regardless of bool's value).
   explicit Vector3(real_t const x1_init, real_t const x2_init, real_t const x3_init);
   //! @brief Construct from three generic d.o.f., identified via \p argFormat.
   explicit Vector3(real_t const, real_t const, real_t const, Vec3from const argFormat);
	virtual ~Vector3() {}

   //! @brief Convert the internal real_t floating-point type.
   template<typename convertTo>
   operator Vector3<convertTo>() const;

   Vector3& operator+=(Vector3 const&);
   Vector3& operator-=(Vector3 const&);
   Vector3& operator*=(real_t const);
   Vector3& operator/=(real_t const);
   
   Vector3& Normalize(); //!< @brief Scale to unit length.

   /*! @brief Access the (T)ransverse (x1, x2) degrees of freedom.
    * 
    *  e.g. vec3.T().Phi() is the phase of the transverse component.
   */ 
   inline Vector2<real_t> &      T()       {return *this;}
   inline Vector2<real_t> const& T() const {return *this;} 
   // Inlining access removes any speed penalties for this indirection

   inline Vector3 operator+(Vector3 const& that) const {return Vector3(*this) += that;}
   inline Vector3 operator-(Vector3 const& that) const {return Vector3(*this) -= that;}
   inline Vector3 operator*(real_t const scale) const {return Vector3(*this) *= scale;}
   inline Vector3 operator/(real_t const scale) const {return Vector3(*this) /= scale;}

   Vector3 operator-() const; //!< @brief Reverse all coordinates (Parity operation).
   
   // Redefine the Vector2 worker functions for 3-space
   // (note that Cross has a different return type)
   real_t Mag()      const; //!< @brief The vector's Euclidean magnitude.
   real_t Mag2()     const; //!< @brief The vector's squared Euclidean magnitude.
   
   /*! @brief Pseudorapidity
    * 
    *  \f$ \eta = {\rm arctanh}\left(\frac{\tt vec.x3}{\tt vec.T().Mag()}\right) \f$
    * 
    *  \return (+/- infinity) for a particle parallel to z-axis 
    *  (the only way to define +/- z^ using Vec3from::LengthEtaPhi).
    *  \return (nan) for stationary particle (but what else is more accurate?).
   */
   real_t Eta() const;
   real_t DeltaEta(Vector3 const&) const;
   real_t Theta() const; //!< @brief The polar angle with the +z-axis
   
   // We need to manually instantiate these because "that"
   // cannot be automatically converted to Vector2 outside the class
   // (whereas "this" can be with the using declaration, e.g. Phi()).
   real_t DeltaPhi(Vector3 const& that) const {return this->DeltaPhi(that);}

   real_t Dot          (Vector3 const&) const; //!< @brief The dot/inner product.
   Vector3 Cross       (Vector3 const&) const; //!< @brief The cross/outer product (itself a Vector3).
   real_t InteriorAngle(Vector3 const&) const; //!< @brief The angle between two vectors (in radians). 
   
};

// Create nicknames for common real types
using Vec3 = Vector3<double>;
using Vec3_f = Vector3<float>;

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    __     __        _             _  _   
    \ \   / /__  ___| |_ ___  _ __| || |  
     \ \ / / _ \/ __| __/ _ \| '__| || |_ 
      \ V /  __/ (__| || (_) | |  |__   _|
       \_/ \___|\___|\__\___/|_|     |_|  
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/  

//! @brief When a Vector4 is constructed from four generic d.o.f., this enum identifies those d.o.f.
enum class Vec4from4 {EnergyEtaPhiMass}; // Add more options as needed

//! @brief When a Vector4 is constructed from a Vector3 and one d.o.f., this enum identifies those d.o.f.
enum class Vec4from2 {Energy, Time, Mass, Length};
	
/*! @brief Cartesian, Minkowski 4-vector
 * 
 * Vector4 is defined in Minkowski space, not Euclidiean space
 * To prevent metric confusion, the Mag and Dot functions are renamed versus Vector3:
 * 
 * 	Vector3::Dot -> Vector4::Contract
 * 	Vector3::Mag -> Vector4::Length
 * 
 * There is also no Cross or InteriorAngle
*/
template <typename real_t>
struct Vector4 : private Vector3<real_t>
{
   // Give the user access to the spatial (x1, x2, x3) degrees of freedom
   using Vector3<real_t>::x1;
   using Vector3<real_t>::x2;
   using Vector3<real_t>::x3;
   real_t x0; //!< @brief The energy/time.
   
   // Give access to unambiguous scalar properties
   using Vector3<real_t>::Phi;
   using Vector3<real_t>::DeltaPhi;
   using Vector3<real_t>::Eta;
   using Vector3<real_t>::DeltaEta;
   using Vector3<real_t>::Theta;   
   
   Vector4(); //!< @brief Construct a null vector (0, 0, 0, 0)
   explicit Vector4(bool); //!< @brief Don't initialize the vector (regardless of bool's value)
   explicit Vector4(real_t const x0, real_t const x1, real_t const x2, real_t const x3);
   //! @brief Construct from four generic d.o.f., identified via \p argFormat
   explicit Vector4(real_t const, real_t const, real_t const, real_t const, Vec4from4 const argFormat);
   //! @brief Construct from a spatial 3-vector and one additional d.o.f., identified via \p argFormat
   explicit Vector4(real_t const w0, Vector3<real_t> const& xORp, Vec4from2 const w0type);
   explicit Vector4(Vector3<real_t> const& x); //!< @brief Construct a light-like 4-vector (e.g. zero interval).
   virtual ~Vector4() {}

   //! @brief Convert the internal real_t floating-point type
   template<typename convertTo>
   operator Vector4<convertTo>() const;

   Vector4& operator+=(Vector4 const&);
   Vector4& operator-=(Vector4 const&);
   Vector4& operator*=(real_t const);
   Vector4& operator/=(real_t const);

   /*! @brief Access the spatial/Euclidean (x1, x2, x3) degrees of freedom
    * 
    *  e.g. vec3.T().Mag() is the length of 3-vector.
   */ 
   inline Vector3<real_t>&        x()       {return *this;}
   inline Vector3<real_t> const&  x() const {return *this;}
   // Inlining access removes any speed penalties for this indirection
	
	/*! @brief Access the spatial/Euclidean (x1, x2, x3) degrees of freedom
    * 
    *  e.g. vec3.T().Mag() is the length of 3-vector.
   */ 
   inline Vector3<real_t>&        p()       {return *this;}
   inline Vector3<real_t> const&  p() const {return *this;}
   // Inlining access removes any speed penalties for this indirection

   inline Vector4 operator+(Vector4 const& that) const {return Vector4(*this) += that;}
   inline Vector4 operator-(Vector4 const& that) const {return Vector4(*this) -= that;}
   inline Vector4 operator*(real_t const scale) const {return Vector4(*this) *= scale;}
   inline Vector4 operator/(real_t const scale) const {return Vector4(*this) /= scale;}
   
   operator std::vector<real_t>() const; //!@brief Convert to a std::vector.

   Vector4 operator-() const; //!< @brief Reverse all signs
   
   // We need to manually instantiate these because "that"
   // cannot be automatically converted to Vector3 outside the class
   // (whereas "this" can be with the using declaration, e.g. Phi()).
   real_t DeltaPhi(Vector4 const& that) const {return this->DeltaPhi(that);}
   real_t DeltaEta(Vector4 const& that) const {return this->DeltaEta(that);}
     
   //! @brief Square root of Length2() (negative if Length2() is negative, instead of imaginary).
   real_t Length() const; 
   //! @brief \f$ x^\mu x_\mu \f$ (\note see SetLengthRelDiffThreshold)
   real_t Length2() const;
   
   //! @brief Square root of Length2_T()  (negative if Length2_T() is negative, instead of imaginary).
   real_t Length_T() const;   
   //! @brief Squared transverse length: \f$ x_0^2-x_3^2 \f$ (\note does not use threshold, see SetLengthRelDiffThreshold()).
   real_t Length2_T() const;
   
   //! @brief Square root of Mass2() (negative if Mass2() is negative, instead of imaginary).
   real_t Mass() const {return Length();}
   //! @brief \f$ p^\mu p_\mu \f$ (\note see SetLengthRelDiffThreshold)
   real_t Mass2() const {return Length2();}
   
   //! @brief Square root of Mass2_T() (negative if Mass2_T() is negative, instead of imaginary).
   real_t Mass_T() const {return Length_T();}
   //! @brief Squared transverse mass: \f$ m^2 + p_T^2 \f$ (\note does not use threshold, see SetLengthRelDiffThreshold()).
   real_t Mass2_T() const {return Length2_T();}   
   
   real_t Contract(Vector4 const&) const; //!< @brief \f$ p^\mu q_\mu \f$
   
   real_t Beta() const; //!< @brief Speed \f$ \beta = |\vec{p}|/E \f$
   Vector3<real_t> BetaVec() const; //!< @brief Speed \f$ \vec{beta} = |vec{p}/E \f$
   static real_t Beta(kdp::Vector3<real_t> const& p3, real_t const mass); //!< @brief More accurate speed from 3-momentum and mass
   
   /*! @brief Rapidity. Difference in rapidity \f$ \Delta y \$ is invariant to longitudinal boosts.
    * 
    *  \f$ y = \text{arctanh}\left(\frac{\tt vec.x3}{\tt vec.x0}\right) \f$
    * 
    *  \return (+/- infinity) for a massless particle parallel to z-axis.
    *  \return (nan) if |x3| > |x0| (but what else is more accurate?).
   */
   real_t Rapidity() const;
   
   /*! @brief An interior "angle" invariant to longitudinal boosts: 
    *  \f$ \Delta R = \Delta y^2 + \Delta \phi^2 \f$
    * 
    *  \warning Uses rapidity y.
   */ 
   real_t DeltaR2_rap(Vector4 const&) const;
   real_t DeltaR_rap(Vector4 const&) const;
   
   /*! @brief An geometric interior "angle" which, for \em massless particles, 
    * is invariant to longitudinal boosts: \f$ \Delta R = \Delta \eta^2 + \Delta \phi^2 \f$
    * 
    *  \warning Uses pseudorapidity \f$ \eta \f$.
   */
   real_t DeltaR2_pseudo(Vector4 const&) const;
   real_t DeltaR_pseudo(Vector4 const&) const;
   
   /*! @brief The Length() of a light-like 4-vector should be zero.
    * 
    *  However, due to rounding error and floating-point cancellation,
    *  we may not find identically zero in practice.
    *  To protect against such effects, we can look at the relative difference:
    * 
    *  	RelDiff  = (E^2-p^2)/(E^2+p^2)
    *  
    *  When this falls below the threshold, we will round it to zero.
    *  This threshold defaults to 16*MachineEpsilon, but can also
    *  be changed by the user. To disable this feature, simply set it to zero.
   */ 
   static void SetLengthRelDiffThreshold(const real_t newThreshold);

   protected:
      static real_t relDiffThreshold;
};

// Create nicknames for common real_t types
using Vec4 = Vector4<double>;
using Vec4_f = Vector4<float>;

/*! @brief An active, Right-handed (RH) rotation of a 3-vector about 
 *  an axis by an angle psi via Rodrigues' rotation formula.
 * 
 *  \note This class has been validated, and 99% of errors are smaller than O(1e-14).
 *  There is probably room for improvement, but we table this for now.
*/ 
template<typename real_t>
class Rotate3
{
	public: 
		using vec3_t = kdp::Vector3<real_t>;
	
	private:
		vec3_t axis;
		real_t cos_psi; 
		real_t oneMcos_psi; //! @brief 1-cos(psi)
		real_t sin_psi;
			
	public: 
		/*! @brief Construct the rotation that actively takes \p u to \p v,
		 *  with a right-handed post-rotation of \p omega degrees about \p v.
		 *  
		 *  The post-rotation is necessary to fully specify the rotation;
		 *  "the rotation which takes u> to v>" is ambiguous because it 
		 *  only defines 2 of 3 Euler angles. There are actually an
		 *  infinite number of rotations which take u> to v>, 
		 *  but only one for a given \p omega. For a detailed mathematical 
		 *  explanation, see Section A.3.1 in:
		 *  "EXPANDING THE HEP FRONTIER WITH BOOSTED B-TAGS AND THE QCD POWER SPECTRUM"
		 * 
		 *  \note The identity operation only occurs when \p u || \p v and sin(omega) == 0
		 *  \note If v and u are antiparallel, then 
		 *  the rotation is ambiguous without supplying the axis; 
		 *  use the other constructor to define psi = Pi with a known axis.
		 * 
		 *  \throws Throws invalid_argument when v and u are antiparallel.
		 */ 
		Rotate3(vec3_t const& u, vec3_t const& v, real_t const omega);
		
		//! @brief Construct a RH rotation which actively rotates vectors 
		//! about \p axis by angle \p psi (\p axis is automatically normalized).
		Rotate3(vec3_t const& axis_in, real_t const psi);
			
		vec3_t& operator()(vec3_t& victim) const; //! @brief Rotate victim
		vec3_t operator()(vec3_t const& b) const; //! @brief Rotate victim
		
		vec3_t const& Axis() const {return axis;} //! @brief The (normalized) axis of rotation
		real_t Angle() const; //! @brief The angle of RH rotation.
};

typedef kdp::Rotate3<double> Rot3;
typedef kdp::Rotate3<float> Rot3_f;

/*! @brief An object that boosts along Axis() by boost factor Gamma()
 *  (the incoming \p axis is automatically normalized).
 * 
 *  This is done using the boost equivalent of the Rodrigues formula.
 *  For a detailed mathematical explanation, see Section A.3.2 in:
 *  "EXPANDING THE HEP FRONTIER WITH BOOSTED B-TAGS AND THE QCD POWER SPECTRUM"
*/
template <typename real_t>
class LorentzBoost
{
	private:		
		Vector3<real_t> axis;
		real_t gamma_m_1;
		real_t betaGamma;
		
		//! @brief One internal boost function; sign = 1 goes forward, sign = -1 goes backkward.
		Vector4<real_t>& Boost_sign(Vector4<real_t>& victim, real_t const sign) const;
		
	public:
		//! @brief Construct the boost from speed vector \p beta
		//! (boost speed beta = |beta>| in the direction of beta^)
		LorentzBoost(Vector3<real_t> const& beta);
		
		//! @brief Construct the boost along \p axis with boost factor \p gamma
		LorentzBoost(Vector3<real_t> const& axis, real_t const gamma);
		
		Vector3<real_t> const& Axis() const {return axis;} //!< @brief normalized boost axis
		real_t Gamma() const {return gamma_m_1 + real_t(1);}
		real_t Beta() const {return betaGamma / Gamma();}
		real_t Rapidity() const {return std::asinh(betaGamma);}
		Vector3<real_t> BetaVec() const {return axis*Beta();}
				
		//! @brief Boost the \p victim by Beta()
		Vector4<real_t>& Forward(Vector4<real_t>& victim) const;
		//! @brief Boost the \p victim by Beta()
		Vector4<real_t> Forward(Vector4<real_t> const& orig) const;
		
		//! @brief Boost the \p victim by -Beta()
		Vector4<real_t>& Backward(Vector4<real_t>& victim) const;
		//! @brief Boost the \p victim by -Beta()
		Vector4<real_t> Backward(Vector4<real_t> const& orig) const;
		
		//! @brief Boost a massless \p victim by Beta()
		Vector3<real_t>& Forward(Vector3<real_t>& victim) const;
		//! @brief Boost a massless \p victim by Beta()
		Vector3<real_t> Forward(Vector3<real_t> const& orig) const;
		
		//! @brief Boost a massless \p victim by -Beta()
		Vector3<real_t>& Backward(Vector3<real_t>& victim) const;
		//! @brief Boost a massless \p victim by -Beta()
		Vector3<real_t> Backward(Vector3<real_t> const& orig) const;
		
};

typedef kdp::LorentzBoost<double> Boost;
typedef kdp::LorentzBoost<float> Boost_f;

GCC_IGNORE_POP

} // End namespace

#endif
