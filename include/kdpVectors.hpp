#ifndef KDP_VECTORS
#define KDP_VECTORS

// includes must be outside of namespace
#include "kdpTools.hpp"

namespace kdp
{

/*! @file kdpVectors.hpp
 *  @brief Defines 2-vectors and 3-vectors (Euclidean metric) and 
 *  4-vectors (Minkowski metric, (+,-,-,-)), along with various operations.
 * 
 *  Vectors are defined as template classes, for easy specification
 *  of the floating-point types storing the components.
 *  
 *  \note There is no equality operation (==) for these vector classes
 *  because exact equality is not safe for vector expressions.
 *  A good test for equality is (aVec - bVec).Mag() / (aVec + bVec).Mag(),
 *  which should be small. Of course, this doesn't work for 4-vectors.
 *  
 *  \note Vec4 inherits Vec3, which inherits Vec2.
 *  However, the inheritence is protected (i.e. hidden from the public), 
 *  which abuses the  Vec3 "is a" Vec2 test for inheritence.
 *  In fact, "composition" (Vector3 "has a" Vector2) makes the most for vectors,
 *  but makes for cumbersome access to the vector components, complicating math.
 *  Instead, we choose inheritence which is protected, not public,
 *  because it prevents unwanted conversion (Vec t2 = Vec3(1, 2, 3)).
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
 *  @author Copyright (C) 2018 Keith Pedersen (Keith.David.Pedersen@gmail.com)
 *  @date JUN 2016, JAN 2018
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
// We declare the enums outside of the class which uese them
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

   Vector2& operator+=(const Vector2&);
   Vector2& operator-=(const Vector2&);
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

   real_t Dot          (Vector2 const&) const; //!< @brief The dot/inner product.
   real_t Cross        (Vector2 const&) const; //!< @brief The cross/outer product.
   real_t InteriorAngle(Vector2 const&) const; //!< @brief The angle between two vectors (in radians).
};

// Create nicknames for common Vector2 types
using Vec2 = Vector2<double>;
using Vec2_f = Vector2<float>;

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

   Vector3& operator+=(const Vector3&);
   Vector3& operator-=(const Vector3&);
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

   Vector3 operator-() const; //! @brief Reverse all coordinates (Parity operation).
   
   // Redefine the Vector2 worker functions for 3-space
   // (note that Cross has a different return type)
   real_t Mag()      const; //!< @brief The vector's Euclidean magnitude.
   real_t Mag2()     const; //!< @brief The vector's squared Euclidean magnitude.
   
   /*! @brief Pseudorapidity
    * 
    *  \f$ \eta = \text{arctanh}\left(\frac{\tt vec.x3}{\tt vec.T().Mag()}\right) \f$
   */
   real_t Eta() const;
   real_t Theta() const; //! @brief Vec3 = [sin(theta)*Vec2, cos(theta)]

   real_t Dot          (Vector3 const&) const; //!< @brief The dot/inner product.
   Vector3 Cross       (Vector3 const&) const; //!< @brief The cross/outer product (itself a Vector3).
   real_t InteriorAngle(Vector3 const&) const; //!< @brief The cross/outer product.

   //~ // Transverse projection coefficients (Vector2::Mag() / Vector3::Mag())
   //~ // These do not come in non-fma forms
   //~ real_t Tprojection()  const;
   //~ real_t Tprojection2() const;
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
enum class Vec4from4 {EnergyEtaPhiM}; // Add more options as needed

//! @brief When a Vector4 is constructed from a Vector3 and one d.o.f., this enum identifies those d.o.f.
enum class Vec4from2 {Energy, Time, Mass, Length, 
	// Take given p3 and calculate energy from it's boost
	Boost_preserve_p3, 
	BoostMinusOne_preserve_p3,
	// Take given p3 and use it's length as the energy, downscale p3 by beta.
	Boost_preserve_E,
	BoostMinusOne_preserve_E};

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
   using Vector2<real_t>::Phi;
   using Vector3<real_t>::Eta;
   using Vector3<real_t>::Theta;
   real_t x0; //!< @brief The energy/time.
   
   Vector4(); //!< @brief Construct a null vector (0, 0, 0, 0)
   explicit Vector4(bool); //!< @brief Don't initialize the vector (regardless of bool's value)
   explicit Vector4(real_t const x0, real_t const x1, real_t const x2, real_t const x3);
   //! @brief Construct from four generic d.o.f., identified via \p argFormat
   explicit Vector4(real_t const, real_t const, real_t const, real_t const, Vec4from4 const argFormat);
   //! @brief Construct from a spatial 3-vector and one additional d.o.f., identified via \p argFormat
   explicit Vector4(real_t const w0, Vector3<real_t> const& xORp, Vec4from2 const w0type);
   explicit Vector4(Vector3<real_t> const& x); //! @brief Construct a light-like 4-vector (e.g. zero interval).
   virtual ~Vector4() {}

   //! @brief Convert the internal real_t floating-point type
   template<typename convertTo>
   operator Vector4<convertTo>() const;

   Vector4& operator+=(const Vector4&);
   Vector4& operator-=(const Vector4&);
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

   Vector4 operator-() const; // reverse all signs
     
   //! @brief Square root of Length2() (return negative length Length2() is negative, instead of imaginary).
   real_t Length()  const; 
   //! @brief \f$ p^\mu p_\mu \f$. \note see SetLengthRelDiffThreshold
   real_t Length2() const;
   
   real_t Contract(Vector4 const&) const; //!< @brief \f$ p^\mu q_\mu \f$
   
   /*! @brief Rapidity
    * 
    *  \f$ y = \text{arctanh}\left(\frac{\tt vec.x3}{\tt vec.x0}\right) \f$
   */
   real_t Rapidity() const;
   
   static real_t BetaFrom_Gamma(real_t const gamma);
   static real_t BetaFrom_GammaMinusOne(real_t const gm1);
   static real_t BetaFrom_Mass_pSquared(real_t const mass, real_t const pSquared);

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

//~ template <typename real_t>
//~ Vector4<real_t> MasslessVec4_EnergyEtaPhi(real_t const E, real_t const eta, real_t const phi);

// Create nicknames for common real_t types
using Vec4 = Vector4<double>;
using Vec4_f = Vector4<float>;

/*! @brief A boost matrix (lambda) in an arbitrary direction.
 * 
 *  Store lambda - 1 to create the shift, which is added to the original vector.
*/
template <typename real_t>
class LorentzBoost
{
	private:
		// The columns/rows of the boost matrix (it's symmetric, so it doesn't matter which).
		// The symmetry also creates a nearly 2-time redundancy in these coefficients,
		// but its faster/easier to store the redundant coefficients.
		Vector4<real_t> lambda_0, lambda_1, lambda_2, lambda_3;
		
	public:
		/*! @brief Constuct the boost matrix.
		 * 
		 *  If \p CMtoTarget is true, return a matrix which takes the CM frame to the target.
		 *  When false, take the target to the CM frame.
		 */
		LorentzBoost(Vector4<real_t> const& target, bool const CMtoTarget = true);
		
		//! @brief Boost the \p victim.
		Vector4<real_t> operator()(Vector4<real_t> const& victim) const;
};

} // End namespace

#endif
