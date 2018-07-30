#cython: language_level=3, optimize.use_switch=False

include "kdpVectors.hpy"

# This file creates wrappers of Vec2, Vec3, Vec4, Rot3, and Boost from "kdpVectors.hpp".
# These wrappers behave identically to the C++ versions (except they cannot be const, 
# and can only be constructed from within Python using the Cartesian constructor).
# An important property is the ability to access sub-vectors with a bound reference. 
# For example, the following code gets and alters the transverse components of a Vec3
#
#		x3 = Vec3(1, 1, 1)
#		print(x3)
#		x2 = x3.T()
#		x2 *= 5
#		print(x3)
#
# Getting this to work was quite a pain in the ass.
# I recount this story for the benefit of myself, in the future, 
# writing Cython code and having forgot all these lessons.

# First, I tried using shared_ptr, so that multiple 
# Python objects are referencing the same underlying object.
# The first time I did it wrong
#
# 		vec2_in_vec3 = shared_ptr[Vec2_c](<Vec2_c>vec3_shared_ptr.get())
# 
# The problem, of course, is that when we intialize from a non-shared pointed
# (calling shared_ptr.get() then casting it to a new type), 
# the shared_ptr count starts at 1, when we definitely need it to start at 2.
# So the Vec2 can delete the Vec3. Not good!

# The most straightforward way to fix this problem is to use static_pointer_cast()
# to cast a shared_ptr<Vec3_c> to a shared_ptr <Vec2_c>, which increments the count. 
# However, this function is not exposed by Cython's libcpp.memory.
# Eventually, I figured out that libcpp.memory is just a declaration file which 
# exposes the C++ functions (very much like my "kdpVectors.hpy").
# So I figured out how to expose static_pointer_cast() to Cython (see below).
# This almost worked.

# Internally, static_pointer_cast uses static_cast, 
# which fails because "Vec2 is an inaccessible base of Vec3."
# Damn; this is my architectural conundrum rearing its ugly head
#		A Vec3 "has a" Vec2, but the "composition" interface is cumbersome, 
# 		so I used private/protected inheritance to abuse OOP.
# Making the inheritence public fixes the problem at hand, 
# but absolutely breaks all the protections I desired by using private inheritance 
# (e.g. Vec2.Dot(Vec3) is suddenly/silently allowed).

# So I'm back to square one, until I stumble upon 
# shared_ptr's aliasing constructor
#
#		template <class U> shared_ptr (const shared_ptr<U>& x, element_type* p) noexcept;
#
#		"The object does not own p, and will not manage its storage. 
#		Instead, it co-owns x's managed object and counts as one additional use of x. 
#		It will also delete x's pointer on release (and not p).
#		It can be used to point to members of objects that are already managed."
# 
# This is EXACTLY what I want, so I try to manually expose shared_ptr
# (instead of using "from libcpp.memory cimport shared_ptr")	
# 
#		cdef extern from "<memory>" namespace "std" nogil:
#			cdef cppclass shared_ptr[T]:
#				shared_ptr()
#				shared_ptr(T*)
#				shared_ptr(shared_ptr[T]&)
#				shared_ptr[U](shared_ptr[U]&, T*) # < But Cython chokes on this
#
#				T* get()
#				T& operator*()
#		
# Unfortunately, Cython chokes on the constructor I need, the templated ctor.
# Curiously, in libcpp.memory, shared_ptr has a templated function which 
# Cython does not choke on ...
#				...
#				bool owner_before[Y](const shared_ptr[Y]&)
#				...
# This gave me an idea; I can wrap the constructor I need in 
# a function I write (shared_ptr_alias()) and define in "shared_ptr_alias.hpp". 
# I can then expose this function to Cython as a template function 
# (versus a template ctor, which it seems to choke on). 
# By the beard of Zeus it works!

from libcpp.memory cimport shared_ptr # First import shared_ptr

# Now import shared_ptr's aliasing constructor via a standalone funcxtion
cdef extern from "shared_ptr_alias.hpp" nogil:
	cdef shared_ptr[T] shared_ptr_alias[T, U](const shared_ptr[U]&, T* other)
	
# This line demonstrates how to expose functionality directly from the standard library.
# However, we do not need static_pointer_cast anymore. 
cdef extern from "<memory>" namespace "std" nogil:
	cdef shared_ptr[T] static_pointer_cast[T, U](shared_ptr[U]&)
	cdef shared_ptr[T] make_shared[T](...)
		
from libcpp cimport bool
from cython.operator cimport dereference as deref 

cdef class Vec2:
	# To allow the equivalent of C++ referencess (e.g. the Vec2 Vec3.T()), 
	# we must use a shared_ptr
#~ 	cdef shared_ptr[Vec2_c] vec
	
	def __cinit__(self, double x1 = 0., double x2 = 0.):
#~ 		self.vec = shared_ptr[Vec2_c](new Vec2_c(x1, x2))
		# It appears that using make_shared is *slightly* faster, but creates a 10% larger binary
		self.vec = make_shared[Vec2_c](x1, x2)
		
	def __dealloc__(self):
		return #shared_ptr knows what to do
		
	def Assign(self, Vec2 other):
		deref(self.vec).assign(deref(other.vec))
		return self
		
	@staticmethod
	cdef Vec2 Factory(const Vec2_c& orig):
		copy = Vec2()
		deref(copy.vec).assign(orig)
		return copy
		
	def Copy(self):
		return Vec2.Factory(deref(self.vec))
		
	def __deepcopy__(self, memo):
		return self.Copy()
		
	def __repr__(self):
		return str([self.x1, self.x2])
	
	# This way of defining properties is "deprecated", 
	# but with the "pythonic" way (@property, etc.), the setter doesn't work.
	property x1:
		# "A doc string can go here."
		def __get__(self):
			return deref(self.vec).x1

		def __set__(self, value):
			deref(self.vec).x1 = value
			
	property x2:
		def __get__(self):
			return deref(self.vec).x2

		def __set__(self, value):
			deref(self.vec).x2 = value

# The broken way to define properties
#~ 	@property
#~ 	def x1(self):
#~ 		return deref(self.vec).x1

#~ 	@x1.setter
#~ 	def x1(self, double value):
#~ 		deref(self.vec).x1 = value	
		
	def __iadd__(self, Vec2 other):
		deref(self.vec).addEquals(deref(other.vec))
		return self
		
	def __isub__(self, Vec2 other):
		deref(self.vec).subEquals(deref(other.vec))
		return self
		
	def __imul__(self, double scale):
		deref(self.vec).mulEquals(scale)
		return self
		
	def __itruediv__(self, double scale):
		deref(self.vec).divEquals(scale)
		return self
		
	def __add__(Vec2 this, Vec2 that):
		return this.Copy().__iadd__(that)
		
	def __sub__(Vec2 this, Vec2 that):
		return this.Copy().__isub__(that)
	
	# Currently, __rmul__ is not supported
	def __mul__(Vec2 this, double scale):
		return this.Copy().__imul__(scale)
	
	# Note this is truediv, not div, which is what Python 3 will call for  a / b
	def __truediv__(Vec2 this, double scale):
		return this.Copy().__itruediv__(scale)

	def __neg__(self):
		ret = Vec2()
		deref(ret.vec).assign(-deref(self.vec))
		return ret
		
	def Normalize(self):
		deref(self.vec).Normalize()
		return self
		
	def Mag(self):
		return deref(self.vec).Mag()
		
	def Mag2(self):
		return deref(self.vec).Mag2()
		
	def Phi(self):
		return deref(self.vec).Phi()
		
	def DeltaPhi(self, Vec2 other):
		return deref(self.vec).DeltaPhi(deref(other.vec))
		
	def Dot(self, Vec2 other):
		return deref(self.vec).Dot(deref(other.vec))
		
	def Cross(self, Vec2 other):
		return deref(self.vec).Cross(deref(other.vec))
		
	def InteriorAngle(self, Vec2 other):
		return deref(self.vec).InteriorAngle(deref(other.vec))

########################################################################

cdef class Vec3:
#~ 	cdef shared_ptr[Vec3_c] vec
	
	def __cinit__(self, double x1 = 0., double x2 = 0., double x3 = 0.):
#~ 		self.vec = shared_ptr[Vec3_c](new Vec3_c(x1, x2, x3))
		self.vec = make_shared[Vec3_c](x1, x2, x3)
	
	def __dealloc__(self):
		return
		
	def Assign(self, Vec3 other):
		deref(self.vec).assign(deref(other.vec))
		
	@staticmethod
	cdef Vec3 Factory(const Vec3_c& orig):
		copy = Vec3()
		deref(copy.vec).assign(orig)
		return copy
		
	@staticmethod
	def Make_LengthThetaPhi(double length, double theta, double phi):
		return Vec3.Factory(Vec3_c(length, theta, phi, V3f_LengthThetaPhi))
		
	def Copy(self):
		return Vec3.Factory(deref(self.vec))
		
	def __deepcopy__(self, memo):
		return self.Copy()
				
	def __repr__(self):
		return str([self.x1, self.x2, self.x3])
	
	property x1:
		def __get__(self):
			return deref(self.vec).x1

		def __set__(self, value):
			deref(self.vec).x1 = value
			
	property x2:
		def __get__(self):
			return deref(self.vec).x2

		def __set__(self, value):
			deref(self.vec).x2 = value
			
	property x3:
		def __get__(self):
			return deref(self.vec).x3

		def __set__(self, value):
			deref(self.vec).x3 = value
		
	def T(self):
		vec2 = Vec2()
		vec2.vec = shared_ptr_alias[Vec2_c, Vec3_c](self.vec, <Vec2_c*>(self.vec.get()))
		return vec2
	
	def __iadd__(self, Vec3 other):
		deref(self.vec).addEquals(deref(other.vec))
		return self
		
	def __isub__(self, Vec3 other):
		deref(self.vec).subEquals(deref(other.vec))
		return self
		
	def __imul__(self, double scale):
		deref(self.vec).mulEquals(scale)
		return self
		
	def __itruediv__(self, double scale):
		deref(self.vec).divEquals(scale)
		return self
		
	def __add__(Vec3 this, Vec3 that):
		return this.Copy().__iadd__(that)
		
	def __sub__(Vec3 this, Vec3 that):
		return this.Copy().__isub__(that)
	
	# Currently, __rmul__ is not supported
	def __mul__(Vec3 this, double scale):
		return this.Copy().__imul__(scale)
	
	# Note this is truediv, not div
	def __truediv__(Vec3 this, double scale):
		return this.Copy().__itruediv__(scale)

	def __neg__(self):
		ret = Vec3()
		deref(ret.vec).assign(-deref(self.vec))
		return ret
		
	def Normalize(self):
		deref(self.vec).Normalize()
		return self
		
	def Mag(self):
		return deref(self.vec).Mag()
		
	def Mag2(self):
		return deref(self.vec).Mag2()
		
	def Phi(self):
		return deref(self.vec).Phi()
		
	def Eta(self):
		return deref(self.vec).Eta()
		
	def Theta(self):
		return deref(self.vec).Theta()
		
	def DeltaEta(self, Vec3 other):
		return deref(self.vec).DeltaEta(deref(other.vec))
		
	def DeltaPhi(self, Vec3 other):
		return deref(self.vec).DeltaPhi(deref(other.vec))
		
	def Dot(self, Vec3 other):
		return deref(self.vec).Dot(deref(other.vec))
		
	def Cross(self, Vec3 other):
		ret = Vec3()
		deref(ret.vec).assign(deref(self.vec).Cross(deref(other.vec)))
		return ret
		
	def InteriorAngle(self, Vec3 other):
		return deref(self.vec).InteriorAngle(deref(other.vec))

########################################################################

cdef class Vec4:
#~ 	cdef shared_ptr[Vec4_c] vec
	
	def __cinit__(self, double x0 = 0., double x1 = 0., double x2 = 0., double x3 = 0.):
#~ 		self.vec = shared_ptr[Vec4_c](new Vec4_c(x0, x1, x2, x3))
		self.vec = make_shared[Vec4_c](x0, x1, x2, x3)
		
	@staticmethod	
	def From_MassP3(double mass, Vec3 p3):
		ret = Vec4()
#~ 		ret.vec = shared_ptr[Vec4_c](new Vec4_c(mass, deref(p3.vec), Mass))
		ret.vec = make_shared[Vec4_c](mass, deref(p3.vec), V4f2_Mass)
		return ret
		
	def __dealloc__(self):
		return
		
	@staticmethod
	cdef Vec4 Factory(const Vec4_c& orig):
		copy = Vec4()
		deref(copy.vec).assign(orig)
		return copy
		
	def Copy(self):
		return Vec4.Factory(deref(self.vec))
		
	def __deepcopy__(self, memo):
		return self.Copy()
		
	def __repr__(self):
		return str([self.x0, [self.x1, self.x2, self.x3]])
	
	property x0:
		# "A doc string can go here."
		def __get__(self):
			return deref(self.vec).x0

		def __set__(self, value):
			deref(self.vec).x0 = value
		
	property x1:
		def __get__(self):
			return deref(self.vec).x1

		def __set__(self, value):
			deref(self.vec).x1 = value
			
	property x2:
		def __get__(self):
			return deref(self.vec).x2

		def __set__(self, value):
			deref(self.vec).x2 = value
			
	property x3:
		def __get__(self):
			return deref(self.vec).x3

		def __set__(self, value):
			deref(self.vec).x3 = value
		
	def x(self):
		vec3 = Vec3()
		vec3.vec = shared_ptr_alias[Vec3_c, Vec4_c](self.vec, <Vec3_c*>(self.vec.get()))
		return vec3
		
	def p(self):
		return self.x()
		
	def __iadd__(self, Vec4 other):
		deref(self.vec).addEquals(deref(other.vec))
		return self
		
	def __isub__(self, Vec4 other):
		deref(self.vec).subEquals(deref(other.vec))
		return self
		
	def __imul__(self, double scale):
		deref(self.vec).mulEquals(scale)
		return self
		
	def __itruediv__(self, double scale):
		deref(self.vec).divEquals(scale)
		return self
		
	def __add__(Vec4 this, Vec4 that):
		return this.Copy().__iadd__(that)
		
	def __sub__(Vec4 this, Vec4 that):
		return this.Copy().__isub__(that)
	
	# Currently, __rmul__ is not supported
	def __mul__(Vec4 this, double scale):
		return this.Copy().__imul__(scale)
	
	# Note this is truediv, not div
	def __truediv__(Vec4 this, double scale):
		return this.Copy().__itruediv__(scale)

	def __neg__(self):
		ret = Vec4()
		deref(ret.vec).assign(-deref(self.vec))
		return ret
		
	def Length(self):
		return deref(self.vec).Length()
		
	def Length2(self):
		return deref(self.vec).Length2()
	
	def Mass(self):
		return deref(self.vec).Mass()
	
	def Mass2(self):
		return deref(self.vec).Mass2()
		
	def Length_T(self):
		return deref(self.vec).Length_T()
		
	def Length2_T(self):
		return deref(self.vec).Length2_T()
		
	def Contract(self, Vec4 other):
		return deref(self.vec).Contract(deref(other.vec))
		
	def Beta(self):
		return deref(self.vec).Beta()
		
	def BetaVec(self):
		return Vec3.Factory(deref(self.vec).BetaVec())
		
	def Rapidity(self):
		return deref(self.vec).Rapidity()
		
	def Eta(self):
		return deref(self.vec).Eta()
		
	def Theta(self):
		return deref(self.vec).Theta()
		
	def Phi(self):
		return deref(self.vec).Phi()
		
	def DeltaEta(self, Vec4 other):
		return deref(self.vec).DeltaEta(deref(other.vec))
		
	def DeltaPhi(self, Vec4 other):
		return deref(self.vec).DeltaPhi(deref(other.vec))
		
	def DeltaR2_rap(self, Vec4 other):
		return deref(self.vec).DeltaR2_rap(deref(other.vec))
		
	def DeltaR_rap(self, Vec4 other):
		return deref(self.vec).DeltaR_rap(deref(other.vec))
		
	def DeltaR2_pseudo(self, Vec4 other):
		return deref(self.vec).DeltaR2_pseudo(deref(other.vec))
		
	def DeltaR_pseudo(self, Vec4 other):
		return deref(self.vec).DeltaR_pseudo(deref(other.vec))
		
	@staticmethod
	def SetLengthRelDiffThreshold(double newThreshold):
		Vec4.SetLengthRelDiffThreshold(newThreshold)
		
########################################################################
		
cdef class Rot3:
	
	def __cinit__(self, Vec3 axis, double angle):
#~ 		self.vec = shared_ptr[Vec4_c](new Vec4_c(x0, x1, x2, x3))
		self.rot = new Rot3_c(deref(axis.vec), angle)
		
	def __dealloc__(self):
		del self.rot
		
	# Have to make this to get the axis for vec to vec rotation
	@staticmethod
	def Make_Vec2Vec(Vec3 u, Vec3 v, double omega):
		# Create a dummy rotation object to solve for axis/angle
		# use a pointer so we don't have to declare default ctor for Cython
		cdef Rot3_c* tmp = new Rot3_c(deref(u.vec), deref(v.vec), omega)
		rot = Rot3(Vec3.Factory(tmp.Axis()), tmp.Angle())
		
		del tmp # must clean
		return rot
	
	def __call__(self, Vec3 victim):
		self.rot.Rotate(deref(victim.vec))
		return victim
		
	def Axis(self):
		return Vec3.Factory(self.rot.Axis())
		
	def Angle(self):
		return self.rot.Angle()		

########################################################################
	
cdef class Boost:
	
	def __cinit__(self, Vec3 beta):
		self.boost = new Boost_c(deref(beta.vec))
		
	def __dealloc__(self):
		del self.boost
		
	def Axis(self):
		return Vec3.Factory(self.boost.Axis())
		
	def BetaVec(self):
		return Vec3.Factory(self.boost.BetaVec())
		
	def Gamma(self):
		return self.boost.Gamma()
		
	def Beta(self):
		return self.boost.Beta()
		
	def Rapidity(self):
		return self.boost.Rapidity()
		
	def Forward(self, Vec4 victim):
		self.boost.Forward(deref(victim.vec))
		return victim
		
	def Backward(self, Vec4 victim):
		self.boost.Backward(deref(victim.vec))
		return victim	
