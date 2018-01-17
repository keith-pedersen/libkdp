include "kdp/kdpVectors.hpy"

from libcpp.memory cimport shared_ptr # First import shared_ptr

cdef class Vec2:
	# To allow the equivalent of C++ referencess (e.g. the Vec2 Vec3.T()), 
	# we must use a shared_ptr
	cdef shared_ptr[Vec2_c] vec
	
	@staticmethod
	cdef Vec2 Factory(const Vec2_c& orig)
		
cdef class Vec3:
	cdef shared_ptr[Vec3_c] vec
	
	@staticmethod
	cdef Vec3 Factory(const Vec3_c& orig)

cdef class Vec4:
	cdef shared_ptr[Vec4_c] vec
	
	@staticmethod
	cdef Vec4 Factory(const Vec4_c& orig)
