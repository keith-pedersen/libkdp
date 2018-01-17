# Testing reveals that the python pqRand and kdp libraries take 15 s to do this task.
# The C++ code takes 1.25 seconds.
# Thus, Cython is about 12 times slower than the equivalent C++ code. 
# All the Python time is spent in IsotropicMasslessVec4().
# Using O2 versus O3 when building the libraries doesn't seem to help much.
# Changing various Cython compilation options has no effect.
# Building statically against kdp.o versus kdp.so had no effect.
# One third of the time is spent initializing a 4-vector without
# checking it's length and doing the sign flip.
# The length check rejects 1/2 of all samples (4 pi / 3) / 8.
# So, in my estimation, 10 times slower is as fast as it's going to get, 
# which is really not that bad for an interpretive language.

from kdp import *
import pYqRand as pqr
import math

def IsotropicMasslessVec4(engine):
	length2 = 2.
	while(length2 >= 1.):
		p3 = Vec3(gen.U_uneven(), gen.U_uneven(), gen.U_uneven())
		length2 = p3.Mag2()
	
	p3.x1 = gen.ApplyRandomSign(p3.x1)
	p3.x2 = gen.ApplyRandomSign(p3.x2)
	p3.x3 = gen.ApplyRandomSign(p3.x3)
		
	return Vec4.From_MassP3(0., p3)
	
v3 = Vec3()
v2 = v3.T()
v2.x1 = 10.
v2.x2 = 10.
v2 /= 10 
print(v3) # test that v2 is bound to v3
v2 = "not a vector" # test that v2 releases its bindings
print(v2)

#Make a bunch of massless jets, then balance momentum
gen = pqr.engine()

#~ vecs = [Vec4(0., gen.U_uneven(), gen.U_uneven(), gen.U_uneven()) for __ in range(int(5e6))] # 5 s 
vecs = [IsotropicMasslessVec4(gen) for __ in range(int(5e6))] # 15 s 
#~ vecs = [gen.U_uneven() for __ in range(int(5e6))]

print("vecs drawn") # it's the drawing that takes forever

total = Vec4()
for vec in vecs:
	total += vec

total.p().Assign(-total.p())
vecs += [total,]

# print the 4-vectors with a new line after each one
#~ print("\n".join(list(map(str, vecs))))
print(vecs[-1])
