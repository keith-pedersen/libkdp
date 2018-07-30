from kdp import *

a = Vec4(1., 0., 0., 1.) + Vec4(1., 0., 0., -0.5)
boost = Boost(a.BetaVec())
print(a)
print(boost.BetaVec())
b = boost.Backward(a)
print(a)
