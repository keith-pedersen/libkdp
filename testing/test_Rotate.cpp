#include "kdpVectors.hpp"
#include "pqRand/pqRand.hpp"
#include <QtCore/QSettings>

kdp::Vec3 IsoVec3(pqRand::engine& gen)
{
	// The easiest way to draw an isotropic 3-vector is rejection sampling.
	// It provides slightly better precision than a theta-phi implementation,
	// and is not much slower.
	
	kdp::Vec3 iso(false); // false means don't initialize the Vec3
	double r2;
	
	do
	{
		// U_S has enough precision for this application
		// because we scale by r2
		iso.x1 = gen.U_even();
		iso.x2 = gen.U_even();
		iso.x3 = gen.U_even();
		
		r2 = iso.Mag2();
	}
	while(r2 > 1.); // Draw only from the unit sphere (in the positive octant)
	
	// It we wanted to scale the length from some distribution, 
	// the current length r**2 is a random variate (an extra DOF since 
	// only two DOF are required for isotropy) which can be 
	// easily converted into a random U(0, 1) via its CDF 
	// 	u = CDF(r2) = (r2)**(3/2)
	// (you can work this out from the differential volume dV/(V * d(r**2))).
	// This U(0,1) can then be plugged into any quantile function Q(u)
	// for the new length. One should check to see if 
	// dividing out by the current length can be rolled into Q(u), 
	// so that you that you use u to draw the scale that the takes 
	// iso to its new length.
	
	// In this application, we simply want unit vectors, 
	// so we throw away the random entropy of r2.
	//~ iso /= std::sqrt(r2);
		
	// Use random signs to move to the other 7 octants
	gen.ApplyRandomSign(iso.x1);
	gen.ApplyRandomSign(iso.x2);
	gen.ApplyRandomSign(iso.x3);
	
	return iso;
}

int main()
{
	pqRand::engine gen;
	
	QSettings settings("test_Rotate.conf", QSettings::NativeFormat);
	
	size_t const sample_rotations = size_t(settings.value("n_rotations", 1e4).toDouble());
	size_t const sample_vecs = size_t(settings.value("n_vecs", 1e4).toDouble());
	double const thresh = std::fabs(settings.value("error_thresh", 8e-16).toDouble());
	
	for(size_t i = 0; i < sample_rotations; ++i)
	{
		auto const u = IsoVec3(gen);
		auto const v = IsoVec3(gen);
		
		kdp::Rotate3<double> rotator(u, v);
		auto const axis = rotator.Axis();
	
		double const axisLengthError = kdp::AbsRelError(axis.Mag(), double(1));
		double const u_ortho = std::fabs(u.Dot(axis)/u.Mag());
		double const v_ortho = std::fabs(v.Dot(axis)/v.Mag());
	
		if((axisLengthError > thresh) or (u_ortho > thresh) or (v_ortho > thresh))
		{
			printf("\n~~~~~~~~~~~~~~~~~~~\n");
			printf("%lu \n", i);
			printf("axis is normalized ... rel error : %.3e\n", axisLengthError);
			printf("axis is orthogonal ... abs dot: %.3e, %.3e\n", u_ortho, v_ortho);
			printf("~~~~~~~~~~~~~~~~~~~\n");
		}
		
		for(size_t j = 0; j < sample_vecs; ++j)
		{
			kdp::Vec3 const orig = IsoVec3(gen);
			auto const rotated = rotator(orig);
			
			double const angle = std::fabs(rotated.Dot(axis) - orig.Dot(axis))/(orig.Mag());
			double const mag  = kdp::AbsRelError(rotated.Mag(), orig.Mag());
			
			if((angle > thresh) or (mag > thresh))
				printf("%lu -- angle, length preserved ... rel error: %.3e, %.3e   ..... %.3e, %.3e\n", 
					i, angle, mag, u.Dot(v)/std::sqrt(u.Mag2()*v.Mag2()), orig.Dot(axis));
		}
	}
}
