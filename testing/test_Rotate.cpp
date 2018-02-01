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
		auto const omega = gen.ApplyRandomSign(M_PI*gen.U_uneven());
		
		// The unambigous rotation from u to v (with a post-rotation about v)
		kdp::Rotate3<double> rotator(u, v, omega);
		kdp::Rotate3<double> identity(u, u*gen.U_uneven(), 0.);
		
		bool errorCaught = false;
		try
		{
			kdp::Rotate3<double> flip(u, -(u*gen.U_uneven()), 0.);
		}
		catch(std::invalid_argument e)
		{
			errorCaught = true;
		}
		
		if(not errorCaught)
			printf("Rotate3 did not throw exception on anti-parallel vectors!\n");
		
		auto const axis = rotator.Axis();
		
		// To test the unambigous rotation, we use two sequential rotations
		// defined by axis/angle
		kdp::Rotate3<double> rotator_1(u.Cross(v).Normalize(), u.InteriorAngle(v));
		kdp::Rotate3<double> rotator_2(v, omega);
	
		double const axisLengthError = kdp::AbsRelError(axis.Mag(), double(1));
		double const uv_same = std::fabs(u.Dot(axis)/u.Mag() - v.Dot(axis)/v.Mag());
	
		if((axisLengthError > thresh) or (uv_same > thresh))
		{
			printf("\n~~~~~~~~~~~~~~~~~~~\n");
			printf("%lu \n", i);
			printf("axis is normalized ... rel error : %.3e\n", axisLengthError);
			printf("axis is same dist  ... abs dot: %.3e %.3e\n", uv_same, u.Dot(axis)/u.Mag());
			printf("~~~~~~~~~~~~~~~~~~~\n");
		}
		
		std::vector<kdp::Vec3> origVec;
		std::vector<kdp::Vec3> rotVec;
		std::vector<kdp::Vec3> rotVec_v2;
		
		for(size_t j = 0; j < sample_vecs; ++j)
		{
			kdp::Vec3 const orig = IsoVec3(gen);
			origVec.push_back(orig);
			
			rotVec.push_back(rotator(orig));
			rotVec_v2.push_back(rotator_2(rotator_1(orig)));
			
			double const identifyError = kdp::AbsRelError(identity(orig).Mag(), orig.Mag());
			if(identifyError > 0.)
				printf("identity is not the identity: %.3e\n", identifyError);
			
			double const angle = std::fabs(rotVec.back().Dot(axis) - orig.Dot(axis))/(orig.Mag());
			double const angle_2 = std::fabs(rotVec_v2.back().Dot(axis) - orig.Dot(axis))/(orig.Mag());
			double const mag  = kdp::AbsRelError(rotVec.back().Mag(), orig.Mag());
			
			auto const orig_perp = orig - rotator_2.Axis() * orig.Dot(rotator_2.Axis());
			auto const rotPart = rotator_2(orig);
			auto const final_perp = rotPart - rotator_2.Axis() * rotPart.Dot(rotator_2.Axis());
			double const angleError = std::fabs(orig_perp.InteriorAngle(final_perp) - std::fabs(rotator_2.Angle()));
			
			double const error = kdp::AbsRelError(rotVec.back().Dot(rotVec_v2.back()), rotVec.back().Mag2());
			//~ double const error = rotVec.back().InteriorAngle(rotVec_v2.back());
			
			if((angle > thresh) or (mag > thresh))
				printf("%lu -- angle, length preserved ... rel error: %.3e, %.3e   ..... %.3e, %.3e\n", 
					i, angle, mag, u.Dot(v)/std::sqrt(u.Mag2()*v.Mag2()), orig.Dot(axis));
					
			if(error > thresh)
				printf("%lu %lu -- error between methods (in radians): %.3e %.3e\n", i, j, error, angle_2);
				
			if(angleError > thresh)
				printf("%lu %lu -- rotation angle does not match target: %.3e, %.3e, %.3e\n", i, j, angleError, orig.Dot(rotator_2.Axis()), rotator_2.Angle());
		}
		
		for(size_t k = 0; k < sample_vecs; ++k)
		{
			for(size_t l = k + 1; l < sample_vecs; ++l)
			{
				double const dotError = std::fabs(rotVec[k].Dot(rotVec[l]) - origVec[k].Dot(origVec[l]));
				
				if(dotError > thresh)
					printf("%lu %lu %lu interior-dot fail: %.3e\n", i, k, l, dotError);
			}
		}
	}
}
