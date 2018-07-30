#include "kdpVectors.hpp"
#include "pqRand/pqRand.hpp"
#include <QtCore/QSettings>

std::mt19937_64 gen;
std::uniform_real_distribution<double> uniform;

double U_0_1()
{
	return uniform(gen);
}

double U_Neg1_1()
{
	return 1. - 2.*U_0_1();
}


kdp::Vec3 IsoVec3()
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
		iso.x1 = U_Neg1_1();
		iso.x2 = U_Neg1_1();
		iso.x3 = U_Neg1_1();
		
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
	
	return iso;
}

int main()
{
	QSettings settings("test_Boost.conf", QSettings::NativeFormat);
	
	size_t const sample_boosts = size_t(settings.value("n_boosts", 1e4).toDouble());
	double const thresh = std::fabs(settings.value("error_thresh", 8e-16).toDouble());
	
	kdp::Vec4 const CM(1., 0., 0., 0.);
		
	for(size_t i = 0; i < sample_boosts; ++i)
	{
		kdp::LorentzBoost<double> boost(IsoVec3());
		kdp::LorentzBoost<double> boost_long(kdp::Vec3(0., 0., U_0_1()));
		
		kdp::Vec4 const photon(IsoVec3());
		kdp::Vec4 const photonTwo(IsoVec3());		
		kdp::Vec4 const massive(U_0_1(), IsoVec3(), kdp::Vec4from2::Mass);
		
		kdp::Vec4 const photon_boosted = boost.Forward(photon);
		kdp::Vec4 const photon_boosted_long = boost_long.Forward(photon);
		kdp::Vec4 const photon_restored = boost.Backward(photon_boosted);
		
		kdp::Vec4 const CM_boosted = boost.Forward(CM);
		kdp::Vec4 const CM_restored = boost.Backward(CM_boosted);
		
		double const deltaR_photons_before = photon.DeltaR_pseudo(photonTwo);
		double const deltaR_photons_after = photon_boosted_long.DeltaR_pseudo(boost_long.Forward(photonTwo));
		
		kdp::Vec4 const CM_boosted_long = boost_long.Forward(CM);
		kdp::Vec4 const CM_restored_long = boost_long.Backward(CM_boosted_long);
		
		kdp::Vec4 const massive_boosted_long = boost_long.Forward(massive);
		kdp::Vec4 const massive_restored_long = boost_long.Backward(massive_boosted_long);
		
		double const deltaR_massive_before = massive_boosted_long.DeltaR_rap(CM_boosted_long);
		double const deltaR_massive_after = massive_restored_long.DeltaR_rap(CM_restored_long);
		
		double const betaError_boosted = kdp::RelError(CM_boosted.Beta(), boost.Beta());
		double const massError_boosted = kdp::RelError(CM_boosted.Mass(), 1.);
		double const betaError_restored = CM_restored.Beta();
		double const massError_restored = kdp::RelError(CM_restored.Mass(), 1.);
		
		double const betaError_photon = 1. - photon_boosted.Beta();
		double const massError_photon = photon_boosted.Mass();
		
		double const deltaR_photon_error = kdp::RelDiff(deltaR_photons_before, deltaR_photons_after);
		double const deltaR_massive_error = kdp::RelDiff(deltaR_massive_before, deltaR_massive_after);
		
		if(betaError_boosted > thresh)
			printf("Beta of boosted CM misses target (rel error): %.3e\n", betaError_boosted);
		if(massError_boosted > thresh)
			printf("Mass of boosted CM misses target (rel error): %.3e\n", massError_boosted);
		
		if(betaError_restored > thresh)
			printf("Beta of restored CM misses target: %.3e\n", betaError_restored);
		if(massError_restored > thresh)
			printf("Mass of restored CM misses target: %.3e\n", massError_restored);
			
		if(deltaR_photon_error > thresh)
			printf("deltaR_pseudo not boost invariant for photons (rel error): %.3e\n", deltaR_photon_error);
		if(deltaR_massive_error > thresh)
			printf("deltaR_rap not boost invariant for massive (rel error): %.3e\n", deltaR_massive_error);
			
		if(betaError_photon > thresh)
			printf("(1 - Beta) of boosted photon misses target: %.3e\n", betaError_photon);
		if(massError_photon > thresh)
			printf("Mass of boosted photon misses target: %.3e\n", massError_photon);
		
	}
		
	return 0;
}
