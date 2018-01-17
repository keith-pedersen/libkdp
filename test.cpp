#include "kdpVectors.hpp"
#include "pqRand/pqRand.hpp"
// This program takes about 1.25 seconds to complete
#include <cmath>

using real_t = pqRand::real_t;

kdp::Vec4 IsotropicMasslessVec4(pqRand::engine& gen)
{
	real_t length2;
	kdp::Vec3 p3;
	//~ kdp::Vec3 p3(false); // Python just can't do this, so it's not a fair comparison
	
	do
	{
		p3.x1 = gen.U_uneven();
		p3.x2 = gen.U_uneven();
		p3.x3 = gen.U_uneven();
		length2 = p3.Mag2();
	}while(length2 >= real_t(1));
	
	gen.ApplyRandomSign(p3.x1);
	gen.ApplyRandomSign(p3.x2);
	gen.ApplyRandomSign(p3.x3);
			
	return kdp::Vec4(real_t(0), p3, kdp::Vec4from2::Mass);
}

int main()
{
	pqRand::engine gen;
	
	size_t const sampleSize = 5e6;
	std::vector<kdp::Vec4> sample;
	//~ sample.reserve(sampleSize);
	
	for(size_t i = 0; i < sampleSize; ++i)
		sample.push_back(IsotropicMasslessVec4(gen));
		
	printf("vecs drawn\n");
	
	kdp::Vec4 total;
		
	for(kdp::Vec4 const& vec : sample)
		total += vec;
		
	total.p() = -total.p();
	sample.push_back(total);
	
	printf("%.3e %.3e %.3e %.3e\n", 
		sample.back().x0,
		sample.back().x1,
		sample.back().x2,
		sample.back().x3);
}
