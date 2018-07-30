// Copyright (C) 2018 by Keith Pedersen (Keith.David.Pedersen@gmail.com)

// Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
// The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

// This file tests the 1D and 2D histograms from kdpHistogram.hpp

#include <random>

#include "kdpHistogram.hpp"
#include "kdpVectors.hpp"

std::mt19937_64 gen;
std::uniform_real_distribution<double> uniform;

using kdp::Vec3;

// Draw an isotropic vector using rejection sampling, but don't normalize the length
Vec3 IsotropicVector()
{
	Vec3 vec(false);
	double r2;
	
	do
	{
		vec.x1 = 1. - 2.*uniform(gen);
		vec.x2 = 1. - 2.*uniform(gen);
		vec.x3 = 1. - 2.*uniform(gen);
		
		r2 = vec.Mag2();
	}
	while(r2 > 1.);
	
	return vec;
}


int main()
{
	kdp::BinSpecs U01_linear("U01", 100, {0., 1.}, kdp::UniformBinType::UNIFORM);
	kdp::BinSpecs U01_less("U01", 20, {0., 1.}, kdp::UniformBinType::UNIFORM);
	kdp::BinSpecs U01_log("U01", 100, {1e-8, 1.}, kdp::UniformBinType::LOG_UNIFORM);
	
	// Bin the square magnitude of the isotropic vectors
	kdp::Histogram1 r2_linear("r2_linear", U01_linear);
	kdp::Histogram1 r2_log("r2_log", U01_log);
	
	// We draw two vectors and bin their square magnitude (r2_A, r2_B) into a 2D histogram
	kdp::Histogram2 pairs("pairs", U01_less, U01_less);
	// We will add the vectors and "accept" them if their sum still lies in the unit sphere.
	kdp::Histogram2 accepted("accepted", U01_less, U01_less);
	
	for(size_t i = 0; i < 1000000; ++i)
	{
		Vec3 vecA = IsotropicVector();
		Vec3 vecB = IsotropicVector();
		
		double const r2_A = vecA.Mag2();
		double const r2_B= vecB.Mag2();
		
		r2_linear.Fill(r2_A);
		r2_log.Fill(r2_A);
		
		r2_linear.Fill(r2_B);
		r2_log.Fill(r2_B);
		
		pairs.Fill(r2_A, r2_B);
		
		
		if((vecA + vecB).Mag2() < 1.)
			accepted.Fill(r2_A, r2_B);
	}
	
	// Dividing accepted by pairs creates a new histogram that shows 
	// the probability of acceptance as a function of the two square magnitudes
	accepted.WriteRatio(pairs, "pAccept");
	
	return 0;
}
