#include <vector>
#include "kdpStdVectorMath.hpp"

int main()
{
	std::vector<double> a(10, 1.);
	std::vector<double> b(9, -1.);
	
	a += b;
	a -= b;
	
	printf("%.3e\n", a.front());
	
	return 0;
}
