#include "myEkf.h"
#include <fstream>
#include <string>
int main()
{
	myEkf filter;

	Matrix<double,10,10> P0;
	for (size_t i = 1; i < 11; i++)
	{
		for (size_t j = 1; j < 11; j++)
			if (i == j)
				P0(i, j) = 0.01; 
	}

	std::ifstream in("Measurement2.txt");
	std::string Measure;
	std::string::size_type sz;     // alias of size_t
	Vector z(6);
	Vector u(3);
	while (std::getline(in, Measure))
	{
		z(1) = std::stod(Measure, &sz);
		Measure = Measure.substr(sz);
		z(2) = std::stod(Measure, &sz);
		Measure = Measure.substr(sz);
		z(3) = std::stod(Measure, &sz);
		Measure = Measure.substr(sz);
		z(4) = std::stod(Measure, &sz);
		Measure = Measure.substr(sz);
		z(5) = std::stod(Measure, &sz);
		Measure = Measure.substr(sz);
		z(6) = std::stod(Measure, &sz);
		Measure = Measure.substr(sz);
		u(1) = std::stod(Measure, &sz);
		Measure = Measure.substr(sz);
		u(2) = std::stod(Measure, &sz);
		Measure = Measure.substr(sz);
		u(3) = std::stod(Measure, &sz);
		filter.step(u, z);
		std::cout << filter.getX() << std::endl;
	}
}