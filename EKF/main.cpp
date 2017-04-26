#include "myEkf.h"
#include <fstream>
#include <string>
int main()
{
	myEkf filter;

	VectorXd P0Vector(10);
	P0Vector << 0.0, 0.0, 0.0, 0.0,1e-4,1e-4,1e-4,3e-8,3e-8,3e-8;
	Matrix<double,10,10> P0(P0Vector.asDiagonal());
	Matrix<double, 10, 1> x0 = Matrix<double, 10, 1>::Zero();
	x0(0) = 1;
	filter.initialize(x0, P0);
	std::ifstream in("Measurement.txt");
	std::ofstream out("out.txt");
	std::string Measure;
	std::string::size_type sz;     // alias of size_t
	Matrix<double,6,1> z;
	Vector3d u;
	while (std::getline(in, Measure))
	{
		z(0) = std::stod(Measure, &sz);
		Measure = Measure.substr(sz);
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
		u(0) = std::stod(Measure, &sz);
		Measure = Measure.substr(sz);
		u(1) = std::stod(Measure, &sz);
		Measure = Measure.substr(sz);
		u(2) = std::stod(Measure, &sz);
		filter.step(u, z);
		out << filter.getX().transpose() << std::endl;
	}
	std::cout << "OK!"<<std::endl;
}