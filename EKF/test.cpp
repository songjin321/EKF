#include "myEkf.h"
#include "Alignment.h"
#include <fstream>
#include <string>
#include <iostream>
/*此main函数可以测试EKF算法的正确性*/

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
/*此main函数测试FAQ算法准确性*/
/*
int main()
{
	Vector3d a;
	
	Vector3d acc(0, 0, -9.81);
	Vector3d mag(0.37, 0, 0.25);
	Quaterniond q(0.5, 0.5, 0.5, 0.5);
	Vector3d mea_acc = q.inverse()._transformVector(acc);
	Vector3d mea_mag = q.inverse()._transformVector(mag);
	Quaterniond mea_q = Alignment::FAQ(mea_acc, mea_mag);
	Vector3d euler = q.toRotationMatrix().eulerAngles(2, 1, 0);
	MatrixXd mea_euler = mea_q.toRotationMatrix().eulerAngles(2, 1, 0);
	Vector3d est_mag = mea_q._transformVector(mea_mag);
	std::cout << euler << std::endl;
	std::cout << mea_euler << std::endl;
	std::cout << est_mag << std::endl;
}
*/

//此函数利用FAQ计算初始值
/*
int main()
{
	Vector3d average_acc(0.0536963378906250,0.0912023437500000,- 9.65576074218750);
	Vector3d average_mag(-0.0573669433593750,- 0.186145019531250,0.624237060546875);
	Quaterniond q0 = Alignment::FAQ(average_acc, average_mag);
	Vector3d estMag = q0._transformVector(average_mag);
	std::cout << estMag << std::endl;
	std::cout << q0.w() << " " << q0.x() << " " << q0.y() << " " << q0.z() << std::endl;
}
*/