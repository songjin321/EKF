#pragma once
#include <kalman\ekfilter.hpp>
#include <Eigen\Dense>
#define PI 3.14159265359;
#define Nsta 10;
#define Mobs  6;
using namespace Eigen;
class myEkf{
public:
	myEkf();

private:
	void makeA();
	void makeH();
	void makeV();
	void makeR();
	void makeW();
	void makeQ();
	void makeProcess();
	void makeMeasure();

	//定义与初始化常量参数

	const double Ts = 0.01,g = 9.81,hx = 0.26,hz = 0.37;
	const double sigma_a = 0.005,sigma_g = 0.4 / 180 * PI
	const double sigma_h = 0.001;
	const double b_sigma_g = 0.01 / 180 * PI 
	const double b_sigma_h = 0.001;

};
