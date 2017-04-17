#pragma once
#include <kalman\ekfilter.hpp>

#define PI 3.14159265359;



class myEkf : public Kalman::EKFilter<double, 1,true,true,true> {
public:
	myEkf();

protected:
	void makeA();
	void makeH();
	void makeV();
	void makeR();
	void makeW();
	void makeQ();
	void makeProcess();
	void makeMeasure();

	double Ts;
	double g, hx, hz;
	double sigma_a, sigma_g, sigma_h, b_sigma_g, b_sigma_h;
};

typedef myEkf::Vector Vector;
typedef myEkf::Matrix Matrix;

