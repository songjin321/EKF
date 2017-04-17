#include "myEkf.h"

myEkf::myEkf()
{
	setDim(10, 3, 9, 6, 6);
	Ts = 0.00025;
	g = 9.81; 
	hx = 0.26;
	hz = 0.37;
	sigma_a = 0.005;
	sigma_g = 0.4 / 180 * PI;
	sigma_h = 0.001;
	b_sigma_g = 0.01 / 180 * PI;
	b_sigma_h = 0.001;
}

void myEkf::makeA()
{
	double q1 = x(1);
	double q2 = x(2);
	double q3 = x(3);
	double q4 = x(4);
	double thx = (x(8) - u(1))*Ts;
	double thy = (x(9) - u(2))*Ts;
	double thz = (x(10) - u(3))*Ts;
	double th = sqrt(thx * thx + thy * thy + thz *thz);
	double cth = cos(th / 2);
	double sth = sin(th / 2);
	double sin_th_x = (sth * thx) / th;
	double sin_th_y = (sth * thy) / th;
	double sin_th_z = (sth * thz) / th;
	A(1, 1) = cth;
	A(1, 2) = sin_th_x;
	A(1, 3) = sin_th_y;
	A(1, 4) = sin_th_z;
	A(1, 5) = 0.0;
	A(1, 6) = 0.0;
	A(1, 7) = 0.0;
	A(1, 8) = (Ts*q2*sth) / th + (Ts*cth*q2*(thx*thx)) / (2 * (th*th)) - (Ts*q2*sth*(thx*thx)) / (th*th*th) -
			  (Ts*q1*sth*thx) / (2 * th) + (Ts*cth*q3*thx*thy) / (2 * (th*th)) + (Ts*cth*q4*thx*thz) / (2 * (th*th)) -
			  (Ts*q3*sth*thx*thy) / (th*th*th) - (Ts*q4*sth*thx*thz) / (th*th*th);
	A(1, 9) = (Ts*q3*sth) / th + (Ts*cth*q3*(thy*thy)) / (2 * (th*th)) - (Ts*q3*sth*(thy*thy)) / (th*th*th) - (Ts*q1*sth*thy) / (2 * th) +
			  (Ts*cth*q2*thx*thy) / (2 * (th*th)) + (Ts*cth*q4*thy*thz) / (2 * (th*th)) - (Ts*q2*sth*thx*thy) / (th*th*th) -
			  (Ts*q4*sth*thy*thz) / (th*th*th);
	A(1, 10) = (Ts*q4*sth) / th + (Ts*cth*q4*(thz*thz)) / (2 * (th*th)) - (Ts*q4*sth*(thz*thz)) / (th*th*th) - (Ts*q1*sth*thz) / (2 * th) +
			  (Ts*cth*q2*thx*thz) / (2 * (th*th)) + (Ts*cth*q3*thy*thz) / (2 * (th*th)) - (Ts*q2*sth*thx*thz) / (th*th*th) -
			  (Ts*q3*sth*thy*thz) / (th*th*th);

	A(2, 1) = -sin_th_x;
	A(2, 2) = cth;
	A(2, 3) = -sin_th_z;
	A(2, 4) = sin_th_y;
	A(2, 5) = 0.0;
	A(2, 6) = 0.0;
	A(2, 7) = 0.0;
	A(2, 8) = (Ts*q1*sth*(thx*thx)) / (th*th*th) - (Ts*cth*q1*(thx*thx)) / (2 * (th*th)) - (Ts*q1*sth) / th - (Ts*q2*sth*thx) / (2 * th) +
			  (Ts*cth*q4*thx*thy) / (2 * (th*th)) - (Ts*cth*q3*thx*thz) / (2 * (th*th)) - (Ts*q4*sth*thx*thy) / (th*th*th) +
			  (Ts*q3*sth*thx*thz) / (th*th*th);
	A(2, 9) = (Ts*q4*sth) / th + (Ts*cth*q4*(thy*thy)) / (2 * (th*th)) - (Ts*q4*sth*(thy*thy)) / (th*th*th) - (Ts*q2*sth*thy) / (2 * th) -
			  (Ts*cth*q1*thx*thy) / (2 * (th*th)) - (Ts*cth*q3*thy*thz) / (2 * (th*th)) + (Ts*q1*sth*thx*thy) / (th*th*th) +
			  (Ts*q3*sth*thy*thz) / (th*th*th);
	A(2, 10) = (Ts*q3*sth*(thz*thz)) / (th*th*th) - (Ts*cth*q3*(thz*thz)) / (2 * (th*th)) - (Ts*q3*sth) / th - (Ts*q2*sth*thz) / (2 * th) -
			  (Ts*cth*q1*thx*thz) / (2 * (th*th)) + (Ts*cth*q4*thy*thz) / (2 * (th*th)) + (Ts*q1*sth*thx*thz) / (th*th*th) -
			  (Ts*q4*sth*thy*thz) / (th*th*th);

	A(3, 1) = -sin_th_y;
	A(3, 2) = sin_th_z;
	A(3, 3) = cth;
	A(3, 4) = -sin_th_x;
	A(3, 5) = 0.0;
	A(3, 6) = 0.0;
	A(3, 7) = 0.0;
	A(3, 8) = (Ts*q4*sth*(thx*thx)) / (th*th*th) - (Ts*cth*q4*(thx*thx)) / (2 * (th*th)) - (Ts*q4*sth) / th - (Ts*q3*sth*thx) / (2 * th) -
			  (Ts*cth*q1*thx*thy) / (2 * (th*th)) + (Ts*cth*q2*thx*thz) / (2 * (th*th)) + (Ts*q1*sth*thx*thy) / (th*th*th) -
			  (Ts*q2*sth*thx*thz) / (th*th*th);
	A(3, 9) = (Ts*q1*sth*(thy*thy)) / (th*th*th) - (Ts*cth*q1*(thy*thy)) / (2 * (th*th)) - (Ts*q1*sth) / th - (Ts*q3*sth*thy) / (2 * th) -
			  (Ts*cth*q4*thx*thy) / (2 * (th*th)) + (Ts*cth*q2*thy*thz) / (2 * (th*th)) + (Ts*q4*sth*thx*thy) / (th*th*th) -
			  (Ts*q2*sth*thy*thz) / (th*th*th);
	A(3, 10) = (Ts*q2*sth) / th + (Ts*cth*q2*(thz*thz)) / (2 * (th*th)) - (Ts*q2*sth*(thz*thz)) / (th*th*th) - (Ts*q3*sth*thz) / (2 * th) -
			  (Ts*cth*q4*thx*thz) / (2 * (th*th)) - (Ts*cth*q1*thy*thz) / (2 * (th*th)) + (Ts*q4*sth*thx*thz) / (th*th*th) +
			  (Ts*q1*sth*thy*thz) / (th*th*th);

	A(4, 1) = -sin_th_z;
	A(4, 2) = -sin_th_y;
	A(4, 3) = sin_th_x;
	A(4, 4) = cth;
	A(4, 5) = 0.0;
	A(4, 6) = 0.0;
	A(4, 7) = 0.0;
	A(4, 8) = (Ts*q3*sth) / th + (Ts*cth*q3*(thx*thx)) / (2 * (th*th)) - (Ts*q3*sth*(thx*thx)) / (th*th*th) - (Ts*q4*sth*thx) / (2 * th) -
			  (Ts*cth*q2*thx*thy) / (2 * (th*th)) - (Ts*cth*q1*thx*thz) / (2 * (th*th)) + (Ts*q2*sth*thx*thy) / (th*th*th) +
			  (Ts*q1*sth*thx*thz) / (th*th*th);
	A(4, 9) = (Ts*q2*sth*(thy*thy)) / (th*th*th) - (Ts*cth*q2*(thy*thy)) / (2 * (th*th)) - (Ts*q2*sth) / th - (Ts*q4*sth*thy) / (2 * th) +
			  (Ts*cth*q3*thx*thy) / (2 * (th*th)) - (Ts*cth*q1*thy*thz) / (2 * (th*th)) - (Ts*q3*sth*thx*thy) / (th*th*th) +
			  (Ts*q1*sth*thy*thz) / (th*th*th);
	A(4, 10) = (Ts*q1*sth*(thz*thz)) / (th*th*th) - (Ts*cth*q1*(thz*thz)) / (2 * (th*th)) - (Ts*q1*sth) / th - (Ts*q4*sth*thz) / (2 * th) +
			  (Ts*cth*q3*thx*thz) / (2 * (th*th)) - (Ts*cth*q2*thy*thz) / (2 * (th*th)) - (Ts*q3*sth*thx*thz) / (th*th*th) +
			  (Ts*q2*sth*thy*thz) / (th*th*th);

	for (size_t i = 5; i < 11; i++)
	{
		for (size_t j = 1; j < 11; j++)
		{
			if (i <= 4 && j <= 4)
				break;
			if (i == j)
				A(i, j) = 1.0;
			else
				A(i, j) = 0.0;
		}
	}
}

void myEkf::makeH()
{
	double q1 = x(1);
	double q2 = x(2);
	double q3 = x(3);
	double q4 = x(4);
	double hbx = x(5);
	double hby = x(6);
	double hbz = x(7);

	H(1, 1) = -2 * g*q3;
	H(1, 2) = 2 * g*q4;
	H(1, 3) = -2 * g*q1;
	H(1, 4) = 2 * g*q2;
	H(1, 5) = 0.0;
	H(1, 6) = 0.0;
	H(1, 7) = 0.0;
	H(1, 8) = 0.0;
	H(1, 9) = 0.0;
	H(1, 10) = 0.0;

	H(2, 1) = 2 * g*q2;
	H(2, 2) = 2 * g*q1;
	H(2, 3) = 2 * g*q4;
	H(2, 4) = 2 * g*q3;
	H(2, 5) = 0.0;
	H(2, 6) = 0.0;
	H(2, 7) = 0.0;
	H(2, 8) = 0.0;
	H(2, 9) = 0.0;
	H(2, 10) = 0.0;

	H(3, 1) = 2 * g*q1;
	H(3, 3) = -2 * g*q2;
	H(3, 3) = -2 * g*q3;
	H(3, 4) = 2 * g*q4;
	H(3, 5) = 0.0;
	H(3, 6) = 0.0;
	H(3, 7) = 0.0;
	H(3, 8) = 0.0;
	H(3, 9) = 0.0;
	H(3, 10) = 0.0;

	H(4, 1) = 2 * q1*(hbx + hx) - 2 * q3*(hbz + hz) + 2 * hby*q4;
	H(4, 4) = 2 * q2*(hbx + hx) + 2 * q4*(hbz + hz) + 2 * hby*q3;
	H(4, 4) = 2 * hby*q2 - 2 * q1*(hbz + hz) - 2 * q3*(hbx + hx);
	H(4, 4) = 2 * q2*(hbz + hz) - 2 * q4*(hbx + hx) + 2 * hby*q1;
	H(4, 5) = 0.0;
	H(4, 6) = 0.0;
	H(4, 7) = 0.0;
	H(4, 8) = q1*q1 + q2*q2 - q3*q3 - q4*q4;
	H(4, 9) = 2 * q2*q3 + 2 * q1*q4;
	H(4, 10) = 2 * q2*q4 - 2 * q1*q3;

	H(5, 1) = 2 * hby*q1 + 2 * q2*(hbz + hz) - 2 * q4*(hbx + hx);
	H(5, 5) = 2 * q3*(hbx + hx) + 2 * q1*(hbz + hz) - 2 * hby*q2;
	H(5, 5) = 2 * q2*(hbx + hx) + 2 * q4*(hbz + hz) + 2 * hby*q3;
	H(5, 5) = 2 * q3*(hbz + hz) - 2 * q1*(hbx + hx) - 2 * hby*q4;
	H(5, 5) = 0.0;
	H(5, 6) = 0.0;
	H(5, 7) = 0.0;
	H(5, 8) = 2 * q2*q3 - 2 * q1*q4;
	H(5, 9) = q1*q1 - q2*q2 + q3*q3 - q4*q4;
	H(5, 10) = 2 * q3*q4 + 2 * q1*q2;

	H(6, 1) = 2 * q3*(hbx + hx) + 2 * q1*(hbz + hz) - 2 * hby*q2;
	H(6, 6) = 2 * q4*(hbx + hx) - 2 * q2*(hbz + hz) - 2 * hby*q1;
	H(6, 6) = 2 * q1*(hbx + hx) - 2 * q3*(hbz + hz) + 2 * hby*q4;
	H(6, 6) = 2 * q2*(hbx + hx) + 2 * q4*(hbz + hz) + 2 * hby*q3;
	H(6, 6) = 0.0;
	H(6, 6) = 0.0;
	H(6, 7) = 0.0;
	H(6, 8) = 2 * q1*q3 + 2 * q2*q4;
	H(6, 9) = 2 * q3*q4 - 2 * q1*q2;
	H(6, 10) = q1*q1 - q2*q2 - q3*q3 + q4*q4;
}

void myEkf::makeV()
{
	V(1, 1) = 1.0;
	V(2, 2) = 1.0;
	V(3, 3) = 1.0;
	V(4, 4) = 1.0;
	V(5, 5) = 1.0;
	V(6, 6) = 1.0;
}

void myEkf::makeR()
{
	R(1, 1) = sigma_a * sigma_a;
	R(2, 2) = sigma_a * sigma_a;
	R(3, 3) = sigma_a * sigma_a;

	R(4, 4) = sigma_h * sigma_h;
	R(5, 5) = sigma_h * sigma_h;
	R(6, 6) = sigma_h * sigma_h;
}

void myEkf::makeW()
{
	double q1 = x(1);
	double q2 = x(2);
	double q3 = x(3);
	double q4 = x(4);
	W(1, 1) = -Ts / 2 * (-q2);
	W(1, 2) = -Ts / 2 * (-q3);
	W(1, 3) = -Ts / 2 * (-q4);
	W(2, 1) = -Ts / 2 * (q1);
	W(2, 2) = -Ts / 2 * (-q4);
	W(2, 3) = -Ts / 2 * (q3);
	W(3, 1) = -Ts / 2 * (q4);
	W(3, 2) = -Ts / 2 * (q1);
	W(3, 3) = -Ts / 2 * (-q2);
	W(4, 1) = -Ts / 2 * (-q3);
	W(4, 2) = -Ts / 2 * (q2);
	W(4, 3) = -Ts / 2 * (q1);
	for (size_t i = 1; i < 11; i++)
	{
		for (size_t j = 1; j < 10; j++)
		{
			if (i <= 4 && j <= 3)
				break;
			if (i == j)
				A(i, j) = sqrt(Ts);
			else
				A(i, j) = 0.0;
		}
	}
}

void myEkf::makeQ()
{
	Q(1, 1) = sigma_g * sigma_g;
	Q(2, 2) = sigma_g * sigma_g;
	Q(3, 3) = sigma_g * sigma_g;

	Q(4, 4) = b_sigma_h * b_sigma_h;
	Q(5, 5) = b_sigma_h * b_sigma_h;
	Q(6, 6) = b_sigma_h * b_sigma_h;

	Q(7, 7) = b_sigma_g * b_sigma_g;
	Q(8, 8) = b_sigma_g * b_sigma_g;
	Q(9, 9) = b_sigma_g * b_sigma_g;
}

void myEkf::makeProcess()
{
	double q1 = x(1);
	double q2 = x(2);
	double q3 = x(3);
	double q4 = x(4);
	double thx = (x(8) - u(1))*Ts;
	double thy = (x(9) - u(2))*Ts;
	double thz = (x(10) - u(3))*Ts;
	double th =  sqrt(thx * thx + thy * thy + thz *thz);
	double cth = cos(th / 2);
	double sth = sin(th / 2);
	double sin_th_x = (sth * thx) / th;
	double sin_th_y = (sth * thy) / th;
	double sin_th_z = (sth * thz) / th;
	Vector x_(x.size());
	x_(1) = cth * q1 + sin_th_x * q2 + sin_th_y * q3 + sin_th_z * q4;
	x_(2) = -sin_th_x * q1 + cth * q2 - sin_th_z * q3 + sin_th_y * q4;
	x_(3) = -sin_th_y * q1 + sin_th_z * q2 + cth * q3 - sin_th_x * q4;
	x_(4) = -sin_th_z * q1 - sin_th_y * q2 + sin_th_x * q3 + cth * q4;
	x_(5) = x(5);
	x_(6) = x(6);
	x_(7) = x(7);
	x_(8) = x(8);
	x_(9) = x(9);
	x_(10) = x(10);
	x.swap(x_);
}

void myEkf::makeMeasure()
{
	double q1 = x(1);
	double q2 = x(2);
	double q3 = x(3);
	double q4 = x(4);
	double hbx = x(5);
	double hby = x(6);
	double hbz = x(7);

	z(1) = -g*(2 * q1*q3 - 2 * q2*q4);
	z(2) = g*(2 * q1*q2 + 2 * q3*q4);
	z(3) = g*(q1*q1 - q2*q2 - q3*q3 + q4*q4);
	z(4) = (hbx + hx)*(q1*q1 + q2*q2 - q3*q3 - q4*q4) + hby*(2 * q1*q4 + 2 * q2*q3) - (hbz + hz)*(2 * q1*q3 - 2 * q2*q4);
	z(5) = hby*(q1*q1 - q2*q2 + q3*q3 - q4*q4) - (hbx + hx)*(2 * q1*q4 - 2 * q2*q3) + (hbz + hz)*(2 * q1*q2 + 2 * q3*q4);
	z(6) = (hbx + hx)*(2 * q1*q3 + 2 * q2*q4) - hby*(2 * q1*q2 - 2 * q3*q4) + (hbz + hz)*(q1*q1 - q2*q2 - q3*q3 + q4*q4);
}
