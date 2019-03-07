#include "myEkf.h"
#include <Eigen/Dense>
myEkf::myEkf()
{
}

void myEkf::initialize(Matrix<double, 10, 1>& x0, Matrix<double, 10, 10>& P0)
{
	x = x0;
	P = P0;
}

void myEkf::step(const Vector3d & u, const Matrix<double, 6, 1> &z)
{
	timeUpdate(u);
	measureUpdate(z);
}

Matrix<double, 10, 1> myEkf::getX()
{
	return x;
}

void myEkf::makeA()
{
	double q1 = x(0);
	double q2 = x(1);
	double q3 = x(2);
	double q4 = x(3);
	double thx = (x(7) - u(0))*Ts;
	double thy = (x(8) - u(1))*Ts;
	double thz = (x(9) - u(2))*Ts;
	double th = sqrt(thx * thx + thy * thy + thz *thz);
	double cth = cos(th / 2);
	double sth = sin(th / 2);
	double sin_th_x = (sth * thx) / th;
	double sin_th_y = (sth * thy) / th;
	double sin_th_z = (sth * thz) / th;
	A(0, 0) = cth;
	A(0, 1) = sin_th_x;
	A(0, 2) = sin_th_y;
	A(0, 3) = sin_th_z;
	A(0, 4) = 0.0;
	A(0, 5) = 0.0;
	A(0, 6) = 0.0;
	A(0, 7) = (Ts*q2*sth) / th - (Ts*q1*thx*sth) / (2 * th) + (Ts*q2*thx*cos(th / 2)*thx) / (2 * th*th) + (Ts*q3*thx*cos(th / 2)*thy) / (2 * th*th) + (Ts*q4*thx*cos(th / 2)*thz) / (2 * th*th) - (Ts*q2*thx*sth*thx) / (th*th*th) - (Ts*q3*thx*sth*thy) / (th*th*th) - (Ts*q4*thx*sth*thz) / (th*th*th);
	A(0, 8) = (Ts*q3*sth) / th - (Ts*q1*thy*sth) / (2 * th) + (Ts*q2*thy*cos(th / 2)*thx) / (2 * th*th) + (Ts*q3*thy*cos(th / 2)*thy) / (2 * th*th) + (Ts*q4*thy*cos(th / 2)*thz) / (2 * th*th) - (Ts*q2*thy*sth*thx) / (th*th*th) - (Ts*q3*thy*sth*thy) / (th*th*th) - (Ts*q4*thy*sth*thz) / (th*th*th);
	A(0, 9) = (Ts*q4*sth) / th - (Ts*q1*thz*sth) / (2 * th) + (Ts*q2*thz*cos(th / 2)*thx) / (2 * th*th) + (Ts*q3*thz*cos(th / 2)*thy) / (2 * th*th) + (Ts*q4*thz*cos(th / 2)*thz) / (2 * th*th) - (Ts*q2*thz*sth*thx) / (th*th*th) - (Ts*q3*thz*sth*thy) / (th*th*th) - (Ts*q4*thz*sth*thz) / (th*th*th);

	A(1, 0) = -sin_th_x;
	A(1, 1) = cth;
	A(1, 2) = -sin_th_z;
	A(1, 3) = sin_th_y;
	A(1, 4) = 0.0;
	A(1, 5) = 0.0;
	A(1, 6) = 0.0;
	A(1, 7) = (Ts*q4*thx*cos(th / 2)*thy) / (2 * th*th) - (Ts*q2*thx*sth) / (2 * th) - (Ts*q1*thx*cos(th / 2)*thx) / (2 * th*th) - (Ts*q1*sth) / th - (Ts*q3*thx*cos(th / 2)*thz) / (2 * th*th) + (Ts*q1*thx*sth*thx) / (th*th*th) - (Ts*q4*thx*sth*thy) / (th*th*th) + (Ts*q3*thx*sth*thz) / (th*th*th);
	A(1, 8) = (Ts*q4*sth) / th - (Ts*q2*thy*sth) / (2 * th) - (Ts*q1*thy*cos(th / 2)*thx) / (2 * th*th) + (Ts*q4*thy*cos(th / 2)*thy) / (2 * th*th) - (Ts*q3*thy*cos(th / 2)*thz) / (2 * th*th) + (Ts*q1*thy*sth*thx) / (th*th*th) - (Ts*q4*thy*sth*thy) / (th*th*th) + (Ts*q3*thy*sth*thz) / (th*th*th);
	A(1, 9) = (Ts*q4*thz*cos(th / 2)*thy) / (2 * th*th) - (Ts*q2*thz*sth) / (2 * th) - (Ts*q1*thz*cos(th / 2)*thx) / (2 * th*th) - (Ts*q3*sth) / th - (Ts*q3*thz*cos(th / 2)*thz) / (2 * th*th) + (Ts*q1*thz*sth*thx) / (th*th*th) - (Ts*q4*thz*sth*thy) / (th*th*th) + (Ts*q3*thz*sth*thz) / (th*th*th);

	A(2, 0) = -sin_th_y;
	A(2, 1) = sin_th_z;
	A(2, 2) = cth;
	A(2, 3) = -sin_th_x;
	A(2, 4) = 0.0;
	A(2, 5) = 0.0;
	A(2, 6) = 0.0;
	A(2, 7) = (Ts*q2*thx*cos(th / 2)*thz) / (2 * th*th) - (Ts*q3*thx*sth) / (2 * th) - (Ts*q4*thx*cos(th / 2)*thx) / (2 * th*th) - (Ts*q1*thx*cos(th / 2)*thy) / (2 * th*th) - (Ts*q4*sth) / th + (Ts*q4*thx*sth*thx) / (th*th*th) + (Ts*q1*thx*sth*thy) / (th*th*th) - (Ts*q2*thx*sth*thz) / (th*th*th);
	A(2, 8) = (Ts*q2*thy*cos(th / 2)*thz) / (2 * th*th) - (Ts*q3*thy*sth) / (2 * th) - (Ts*q4*thy*cos(th / 2)*thx) / (2 * th*th) - (Ts*q1*thy*cos(th / 2)*thy) / (2 * th*th) - (Ts*q1*sth) / th + (Ts*q4*thy*sth*thx) / (th*th*th) + (Ts*q1*thy*sth*thy) / (th*th*th) - (Ts*q2*thy*sth*thz) / (th*th*th);
	A(2, 9) = (Ts*q2*sth) / th - (Ts*q3*thz*sth) / (2 * th) - (Ts*q4*thz*cos(th / 2)*thx) / (2 * th*th) - (Ts*q1*thz*cos(th / 2)*thy) / (2 * th*th) + (Ts*q2*thz*cos(th / 2)*thz) / (2 * th*th) + (Ts*q4*thz*sth*thx) / (th*th*th) + (Ts*q1*thz*sth*thy) / (th*th*th) - (Ts*q2*thz*sth*thz) / (th*th*th);

	A(3, 0) = -sin_th_z;
	A(3, 1) = -sin_th_y;
	A(3, 2) = sin_th_x;
	A(3, 3) = cth;
	A(3, 4) = 0.0;
	A(3, 5) = 0.0;
	A(3, 6) = 0.0;
	A(3, 7) = (Ts*q3*sth) / th - (Ts*q4*thx*sth) / (2 * th) + (Ts*q3*thx*cos(th / 2)*thx) / (2 * th*th) - (Ts*q2*thx*cos(th / 2)*thy) / (2 * th*th) - (Ts*q1*thx*cos(th / 2)*thz) / (2 * th*th) - (Ts*q3*thx*sth*thx) / (th*th*th) + (Ts*q2*thx*sth*thy) / (th*th*th) + (Ts*q1*thx*sth*thz) / (th*th*th);
	A(3, 8) = (Ts*q3*thy*cos(th / 2)*thx) / (2 * th*th) - (Ts*q4*thy*sth) / (2 * th) - (Ts*q2*sth) / th - (Ts*q2*thy*cos(th / 2)*thy) / (2 * th*th) - (Ts*q1*thy*cos(th / 2)*thz) / (2 * th*th) - (Ts*q3*thy*sth*thx) / (th*th*th) + (Ts*q2*thy*sth*thy) / (th*th*th) + (Ts*q1*thy*sth*thz) / (th*th*th);
	A(3, 9) = (Ts*q3*thz*cos(th / 2)*thx) / (2 * th*th) - (Ts*q4*thz*sth) / (2 * th) - (Ts*q1*sth) / th - (Ts*q2*thz*cos(th / 2)*thy) / (2 * th*th) - (Ts*q1*thz*cos(th / 2)*thz) / (2 * th*th) - (Ts*q3*thz*sth*thx) / (th*th*th) + (Ts*q2*thz*sth*thy) / (th*th*th) + (Ts*q1*thz*sth*thz) / (th*th*th);

	for (size_t i = 4; i < 10; i++)
	{
		for (size_t j = 0; j < 10; j++)
		{
			if (i == j)
				A(i, j) = 1.0;
			else
				A(i, j) = 0.0;
		}
	}
}

void myEkf::makeH()
{
	double q1 = x(0);
	double q2 = x(1);
	double q3 = x(2);
	double q4 = x(3);
	double hbx = x(4);
	double hby = x(5);
	double hbz = x(6);

	H(0, 0) = -2 * g*q3;
	H(0, 1) = 2 * g*q4;
	H(0, 2) = -2 * g*q1;
	H(0, 3) = 2 * g*q2;
	H(0, 4) = 0.0;
	H(0, 5) = 0.0;
	H(0, 6) = 0.0;
	H(0, 7) = 0.0;
	H(0, 8) = 0.0;
	H(0, 9) = 0.0;

	H(1, 0) = 2 * g*q2;
	H(1, 1) = 2 * g*q1;
	H(1, 2) = 2 * g*q4;
	H(1, 3) = 2 * g*q3;
	H(1, 4) = 0.0;
	H(1, 5) = 0.0;
	H(1, 6) = 0.0;
	H(1, 7) = 0.0;
	H(1, 8) = 0.0;
	H(1, 9) = 0.0;

	H(2, 0) = 2 * g*q1;
	H(2, 1) = -2 * g*q2;
	H(2, 2) = -2 * g*q3;
	H(2, 3) = 2 * g*q4;
	H(2, 4) = 0.0;
	H(2, 5) = 0.0;
	H(2, 6) = 0.0;
	H(2, 7) = 0.0;
	H(2, 8) = 0.0;
	H(2, 9) = 0.0;

	H(3, 0) = 2 * q1*(hbx + hx) - 2 * q3*(hbz + hz) + 2 * hby*q4;
	H(3, 1) = 2 * q2*(hbx + hx) + 2 * q4*(hbz + hz) + 2 * hby*q3;
	H(3, 2) = 2 * hby*q2 - 2 * q1*(hbz + hz) - 2 * q3*(hbx + hx);
	H(3, 3) = 2 * q2*(hbz + hz) - 2 * q4*(hbx + hx) + 2 * hby*q1;
	H(3, 4) = q1*q1 + q2*q2 - q3*q3 - q4*q4;
	H(3, 5) = 2 * q2*q3 + 2 * q1*q4;
	H(3, 6) = 2 * q2*q4 - 2 * q1*q3;
	H(3, 7) = 0.0;
	H(3, 8) = 0.0;
	H(3, 9) = 0.0;

	H(4, 0) = 2 * hby*q1 + 2 * q2*(hbz + hz) - 2 * q4*(hbx + hx);
	H(4, 1) = 2 * q3*(hbx + hx) + 2 * q1*(hbz + hz) - 2 * hby*q2;
	H(4, 2) = 2 * q2*(hbx + hx) + 2 * q4*(hbz + hz) + 2 * hby*q3;
	H(4, 3) = 2 * q3*(hbz + hz) - 2 * q1*(hbx + hx) - 2 * hby*q4;
	H(4, 4) = 2 * q2*q3 - 2 * q1*q4;
	H(4, 5) = q1*q1 - q2*q2 + q3*q3 - q4*q4;
	H(4, 6) = 2 * q3*q4 + 2 * q1*q2;
	H(4, 7) = 0.0;
	H(4, 8) = 0.0;
	H(4, 9) = 0.0;

	H(5, 0) = 2 * q3*(hbx + hx) + 2 * q1*(hbz + hz) - 2 * hby*q2;
	H(5, 1) = 2 * q4*(hbx + hx) - 2 * q2*(hbz + hz) - 2 * hby*q1;
	H(5, 2) = 2 * q1*(hbx + hx) - 2 * q3*(hbz + hz) + 2 * hby*q4;
	H(5, 3) = 2 * q2*(hbx + hx) + 2 * q4*(hbz + hz) + 2 * hby*q3;
	H(5, 4) = 2 * q1*q3 + 2 * q2*q4;
	H(5, 5) = 2 * q3*q4 - 2 * q1*q2;
	H(5, 6) = q1*q1 - q2*q2 - q3*q3 + q4*q4;
	H(5, 7) = 0.0;
	H(5, 8) = 0.0;
	H(5, 9) = 0.0;
}

void myEkf::makeV()
{
	V(0, 0) = 1.0;
	V(1, 1) = 1.0;
	V(2, 2) = 1.0;
	V(3, 3) = 1.0;
	V(4, 4) = 1.0;
	V(5, 5) = 1.0;
}

void myEkf::makeR()
{
	R(0, 0) = sigma_a * sigma_a;
	R(1, 1) = sigma_a * sigma_a;
	R(2, 2) = sigma_a * sigma_a;

	R(3, 3) = sigma_h * sigma_h;
	R(4, 4) = sigma_h * sigma_h;
	R(5, 5) = sigma_h * sigma_h;
}

void myEkf::makeW()
{
	double q1 = x(0);
	double q2 = x(1);
	double q3 = x(2);
	double q4 = x(3);
	W(0, 0) = -Ts / 2 * (-q2);
	W(0, 1) = -Ts / 2 * (-q3);
	W(0, 2) = -Ts / 2 * (-q4);
	W(1, 0) = -Ts / 2 * (q1);
	W(1, 1) = -Ts / 2 * (-q4);
	W(1, 2) = -Ts / 2 * (q3);
	W(2, 0) = -Ts / 2 * (q4);
	W(2, 1) = -Ts / 2 * (q1);
	W(2, 2) = -Ts / 2 * (-q2);
	W(3, 0) = -Ts / 2 * (-q3);
	W(3, 1) = -Ts / 2 * (q2);
	W(3, 2) = -Ts / 2 * (q1);
	W(4, 3) = sqrt(Ts);
	W(5, 4) = sqrt(Ts);
	W(6, 5) = sqrt(Ts);
	W(7, 6) = sqrt(Ts);
	W(8, 7) = sqrt(Ts);
	W(9, 8) = sqrt(Ts);

}

void myEkf::makeQ()
{												
	Q(0, 0) = sigma_g * sigma_g;      
	Q(1, 1) = sigma_g * sigma_g;	  
	Q(2, 2) = sigma_g * sigma_g;	  
									  
	Q(3, 3) = b_sigma_h * b_sigma_h;  
	Q(4, 4) = b_sigma_h * b_sigma_h;  
	Q(5, 5) = b_sigma_h * b_sigma_h;  
									  
	Q(6, 6) = b_sigma_g * b_sigma_g;  
	Q(7, 7) = b_sigma_g * b_sigma_g;  
	Q(8, 8) = b_sigma_g * b_sigma_g;
}

void myEkf::makeProcess()
{
	double q1 = x(0);
	double q2 = x(1);
	double q3 = x(2);
	double q4 = x(3);
	double thx = (x(7) - u(0))*Ts;
	double thy = (x(8) - u(1))*Ts;
	double thz = (x(9) - u(2))*Ts;
	double th =  sqrt(thx * thx + thy * thy + thz *thz);
	double cth = cos(th / 2);
	double sth = sin(th / 2);
	double sin_th_x = (sth * thx) / th;
	double sin_th_y = (sth * thy) / th;
	double sin_th_z = (sth * thz) / th;
	Matrix<double,10,1> x_;
	x_(0) = cth * q1 + sin_th_x * q2 + sin_th_y * q3 + sin_th_z * q4;
	x_(1) = -sin_th_x * q1 + cth * q2 - sin_th_z * q3 + sin_th_y * q4;
	x_(2) = -sin_th_y * q1 + sin_th_z * q2 + cth * q3 - sin_th_x * q4;
	x_(3) = -sin_th_z * q1 - sin_th_y * q2 + sin_th_x * q3 + cth * q4;
	x_(4) = x(4);
	x_(5) = x(5);
	x_(6) = x(6);
	x_(7) = x(7);
	x_(8) = x(8);
	x_(9) = x(9);
	double normq = sqrt(x_(0)*x_(0) + x_(1)*x_(1) + x_(2)*x_(2) + x_(3)*x_(3));
	x_(0) = x_(0) / normq;
	x_(1) = x_(1) / normq;
	x_(2) = x_(2) / normq;
	x_(3) = x_(3) / normq;
	x = x_;
}

void myEkf::makeMeasure()
{
	double q1 = x(0);
	double q2 = x(1);
	double q3 = x(2);
	double q4 = x(3);
	double hbx = x(4);
	double hby = x(5);
	double hbz = x(6);

	h(0) = -g*(2 * q1*q3 - 2 * q2*q4);
	h(1) = g*(2 * q1*q2 + 2 * q3*q4);
	h(2) = g*(q1*q1 - q2*q2 - q3*q3 + q4*q4);
	h(3) = (hbx + hx)*(q1*q1 + q2*q2 - q3*q3 - q4*q4) + hby*(2 * q1*q4 + 2 * q2*q3) - (hbz + hz)*(2 * q1*q3 - 2 * q2*q4);
	h(4) = hby*(q1*q1 - q2*q2 + q3*q3 - q4*q4) - (hbx + hx)*(2 * q1*q4 - 2 * q2*q3) + (hbz + hz)*(2 * q1*q2 + 2 * q3*q4);
	h(5) = (hbx + hx)*(2 * q1*q3 + 2 * q2*q4) - hby*(2 * q1*q2 - 2 * q3*q4) + (hbz + hz)*(q1*q1 - q2*q2 - q3*q3 + q4*q4);
}

void myEkf::timeUpdate(const Vector3d &_u)
{
	u = _u;


	makeA();
	makeW();
	makeQ();
	makeProcess();

	P = A*P*A.transpose() + W*Q*W.transpose();
}


void myEkf::measureUpdate(const Matrix<double, 6, 1>& _z)
{
	z = _z;

	makeMeasure();
	makeH();
	makeV();
	makeR();


	Matrix<double, 10, 6> K = P * H.transpose() * (H*P*H.transpose() + R).inverse();
	x = x + K * (z - h);
	Vector4d q(x(0), x(1), x(2), x(3));
	x(0) = x(0) / q.norm();
	x(1) = x(1) / q.norm();
    x(2) = x(2) / q.norm();
	x(3) = x(3) / q.norm();
	P = (MatrixXd::Identity(10, 10) - K * H)*P;
}
