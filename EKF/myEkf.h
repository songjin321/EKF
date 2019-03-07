#pragma once
#include <Eigen\Dense>
#define PI 3.14159265359;
using namespace Eigen;
class myEkf{
public:
	myEkf();

	/**
	*��ʼ��x��P
	*@param x0 ״̬��ʼֵ
	*@param P0 ������ʼֵ
	*/
	void initialize(Matrix<double, 10, 1> &x0, Matrix<double, 10, 10> &P0);

	/**
	*����ʱ����º��������
	*/
	void step(const Vector3d &u, const Matrix<double, 6, 1> &z);

	/**
	*����ϵͳ״̬��x
	*/
	Matrix<double, 10, 1> getX();
private:
	void makeA();
	void makeH();
	void makeV();
	void makeR();
	void makeW();
	void makeQ();
	void makeProcess();
	void makeMeasure();

	void timeUpdate(const Vector3d &u);
	void measureUpdate(const Matrix<double, 6, 1> &z);
	//�������ʼ����������

	const double Ts = 0.01,g = 9.81,hx = 0.26,hz = 0.37;
	const double sigma_a = 0.05,sigma_g = 0.4 / 180 * PI
	const double sigma_h = 0.001;
	const double b_sigma_g = 0.01 / 180 * PI 
	const double b_sigma_h = 0.01;

	//�˴���10�Ƿ������Nsta����

	Matrix<double,10,1> x = MatrixXd::Zero(10,1);
	
	Matrix<double, 10, 10> P = MatrixXd::Zero(10, 10);

	Matrix<double, 10, 10> A = MatrixXd::Zero(10, 10);

	Matrix<double, 9, 9> Q = MatrixXd::Zero(9,9);

	Matrix<double, 10, 9> W = MatrixXd::Zero(10, 9);

	Matrix<double, 6, 10> H = MatrixXd::Zero(6, 10);

	Matrix<double, 6, 6> R = MatrixXd::Zero(6, 6);
	
	Matrix<double, 6, 6> V = MatrixXd::Zero(6, 6);

	Vector3d u;

	Matrix<double, 6, 1> z;

	Matrix<double, 6, 1> h;
};

