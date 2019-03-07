#include "Alignment.h"



Alignment::Alignment()
{
}

Eigen::Quaterniond Alignment::FAQ(const Eigen::Vector3d &acc, const Eigen::Vector3d &bm)
{
	Eigen::Vector3d _acc = acc.normalized();
	double ax = _acc(0);
	double ay = _acc(1);
	double az = _acc(2);

	//计算qe
	double sin_th = ax;
	double cos_th = sqrt(1 - sin_th * sin_th);
	double sin_th_div2(0);
	if (sin_th >= 0)
	{
		sin_th_div2 = sqrt((1 - cos_th) / 2);
	}
	else
	{
		sin_th_div2 = -sqrt((1 - cos_th) / 2);
	}
	double cos_th_div2 = sqrt((1 + cos_th) / 2);
	Eigen::Quaterniond qe(cos_th_div2,0,sin_th_div2,0);

	//计算qr
	Eigen::Quaterniond qr;
	double sin_phi = -ay / cos_th;
	double cos_phi = -az / cos_th;
	double sin_phi_div2(0);
	if (ay == 0 && az == 0)
	{
		qr = { 1,0,0,0 };
	}
	else if (sin_phi >= 0)
	{
		sin_phi_div2 = sqrt((1 - cos_phi) / 2);
	}
	else
	{
		sin_phi_div2 = -sqrt((1 - cos_phi) / 2);
	}
	double cos_phi_div2 = sqrt((1 + cos_phi) / 2);
	qr = { cos_phi_div2, sin_phi_div2, 0, 0 };

	//计算qa 忽略较小的磁偏角
	Eigen::Vector3d em = qe._transformVector(qr._transformVector(bm));
	Eigen::Vector2d mxy(em(0), em(1));
	mxy.normalize();
	double cos_psi = mxy(0);
	double sin_psi = -mxy(1);
	double sin_psi_div2(0);
	if (sin_psi >= 0)
	{
		sin_psi_div2 = sqrt((1 - cos_psi) / 2);
	}
	else
	{
		sin_psi_div2 = -sqrt((1 - cos_psi) / 2);
	}
	double cos_psi_div2 = sqrt((1 + cos_psi) / 2);
	Eigen::Quaterniond qa = { cos_psi_div2, 0, 0, sin_psi_div2 };

	Eigen::Quaterniond q = qa * qe * qr;
	return q;
}


Alignment::~Alignment()
{
}
