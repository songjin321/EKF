/*
 Name:		main.ino
 Created:	2017/4/9 13:29:49
 Author:	SJ
*/

// the setup function runs once when you press reset or power the board

#include "myEkf.h"
#include "Alignment.h"

int16_t accelData[3], gyroData[3], magData[3]; //加速度、陀螺仪、磁力计值

myEkf filter;

//定义采样间隔，单位ms
const int Ts = 10;
//程序执行时间
unsigned long times = 0;
//初始化过程需要的时间,单位ms
const int initialTime = 15000;
//SYSTEM_MODE(MODE_SEMI_AUTOMATIC);
void setup()
{	
	//初始化串口
	Serial.begin(256000);
	//初始化传感器
	BMI160.begin();
//	IntoRobot.disconnect();
}
void loop()
{
	Vector3d acc, gyro, mag;
	Vector3d accSum(0, 0, 0);
	Vector3d magSum(0, 0, 0);
	unsigned int count = 0;
	while (1)
	{
		times = millis();
		if (times % Ts == 0)
		{
			Serial.printf("%d\r\n", times);
			BMI160.getAccelGyroMagData(accelData, gyroData, magData);
			unitTransform(accelData, gyroData, magData, acc, gyro, mag);

			if (times < initialTime)
			{
				accSum += acc;
				magSum += mag;
				count++;
			}
			if (times == initialTime)
			{
				//执行初始化过程，初始化滤波器
				Vector3d average_acc = accSum / count;
				Vector3d average_mag = magSum / count;
				Quaterniond q0 = Alignment::FAQ(average_acc, average_mag);
				Vector3d estMag = q0._transformVector(average_mag);

				Serial.print(*estMag.data());

				VectorXd P0Vector(10);
				P0Vector << 0.001, 0.001, 0.001, 0.001, 1e-4, 1e-4, 1e-4, 0.01, 0.01, 0.01;
				Matrix<double, 10, 10> P0(P0Vector.asDiagonal());
				Matrix<double, 10, 1> x0 = Matrix<double, 10, 1>::Zero();
				x0(0) = q0.w();
				x0(1) = q0.x();
				x0(2) = q0.y();
				x0(3) = q0.z();
				filter.setMag(estMag);
				filter.initialize(x0, P0);
			}
			if (times > initialTime)
			{
				VectorXd z;
				VectorXd u;
				z << acc, mag;
				u << gyro;
				//printRawData(u, z);
				filter.step(u, z);
				Matrix<double, 10, 1> X = filter.getX();
				Quaterniond q = {   };
				Vector3d euler = q.toRotationMatrix().eulerAngles(2, 1, 0);
				printResult(euler);
			}
			
		}
	}
}
void printResult(const Vector3d &euler)
{
	double Yaw = euler(0) * 180 / PI;
	double Pitch = euler(1) * 180 / PI;
	double Roll = euler(2) * 180 / PI;
//	Serial.printf("Yaw = %f Pitch = %f Roll = %f", Yaw, Pitch, Roll);
	Serial.print(Yaw);
	Serial.printf("\r\n");
}
void printRawDataUSB(Matrix<double, 6, 1> &z,Vector3d &u)
{
	int time = millis();
	SerialUSB.print(time);
	SerialUSB.printf(" ");
	SerialUSB.print(z(0),6);
	SerialUSB.printf(" ");
	SerialUSB.print(z(1),6);
	SerialUSB.printf(" ");
	SerialUSB.print(z(2), 6);
	SerialUSB.printf(" ");
	SerialUSB.print(z(3), 6);
	SerialUSB.printf(" ");
	SerialUSB.print(z(4), 6);
	SerialUSB.printf(" ");
	SerialUSB.print(z(5), 6);
	SerialUSB.printf(" ");
	SerialUSB.print(u(0), 6);
	SerialUSB.printf(" ");
	SerialUSB.print(u(1), 6);
	SerialUSB.printf(" ");
	SerialUSB.print(u(2), 6);
	SerialUSB.printf("\n");
}
void printRawData(VectorXd &u, VectorXd &z)
{
	double * _z = z.data();
	double * _u = u.data();
	Serial.printf("%d %f %f %f %f %f %f %f %f %f\r\n", times, *_z, *(_z+1), *(_z + 2), *(_z + 3), *(_z + 4)
		, *(_z + 5), *_u, *(_u+1), *(_u+2));
}
void printInitial(Quaterniond &q0, Vector3d &mag)
{

}

/** \breif 转化单位
  *
  */
void unitTransform(int16_t accelData[3], int16_t gyroData[3], int16_t magData[3],
	Vector3d &acc, Vector3d &gyro, Vector3d &mag)
{
	//转化为m/s/s
	const double accUnit = 9.81 / 2048;

	//转化为Gauss
	const double magUnitX = 2600 / 32768 * 0.1;
	const double magUnitY = 2600 / 32768 * 0.1;
	const double magUnitZ = 5000 / 32768 * 0.1;

	//转化为rad/s
	const double gyroUnit = M_PI/ (16.4 * 180);

	acc[0] = accelData[0] * accUnit;
	acc[1] = accelData[1] * accUnit;
	acc[2] = accelData[2] * accUnit;
	gyro[0] = gyroData[0] * gyroUnit;
	gyro[1] = gyroData[1] * gyroUnit;
	gyro[2] = gyroData[2] * gyroUnit;
	mag[0] = magData[0] * magUnitX;
	mag[1] = magData[1] * magUnitY;
	mag[2] = magData[2] * magUnitZ;
}


