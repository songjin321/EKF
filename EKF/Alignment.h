#pragma once
/**ʹ�� FQA���������ʼ��̬
*/
#include <Eigen\Dense>
class Alignment
{
public:
	Alignment();
	static Eigen::Quaterniond FAQ(const Eigen::Vector3d & acc, const Eigen::Vector3d & mag);
	~Alignment();
};

