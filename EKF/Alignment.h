#pragma once
/**使用 FQA方法解算初始姿态
*/
#include <Eigen\Dense>
class Alignment
{
public:
	Alignment();
	static Eigen::Quaterniond FAQ(const Eigen::Vector3d & acc, const Eigen::Vector3d & mag);
	~Alignment();
};

