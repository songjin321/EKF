使用EKF融合惯性测量单元和磁力计估算姿态角
==========================
# 文件说明

### myEkf.cpp 

Ekf类源文件，初始化后，给定传感器数据即可估计姿态

### Alignment.cpp 

包含一个采用FAQ算法初始对准函数，估算初始姿态和地磁矢量大小

### test.cpp 

利用仿真数据进行算法正确性测试

### main.ino  

Arduino主函数，读取传感器数据，调用Alignment初始化，然后调用EKF解算姿态

# 原理说明

参见本人[毕业设计论文](宋瑾-毕业论文-最终.pdf)3.1小结

# 演示视频

<a href="https://www.youtube.com/watch?v=5l8blxg7x_w"> <img src="https://img.youtube.com/vi/5l8blxg7x_w/maxresdefault.jpg" alt="video" width="600"/>
</a>