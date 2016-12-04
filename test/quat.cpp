// gcc -o quat quat.cpp -std=c++11 -lstdc++ -lm

#include "../math/quat.h"
#include <iostream>

using namespace tl;
typedef float T;

int main()
{
	T angle = 1.23;
	ublas::matrix<T> mat = rotation_matrix_3d_x(angle);
	std::cout << "mat = " << mat << std::endl;

	math::quaternion<T> quat = rot3_to_quat(mat);
	std::cout << "quat = " << quat << std::endl;

	mat = quat_to_rot3(quat);
	std::cout << "mat = " << mat << std::endl;

	std::cout << quat_to_cmat(quat) << std::endl;


	ublas::vector<T> vec = make_vec({1., 2., 3.});
	std::cout << ublas::prod(mat, vec) << std::endl;
	std::cout << quat_vec_prod(quat, vec) << std::endl;


	auto vecEuler = quat_to_euler(quat);
	std::cout << "q from euler: " 
		<< euler_to_quat(vecEuler[0],vecEuler[1],vecEuler[2]) 
		<< std::endl;


	ublas::matrix<T> mat2 = rotation_matrix_3d_y(angle);
	math::quaternion<T> quat2 = rot3_to_quat(mat2);

	std::cout << "prod1: " << quat_prod(quat, quat2) << std::endl;
	std::cout << "prod2: " << quat * quat2 << std::endl;
	return 0;
}
