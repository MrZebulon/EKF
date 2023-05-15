//
// Created by Samuel on 15/05/2023.
//

#ifndef SHAI_UTILS_HPP
#define SHAI_UTILS_HPP

#include <eigen3/Eigen/Core>

inline Eigen::Matrix4d hamilton_product_as_matrix(const Eigen::Vector4d& q) {
	Eigen::Matrix4d mat;
	mat << q(0), -q(1), -q(2), -q(3),
			q(1), q(0), -q(3), q(2),
			q(2), q(3), q(0), -q(1),
			q(3), -q(2), q(1), q(0);
	return mat;
}


#endif //SHAI_UTILS_HPP
