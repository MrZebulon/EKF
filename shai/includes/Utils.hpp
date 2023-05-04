//
// Created by Samuel on 04/05/2023.
//

#ifndef SHAI_UTILS_HPP
#define SHAI_UTILS_HPP

#include "eigen3/Eigen/Eigen"

void body_to_inertial(const Eigen::Vector4d& q, const Eigen::Vector3d& vect_b, Eigen::Vector3d& out){
	Eigen::Matrix3d rot_mat;
	
	rot_mat <<  q(0)*q(0)+q(1)*q(1)-q(2)*q(2)-q(3)*q(3), 2*(q(1)*q(2)-q(0)*q(3)), 2*(q(1)*q(3)+q(0)*q(2)),
				2*(q(1)*q(2)+q(0)*q(3)), q(0)*q(0)-q(1)*q(1)+q(2)*q(2)-q(3)*q(3), 2*(q(2)*q(3)-q(0)*q(1)),
				2*(q(1)*q(3)-q(0)*q(2)), 2*(q(2)*q(3)+q(0)*q(1)), q(0)*q(0)-q(1)*q(1)-q(2)*q(2)+q(3)*q(3);

	out = rot_mat * vect_b;
}

void mult_quat(const Eigen::Vector4d& q1, const Eigen::Vector4d& q2, Eigen::Vector4d& out) {
	out <<  q1(0)*q2(0) - q1(1)*q2(1) - q1(2)*q2(2) - q1(3)*q2(3),
			q1(0)*q2(1) + q2(0)*q1(1) + q1(2)*q2(3) - q2(2)*q1(3),
			q1(0)*q2(2) + q2(0)*q1(2) - q1(1)*q2(3) + q2(1)*q1(3),
			q1(0)*q2(3) + q2(0)*q1(3) + q1(1)*q2(2) - q2(1)*q1(2);
}

template<std::size_t e>
double fuse_data_on_u(const Eigen::Vector<double, e>& data_sources, const Eigen::Vector<double, e>& noises) {
	double out = 0;
	double tot_noise = 0;

	for (size_t i = 0; i < e; i++){
		out += data_sources(i) * 1 / noises(i);
		tot_noise += 1 / noises(i);
	}

	return out/tot_noise;
}

#endif //SHAI_UTILS_HPP
