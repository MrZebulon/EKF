//
// Created by Samuel on 04/05/2023.
//

#ifndef SHAI_UTILS_HPP
#define SHAI_UTILS_HPP

#include "eigen3/Eigen/Eigen"

void body_to_inertial(const Eigen::Vector4d& q, const Eigen::Vector3d& vect_b, Eigen::Vector3d& out){
	Eigen::Matrix3d rot_mat;
	
	rot_mat <<  q(0)*q(0)+q.w()*q.w()-q.x()*q.x()-q.y()*q.y(), 2*(q.w()*q.x()-q(0)*q.y()), 2*(q.w()*q.y()+q(0)*q.x()),
				2*(q.w()*q.x()+q(0)*q.y()), q(0)*q(0)-q.w()*q.w()+q.x()*q.x()-q.y()*q.y(), 2*(q.x()*q.y()-q(0)*q.w()),
				2*(q.w()*q.y()-q(0)*q.x()), 2*(q.x()*q.y()+q(0)*q.w()), q(0)*q(0)-q.w()*q.w()-q.x()*q.x()+q.y()*q.y();

	out = rot_mat * vect_b;
}

void mult_quat(const Eigen::Vector4d& q.x(), const Eigen::Vector4d& q.y(), Eigen::Vector4d& out) {
	out <<  q.x()(0)*q.y()(0) - q.x()(1)*q.y()(1) - q.x()(2)*q.y()(2) - q.x()(3)*q.y()(3),
			q.x()(0)*q.y()(1) + q.y()(0)*q.x()(1) + q.x()(2)*q.y()(3) - q.y()(2)*q.x()(3),
			q.x()(0)*q.y()(2) + q.y()(0)*q.x()(2) - q.x()(1)*q.y()(3) + q.y()(1)*q.x()(3),
			q.x()(0)*q.y()(3) + q.y()(0)*q.x()(3) + q.x()(1)*q.y()(2) - q.y()(1)*q.x()(2);
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

void hamilton_product_as_matrix(const Eigen::Quaternion<double>& q, Eigen::MatrixXd& out){
	out <<  q.w(), -q.x(), -q.y(), -q.z(),
			q.x(), q.w(), -q.z(), q.y(),
			q.y(), q.z(), q.w(), -q.x(),
			q.z(), -q.y(), q.x(), q.w();
}

void rotation_matrix(const Eigen::Quaternion<double>& q, Eigen::MatrixXd& out){

	out <<  q.w()*q.w()+q.x()*q.x()-q.y()*q.y()-q.z()*q.z(), 2*(q.x()*q.y()-q.w()*q.z()), 2*(q.x()*q.z()+q.w()*q.y()),
			2*(q.x()*q.y()+q.w()*q.z()), q.w()*q.w()-q.x()*q.x()+q.y()*q.y()-q.z()*q.z(), 2*(q.y()*q.z()-q.w()*q.x()),
			2*(q.x()*q.z()-q.w()*q.y()), 2*(q.y()*q.z()+q.w()*q.x()), q.w()*q.w()-q.x()*q.x()-q.y()*q.y()+q.z()*q.z();
}

#endif //SHAI_UTILS_HPP
