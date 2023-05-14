//
// Created by Samuel on 04/05/2023.
//

#ifndef SHAI_UTILS_HPP
#define SHAI_UTILS_HPP

#include <eigen3/Eigen/Eigen>

inline void body_to_inertial(const Eigen::Vector4d& q, const Eigen::Vector3d& vect_b, Eigen::Vector3d& out){
	Eigen::Quaterniond quat(q(0), q(1), q(2), q(3));
	quat.normalize();
	Eigen::Matrix3d rot_matrix = quat.toRotationMatrix();
	out << rot_matrix * vect_b;
}

inline void rotmat(const Eigen::Vector4d& q, Eigen::Matrix3d& out) {
	double q0 = q(0);
	double q1 = q(1);
	double q2 = q(2);
	double q3 = q(3);

	out << q0*q0+q1*q1-q2*q2-q3*q3, 2*(q1*q2-q0*q3), 2*(q1*q3+q0*q2),
			2*(q1*q2+q0*q3), q0*q0-q1*q1+q2*q2-q3*q3, 2*(q2*q3-q0*q1),
			2*(q1*q3-q0*q2), 2*(q2*q3+q0*q1), q0*q0-q1*q1-q2*q2+q3*q3;
}

inline void mult_quat(const Eigen::Vector4d& q1, const Eigen::Vector4d& q2, Eigen::Vector4d& out){
	double qnew_0 = q1(0) * q2(0) - q1(1) * q2(1) - q1(2) * q2(2) - q1(3) * q2(3);
	double qnew_1 = q1(0) * q2(1) + q2(0) * q1(1) + q1(2) * q2(3) - q2(2) * q1(3);
	double qnew_2 = q1(0) * q2(2) + q2(0) * q1(2) - q1(1) * q2(3) + q2(1) * q1(3);
	double qnew_3 = q1(0) * q2(3) + q2(0) * q1(3) + q1(1) * q2(2) - q2(1) * q1(2);

	out << qnew_0, qnew_1, qnew_2, qnew_3;
}

template<std::size_t e>
inline double fuse_data_on_u(const Eigen::Vector<double, e>& data_sources, const Eigen::Vector<double, e>& noises) {
	double out = 0;
	double tot_noise = 0;

	for (size_t i = 0; i < e; i++){
		out += data_sources(i) * 1 / noises(i);
		tot_noise += 1 / noises(i);
	}

	return out/tot_noise;
}

inline void hamilton_product_as_matrix(const Eigen::Quaternion<double>& q, Eigen::MatrixXd& out){
	out <<  q.w(), -q.x(), -q.y(), -q.z(),
			q.x(), q.w(), -q.z(), q.y(),
			q.y(), q.z(), q.w(), -q.x(),
			q.z(), -q.y(), q.x(), q.w();
}

#endif //SHAI_UTILS_HPP
