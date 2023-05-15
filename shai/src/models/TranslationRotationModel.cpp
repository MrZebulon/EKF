//
// Created by Samuel on 15/05/2023.
//

#include "models/TranslationRotationModel.hpp"
#include "Utils.hpp"

#define NX 17
#define NW 7
#define NZ 1

using namespace shai::models;

VectorXd TranslationRotationModel::get_init_state() {
	VectorXd x_init = VectorXd::Zero(NX);
	x_init(6) = 1;
	x_init.segment<3>(10) << accel.bias;
	x_init.segment<3>(13) << gyro.bias;
	x_init.segment<1>(16) << baro.bias;
	return x_init;
}

MatrixXd TranslationRotationModel::get_init_cov() {
	MatrixXd P_init = 1e-9 * MatrixXd::Identity(NX, NX);
	return P_init;
}

VectorXd TranslationRotationModel::compute_x_new(const VectorXd &x, const VectorXd &u) {
	VectorXd x_new = VectorXd::Zero(NX);
	Quaterniond q{x.segment<4>(6)};
	Quaterniond dq{{0, 0.5 * u.segment<3>(3) * _dt - x.segment<3>(13)}};

	x_new.segment<3>(0) = x.segment<3>(0) + x.segment<3>(3) * _dt;
	x_new.segment<3>(3) = q.toRotationMatrix() * (u.segment<3>(0) * _dt - x.segment<3>(10));
	x_new.segment<7>(10)  = x.segment<7>(10);

	q *= dq;
	q.normalize();
	x_new.segment<4>(6) = q.coeffs();

	return x_new;
}

Matrix<double, 3, 4> acc_rotmat_jac(const Vector4d& q, const Vector3d& a) {
	Matrix<double, 3, 4> mat = Matrix<double, 3, 4>::Zero();
	double q0 = q(0), q1 = q(1), q2 = q(2), q3 = q(3);
	double ax = a(0), ay = a(1), az = a(2);
	mat <<  2*ax*q0 - 2*ay*q3 + 2*az*q2, 2*ax*q1 + 2*ay*q2 + 2*az*q3, 2*ay*q1 - 2*ax*q2 + 2*az*q0, 2*az*q1 - 2*ay*q0 - 2*ax*q3,
			2*ax*q3 + 2*ay*q0 - 2*az*q1, 2*ax*q2 - 2*ay*q1 - 2*az*q0, 2*ax*q1 + 2*ay*q2 + 2*az*q3, 2*ax*q0 - 2*ay*q3 + 2*az*q2,
			2*ay*q1 - 2*ax*q2 + 2*az*q0, 2*ax*q3 + 2*ay*q0 - 2*az*q1, 2*ay*q3 - 2*ax*q0 - 2*az*q2, 2*ax*q1 + 2*ay*q2 + 2*az*q3;

	return mat;
}

MatrixXd TranslationRotationModel::get_F_matrix(const VectorXd &x, const VectorXd &u) {
	MatrixXd F = MatrixXd::Zero(NX, NX);
	Quaterniond q{x.segment<4>(6)};
	double w_x = u(3) * _dt, w_y = u(4) * _dt, w_z = u(5) * _dt;

	F.block<3, 3>(0, 3) << MatrixXd::Identity(3, 3) * _dt;
	F.block<3, 4>(3, 6) << acc_rotmat_jac(q.coeffs(), u.segment<3>(0) * _dt);
	F.block<3, 3>(3, 10)<< - (q.toRotationMatrix());

	F.row(6) << 0, 0, 0, 0, 0, 0, 0, x(13)/2 - w_x/2, x(14)/2 - w_y/2, x(15)/2 - w_z/2, 0, 0, 0, x(7)/2, x(8)/2, x(9)/2, 0;
	F.row(7) << 0, 0, 0, 0, 0, 0, w_x/2 - x(13)/2, 0, w_z/2 - x(15)/2, x(14)/2 - w_y/2, 0, 0, 0, -x(6)/2, x(9)/2, -x(8)/2, 0;
	F.row(8) << 0, 0, 0, 0, 0, 0, w_y/2 - x(14)/2, x(15)/2 - w_z/2, 0, w_x/2 - x(13)/2, 0, 0, 0, -x(9)/2, -x(6)/2, x(7)/2, 0;
	F.row(9) << 0, 0, 0, 0, 0, 0, w_z/2 - x(15)/2, w_y/2 - x(14)/2, x(13)/2 - w_x/2, 0, 0, 0, 0, x(8)/2, -x(7)/2, -x(6)/2, 0;
	return F;
}

MatrixXd TranslationRotationModel::get_Q_matrix(const VectorXd &x, const VectorXd &u, const VectorXd &w) {
	MatrixXd G = MatrixXd::Zero(NX, NW);
	Quaterniond q{x.segment<4>(6)};

	G.block<3, 3>(3, 0) << q.toRotationMatrix();
	G.block<4, 4>(6, 3) = 0.5 * hamilton_product_as_matrix(x.segment<4>(6));

	MatrixXd W = MatrixXd::Zero(NW, NW);
	W.diagonal() << w.segment<3>(0), 0, w.segment<3>(3);
	return  G * W * G.transpose();
}

VectorXd TranslationRotationModel::get_measurement_estimate(const VectorXd &x) {
	VectorXd z = VectorXd::Zero(NZ);
	z << x(2);
	return z;
}

MatrixXd TranslationRotationModel::get_H_matrix() {
	MatrixXd H = MatrixXd::Zero(NZ, NX);
	H(0, 2) = 1;
	return H;
}

MatrixXd TranslationRotationModel::get_R_matrix() {
	MatrixXd R = MatrixXd::Zero(NZ, NZ);
	R << 0.1;
	return R;
}

VectorXd TranslationRotationModel::get_noise_vect() {
	VectorXd w = VectorXd::Zero(NW);
	double Fs = 1/_dt;
	double s = 0.5*(1./(Fs*Fs));

	w.segment<3>(0) << Vector3d::Constant(s * accel.noise);
	w.segment<3>(3) << Vector3d::Constant(s * gyro.noise);
	w.segment<1>(6) << s * baro.noise;
	return w;
}

MatrixXd TranslationRotationModel::get_Qs_matrix() {
	MatrixXd Qs = MatrixXd::Zero(NX, NX);
	double Fs = 1/_dt;
	double s = 0.5*(1./(Fs*Fs));

	Qs.diagonal().segment<10>(0) << Vector<double, 10>::Constant(additive_noise);
	Qs.diagonal().segment<3>(10) << Vector3d::Constant(s * accel.drift);
	Qs.diagonal().segment<3>(13) << Vector3d::Constant(s * gyro.drift);
	Qs.diagonal().segment<1>(16) << s * baro.drift;

	return Qs;
}
