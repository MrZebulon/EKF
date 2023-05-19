//
// Created by Samuel on 15/05/2023.
//

#include <cmath>
#include <eigen3/Eigen/Geometry>

#include "models/TranslationRotationModel.hpp"

#define NX 17
#define NW 7
#define NZ 1

using namespace shai::models;

Matrix4d hamilton_product_as_matrix(const Vector4d& q){
	Matrix4d mat;
	mat <<  q(0), -q(1), -q(2), -q(3),
			q(1),  q(0), -q(3), q(2),
			q(2), q(3), q(0), -q(1),
			q(3), -q(2), q(1), q(0);

	return mat;
}

Matrix3d rotmat(const Vector4d& q){
	Matrix3d mat;

	double q0 = q(0), q1 = q(1), q2 = q(2), q3 = q(3);

	mat <<  q0*q0+q1*q1-q2*q2-q3*q3, 2*(q1*q2-q0*q3), 2*(q1*q3+q0*q2),
			2*(q1*q2+q0*q3), q0*q0-q1*q1+q2*q2-q3*q3, 2*(q2*q3-q0*q1),
			2*(q1*q3-q0*q2), 2*(q2*q3+q0*q1), q0*q0-q1*q1-q2*q2+q3*q3;
	return mat;
}

Vector4d multquat(const Vector4d& q1, const Vector4d& q2){
	Vector4d q_new;
	q_new(0) = q1(0)*q2(0)    -q1(1)*q2(1)    -q1(2)*q2(2)    -q1(3)*q2(3);
	q_new(0) = q1(0)*q2(1)    +q2(0)*q1(1)    +q1(2)*q2(3)    -q2(2)*q1(3);
	q_new(0) = q1(0)*q2(2)    +q2(0)*q1(2)    -q1(1)*q2(3)    +q2(1)*q1(3);
	q_new(0) = q1(0)*q2(3)    +q2(0)*q1(3)    +q1(1)*q2(2)    -q2(1)*q1(2);
	return q_new;
}

VectorXd TranslationRotationModel::get_init_state() {
	VectorXd x_init = VectorXd::Zero(NX);
	x_init << 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, accel.bias, gyro.bias, baro.bias;
	return x_init;
}

MatrixXd TranslationRotationModel::get_init_cov() {
	MatrixXd P_init = 1e-9 * MatrixXd::Identity(NX, NX);
	return P_init;
}

VectorXd TranslationRotationModel::compute_x_new(const VectorXd &x, const VectorXd &u) {
	Vector3d delta_acc = rotmat(x.segment<4>(6)) * (_dt * u.segment<3>(0) - x.segment<3>(10));
	Vector3d axis = (_dt * u.segment<3>(3) - x.segment<3>(13)) / 2.;

	VectorXd dx(x.size());
	dx << x(3) * _dt,
			x(4) * _dt,
			x(5) * _dt,
			delta_acc(0),
			delta_acc(1),
			delta_acc(2),
			0,
			0,
			0,
			0,
			0,
			0,
			0,
			0,
			0,
			0,
			0;

	VectorXd x_new = x + dx;
	Vector4d delta_q{1, axis(0), axis(1), axis(2)};
	Vector4d q = multquat(x.segment<4>(6), delta_q);
	q.normalize();

	x_new.segment<4>(6) = q;

	return x_new;
}

MatrixXd TranslationRotationModel::get_F_matrix(const VectorXd &x, const VectorXd &u) {
	MatrixXd F = MatrixXd::Zero(NX, NX);
	double q0 = x(6), q1 = x(7), q2 = x(8), q3 = x(9);

	double a_x = u(0) * _dt, a_y = u(1) * _dt, a_z = u(2) * _dt;
	double w_x = u(3) * _dt, w_y = u(4) * _dt, w_z = u(5) * _dt;
	double b_acc_x = x(10), b_acc_y = x(11),b_acc_z = x(12);
	double b_gyro_x = x(13), b_gyro_y = x(14), b_gyro_z = x(15);

	F.block<3, 3>(0, 3) << Matrix3d::Identity() * _dt;

	F.row(3) << 0, 0, 0, 0, 0, 0, 2*q0*(a_x - b_acc_x) - 2*q3*(a_y - b_acc_y) + 2*q2*(a_z - b_acc_z), 2*q1*(a_x - b_acc_x) + 2*q2*(a_y - b_acc_y) + 2*q3*(a_z - b_acc_z), 2*q1*(a_y - b_acc_y) - 2*q2*(a_x - b_acc_x) + 2*q0*(a_z - b_acc_z), 2*q1*(a_z - b_acc_z) - 2*q0*(a_y - b_acc_y) - 2*q3*(a_x - b_acc_x), 0, 0, 0, -q0*q0 - q1*q1 + q2*q2 + q3*q3, 2*q0*q3 - 2*q1*q2, -2*q0*q2 - 2*q1*q3, 0;
	F.row(4) << 0, 0, 0, 0, 0, 0, 2*q3*(a_x - b_acc_x) + 2*q0*(a_y - b_acc_y) - 2*q1*(a_z - b_acc_z), 2*q2*(a_x - b_acc_x) - 2*q1*(a_y - b_acc_y) - 2*q0*(a_z - b_acc_z), 2*q1*(a_x - b_acc_x) + 2*q2*(a_y - b_acc_y) + 2*q3*(a_z - b_acc_z), 2*q0*(a_x - b_acc_x) - 2*q3*(a_y - b_acc_y) + 2*q2*(a_z - b_acc_z), 0, 0, 0, -2*q0*q3 - 2*q1*q2, -q0*q0 + q1*q1 - q2*q2 + q3*q3, 2*q0*q1 - 2*q2*q3, 0;
	F.row(5) << 0, 0, 0, 0, 0, 0, 2*q1*(a_y - b_acc_y) - 2*q2*(a_x - b_acc_x) + 2*q0*(a_z - b_acc_z), 2*q3*(a_x - b_acc_x) + 2*q0*(a_y - b_acc_y) - 2*q1*(a_z - b_acc_z), 2*q3*(a_y - b_acc_y) - 2*q0*(a_x - b_acc_x) - 2*q2*(a_z - b_acc_z), 2*q1*(a_x - b_acc_x) + 2*q2*(a_y - b_acc_y) + 2*q3*(a_z - b_acc_z), 0, 0, 0, 2*q0*q2 - 2*q1*q3, -2*q0*q1 - 2*q2*q3, -q0*q0 + q1*q1 + q2*q2 - q3*q3, 0;

	F.row(6) << 0, 0, 0, 0, 0, 0, 0, b_gyro_x/2 - w_x/2, b_gyro_y/2 - w_y/2, b_gyro_z/2 - w_z/2, 0, 0, 0, q1/2, q2/2, q3/2, 0;
	F.row(7) << 0, 0, 0, 0, 0, 0, w_x/2 - b_gyro_x/2, 0, w_z/2 - b_gyro_z/2, b_gyro_y/2 - w_y/2, 0, 0, 0, -q0/2, q3/2, -q2/2, 0;
	F.row(8) << 0, 0, 0, 0, 0, 0, w_y/2 - b_gyro_y/2, b_gyro_z/2 - w_z/2, 0, w_x/2 - b_gyro_x/2, 0, 0, 0, -q3/2, -q0/2, q1/2, 0;
	F.row(9) << 0, 0, 0, 0, 0, 0, w_z/2 - b_gyro_z/2, w_y/2 - b_gyro_y/2, b_gyro_x/2 - w_x/2, 0, 0, 0, 0, q2/2, -q1/2, -q0/2, 0;

	return F;
}

MatrixXd TranslationRotationModel::get_Q_matrix(const VectorXd &x, const VectorXd &u, const VectorXd &w) {
	MatrixXd Q = MatrixXd::Zero(NX, NX);

	double q0 = x(6), q1 = x(7), q2 = x(8), q3 = x(9);

	double dvxCov = w(0), dvyCov = w(1), dvzCov = w(2);
	double daxCov = w(3), dayCov = w(4), dazCov = w(5);

	MatrixXd G = MatrixXd::Zero(NX, NW);
	G.block<3, 3>(3, 0) = rotmat({q0, q1, q2, q3});
	G.block<4, 4>(6, 3) = 0.5*hamilton_product_as_matrix(x.segment<4>(6));
	Q = G * Vector<double, NW>{dvxCov, dvyCov, dvzCov, 0, daxCov, dayCov, dazCov}.asDiagonal() * G.transpose();

	return Q;
}

VectorXd TranslationRotationModel::get_measurement_estimate(const VectorXd &x) {
	VectorXd z = VectorXd::Zero(NZ);
	z << x(2) + x(16);
	return z;
}

MatrixXd TranslationRotationModel::get_H_matrix() {
	MatrixXd H = MatrixXd::Zero(NZ, NX);
	H << 0, 0, 1,  0, 0, 0,  0, 0, 0, 0,  0, 0, 0,  0, 0, 0,  0;
	return H;
}

MatrixXd TranslationRotationModel::get_R_matrix() {
	MatrixXd R = MatrixXd::Zero(NZ, NZ);
	R << baro.noise;
	return R;
}

VectorXd TranslationRotationModel::get_noise_vect() {
	VectorXd w = VectorXd::Zero(NW);
	double Fs = 1./_dt;
	double s =  0.5 * (1. / (Fs * Fs));

	w << Vector3d::Constant(s * accel.noise), Vector3d::Constant(s * gyro.noise), s * baro.noise;
	return w;
}

MatrixXd TranslationRotationModel::get_Qs_matrix() {
	double Fs = 1./_dt;
	double scale_var = 0.5 *(1./ (Fs * Fs));

	Vector3d vel_delta_bias_sigma = Vector3d::Constant(scale_var * accel.drift);
	Vector3d ang_vel_delta_bias_sigma = Vector3d::Constant(scale_var * gyro.drift);
	double pos_delta_bias_sigma = scale_var * baro.drift;

	Eigen::Vector<double, NX> Qs;
	Qs << Vector<double, 10>::Constant(additive_noise), vel_delta_bias_sigma, ang_vel_delta_bias_sigma, pos_delta_bias_sigma;
	return Qs.asDiagonal();
}
