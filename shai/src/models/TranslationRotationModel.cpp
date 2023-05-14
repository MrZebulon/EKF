//
// Created by Samuel on 14/05/2023.
//

#include "models/TranslationRotationModel.hpp"
#include "Utils.hpp"

using namespace shai::models;

Eigen::VectorXd TranslationRotationModel::get_init_state() {
	Eigen::VectorXd x_init(nx);
	x_init <<   0, 0, 0,
				0, 0, 0,
				1, 0, 0, 0,
				_accel_params.bias[0], _accel_params.bias[1], _accel_params.bias[2],
				_gyro_params.bias[0], _gyro_params.bias[1], _gyro_params.bias[2],
				_baro_params.bias[0];

	return x_init;
}

Eigen::MatrixXd TranslationRotationModel::get_init_cov() {
	return 1e-9 * Eigen::MatrixXd::Identity(nx, nx);
}

void TranslationRotationModel::compute_x_new(const Eigen::VectorXd &x,
											 const Eigen::VectorXd &u,
											 Eigen::VectorXd &out) {

	// p_new = p_old + v_old * dt
	Eigen::Vector3d p_new = x.segment<3>(0) + x.segment<3>(3) * _dt;

	//v_new = v_old + R(q) * (a*dt - bias)
	Eigen::Vector3d v_new;
	body_to_inertial(x.segment<4>(6), u.segment<3>(0) * _dt - x.segment<3>(10), v_new);
	v_new += x.segment<3>(3);

	//q_new = 1/2 * q_old * dq
	Eigen::Vector4d dq;
	dq << 1., 0.5 * (u.segment<3>(3)*_dt - x.segment<3>(13));
	Eigen::Vector4d q_new;
	mult_quat(x.segment<4>(6), dq, q_new);
	q_new.normalize();

	Eigen::Vector<double, 7> bias_new = x.segment<7>(10);
	out << p_new, v_new, q_new, bias_new;

}

void TranslationRotationModel::get_F_matrix(const Eigen::VectorXd &x,
											const Eigen::VectorXd &u,
											Eigen::MatrixXd &out) {

	double q0 = x(6), q1 = x(7), q2 = x(8), q3 = x(9);
	double b_acc_x = x(10), b_acc_y = x(11), b_acc_z = x(12);
	double b_gyro_x = x(13), b_gyro_y = x(14), b_gyro_z = x(15);

	double a_x = u(0) * _dt, a_y = u(1) * _dt, a_z = u(2) * _dt;
	double w_x = u(3) * _dt, w_y = u(4) * _dt, w_z = u(5) * _dt;

	out.block<3, 3>(0, 3) = Eigen::Matrix3d::Identity() * _dt;

	out.row(3) << 0, 0, 0, 0, 0, 0,
			2 * q0 * (a_x - b_acc_x) - 2 * q3 * (a_y - b_acc_y) + 2 * q2 * (a_z - b_acc_z),
			2 * q1 * (a_x - b_acc_x) + 2 * q2 * (a_y - b_acc_y) + 2 * q3 * (a_z - b_acc_z),
			2 * q1 * (a_y - b_acc_y) - 2 * q2 * (a_x - b_acc_x) + 2 * q0 * (a_z - b_acc_z),
			2 * q1 * (a_z - b_acc_z) - 2 * q0 * (a_y - b_acc_y) - 2 * q3 * (a_x - b_acc_x),
			0, 0, 0,
			-q0 * q0 - q1 * q1 + q2 * q2 + q3 * q3,
			2 * q0 * q3 - 2 * q1 * q2,
			-2 * q0 * q2 - 2 * q1 * q3,
			0;

	out.row(4) << 0, 0, 0, 0, 0, 0,
			2 * q0 * (a_x - b_acc_x) - 2 * q3 * (a_y - b_acc_y) + 2 * q2 * (a_z - b_acc_z),
			2 * q1 * (a_x - b_acc_x) + 2 * q2 * (a_y - b_acc_y) + 2 * q3 * (a_z - b_acc_z),
			2 * q1 * (a_y - b_acc_y) - 2 * q2 * (a_x - b_acc_x) + 2 * q0 * (a_z - b_acc_z),
			2 * q1 * (a_z - b_acc_z) - 2 * q0 * (a_y - b_acc_y) - 2 * q3 * (a_x - b_acc_x),
			0, 0, 0,
			-q0 * q0 - q1 * q1 + q2 * q2 + q3 * q3,
			2 * q0 * q3 - 2 * q1 * q2,
			-2 * q0 * q2 - 2 * q1 * q3, 0;

	out.row(5) << 0, 0, 0, 0, 0, 0,
			2 * q1 * (a_y - b_acc_y) - 2 * q2 * (a_x - b_acc_x) + 2 * q0 * (a_z - b_acc_z),
			2 * q3 * (a_x - b_acc_x) + 2 * q0 * (a_y - b_acc_y) - 2 * q1 * (a_z - b_acc_z),
			2 * q3 * (a_y - b_acc_y) - 2 * q0 * (a_x - b_acc_x) - 2 * q2 * (a_z - b_acc_z),
			2 * q1 * (a_x - b_acc_x) + 2 * q2 * (a_y - b_acc_y) + 2 * q3 * (a_z - b_acc_z),
			0, 0, 0,
			2 * q0 * q2 - 2 * q1 * q3,
			-2 * q0 * q1 - 2 * q2 * q3,
			-q0 * q0 + q1 * q1 + q2 * q2 - q3 * q3, 0;

	out.row(6) << 0, 0, 0,  0, 0, 0,
			0,
			w_x / 2 - b_gyro_x / 2,
			b_gyro_y / 2 - w_y / 2,
			b_gyro_z / 2 - w_z / 2,
			0, 0, 0,
			q1 / 2, q2 / 2, q3 / 2,
			0;

	out.row(7) << 0, 0, 0,  0, 0, 0,
			w_x / 2 - b_gyro_x / 2,
			0,
			w_z / 2 - b_gyro_z / 2,
			b_gyro_y / 2 - w_y / 2,
			0, 0, 0,
			-q0 / 2, q3 / 2, -q2 / 2,
			0;

	out.row(8) << 0, 0, 0, 0, 0, 0,
			w_y / 2 - b_gyro_y / 2,
			b_gyro_z / 2 - w_z / 2,
			0,
			w_x / 2 - b_gyro_x / 2,
			0, 0, 0,
			-q3 / 2, -q0 / 2, q1 / 2,
			0;

	out.row(9) << 0, 0, 0, 0, 0, 0,
			w_z / 2 - b_gyro_z / 2,
			w_y / 2 - b_gyro_y / 2,
			b_gyro_x / 2 - w_x / 2,
			0,
			0, 0, 0,
			q2 / 2, -q1 / 2, -q0 / 2,
			0;
}

void TranslationRotationModel::get_Q_matrix(const Eigen::VectorXd &x,
											const Eigen::VectorXd &u,
											const Eigen::VectorXd &w,
											Eigen::MatrixXd &out) {
	double q0 = x(6), q1 = x(7), q2 = x(8), q3 = x(9);

	double dvxCov = w(0), dvyCov = w(1), dvzCov = w(2);
	double daxCov = w(3), dayCov = w(4), dazCov = w(5);

	out.row(3) << 0, 0, 0,
			dvyCov * std::pow(2 * q0 * q3 - 2 * q1 * q2, 2) + dvzCov * std::pow(2 * q0 * q2 + 2 * q1 * q3, 2) + dvxCov * std::pow(q0 * q0 + q1 * q1 - q2 * q2 - q3 * q3, 2),
			dvxCov * (2 * q0 * q3 + 2 * q1 * q2) * (q0 * q0 + q1 * q1 - q2 * q2 - q3 * q3) - dvyCov * (2 * q0 * q3 - 2 * q1 * q2) * (q0 * q0 - q1 * q1 + q2 * q2 - q3 * q3) - dvzCov * (2 * q0 * q1 - 2 * q2 * q3) * (2 * q0 * q2 + 2 * q1 * q3),
			dvzCov * (2 * q0 * q2 + 2 * q1 * q3) * (q0 * q0 - q1 * q1 - q2 * q2 + q3 * q3) - dvxCov * (2 * q0 * q2 - 2 * q1 * q3) * (q0 * q0 + q1 * q1 - q2 * q2 - q3 * q3) - dvyCov * (2 * q0 * q1 + 2 * q2 * q3) * (2 * q0 * q3 - 2 * q1 * q2),
			0, 0, 0, 0,  0, 0, 0,  0, 0, 0,  0;

	out.row(4) << 0, 0, 0,
			dvxCov * (2 * q0 * q3 + 2 * q1 * q2) * (q0 * q0 + q1 * q1 - q2 * q2 - q3 * q3) - dvyCov * (2 * q0 * q3 - 2 * q1 * q2) * (q0 * q0 - q1 * q1 + q2 * q2 - q3 * q3) - dvzCov * (2 * q0 * q1 - 2 * q2 * q3) * (2 * q0 * q2 + 2 * q1 * q3),
			dvxCov * (2 * q0 * q3 + 2 * q1 * q2) * (2 * q0 * q3 + 2 * q1 * q2) + dvzCov * (2 * q0 * q1 - 2 * q2 * q3) * (2 * q0 * q1 - 2 * q2 * q3) + dvyCov * (q0 * q0 - q1 * q1 + q2 * q2 - q3 * q3) * (q0 * q0 - q1 * q1 + q2 * q2 - q3 * q3),
			dvyCov * (2 * q0 * q1 + 2 * q2 * q3) * (q0 * q0 - q1 * q1 + q2 * q2 - q3 * q3) - dvzCov * (2 * q0 * q1 - 2 * q2 * q3) * (q0 * q0 - q1 * q1 - q2 * q2 + q3 * q3) - dvxCov * (2 * q0 * q2 - 2 * q1 * q3) * (2 * q0 * q3 + 2 * q1 * q2),
			0, 0, 0, 0,  0, 0, 0,  0, 0, 0,  0;

	out.row(5) << 0, 0, 0,
			dvzCov * (2 * q0 * q2 + 2 * q1 * q3) * (q0 * q0 - q1 * q1 - q2 * q2 + q3 * q3) - dvxCov * (2 * q0 * q2 - 2 * q1 * q3) * (q0 * q0 + q1 * q1 - q2 * q2 - q3 * q3) - dvyCov * (2 * q0 * q1 + 2 * q2 * q3) * (2 * q0 * q3 - 2 * q1 * q2),
			dvyCov * (2 * q0 * q1 + 2 * q2 * q3) * (q0 * q0 - q1 * q1 + q2 * q2 - q3 * q3) - dvzCov * (2 * q0 * q1 - 2 * q2 * q3) * (q0 * q0 - q1 * q1 - q2 * q2 + q3 * q3) - dvxCov * (2 * q0 * q2 - 2 * q1 * q3) * (2 * q0 * q3 + 2 * q1 * q2),
			dvxCov * std::pow(2 * q0 * q2 - 2 * q1 * q3, 2) + dvyCov * std::pow(2 * q0 * q1 + 2 * q2 * q3, 2) + dvzCov * std::pow(q0 * q0 - q1 * q1 - q2 * q2 + q3 * q3, 2),
			0, 0, 0, 0,  0, 0, 0,  0, 0, 0,  0;

	out.row(6) << 0, 0, 0,  0, 0, 0,
			(daxCov * q1 * q1) / 4 + (dayCov * q2 * q2) / 4 + (dazCov * q3 * q3) / 4,
			(dayCov * q2 * q3) / 4 - (daxCov * q0 * q1) / 4 - (dazCov * q2 * q3) / 4,
			(dazCov * q1 * q3) / 4 - (dayCov * q0 * q2) / 4 - (daxCov * q1 * q3) / 4,
			(daxCov * q1 * q2) / 4 - (dayCov * q1 * q2) / 4 - (dazCov * q0 * q3) / 4,
			0, 0, 0,  0, 0, 0,  0;

	out.row(7) << 0, 0, 0, 0, 0, 0,
			(dayCov * q2 * q3) / 4 - (daxCov * q0 * q1) / 4 - (dazCov * q2 * q3) / 4,
			(daxCov * q0 * q0) / 4 + (dazCov * q2 * q2) / 4 + (dayCov * q3 * q3) / 4,
			(daxCov * q0 * q3) / 4 - (dayCov * q0 * q3) / 4 - (dazCov * q1 * q2) / 4,
			(dazCov * q0 * q2) / 4 - (dayCov * q1 * q3) / 4 - (daxCov * q0 * q2) / 4,
			0, 0, 0,  0, 0, 0,  0;

	out.row(8) << 0, 0, 0, 0, 0, 0, (dazCov * q1 * q3) / 4 - (dayCov * q0 * q2) / 4 - (daxCov * q1 * q3) / 4,
			(daxCov * q0 * q3) / 4 - (dayCov * q0 * q3) / 4 - (dazCov * q1 * q2) / 4,
			(dayCov * q0 * q0) / 4 + (dazCov * q1 * q1) / 4 + (daxCov * q3 * q3) / 4,
			(dayCov * q0 * q1) / 4 - (daxCov * q2 * q3) / 4 - (dazCov * q0 * q1) / 4,
			0, 0, 0,  0, 0, 0,  0;

	out.row(9) << 0, 0, 0, 0, 0, 0, (daxCov * q1 * q2) / 4 - (dayCov * q1 * q2) / 4 - (dazCov * q0 * q3) / 4,
			(dazCov * q0 * q2) / 4 - (dayCov * q1 * q3) / 4 - (daxCov * q0 * q2) / 4,
			(dayCov * q0 * q1) / 4 - (daxCov * q2 * q3) / 4 - (dazCov * q0 * q1) / 4,
			(dazCov * q0 * q0) / 4 + (dayCov * q1 * q1) / 4 + (daxCov * q2 * q2) / 4,
			0, 0, 0,  0, 0, 0,  0;
}

void TranslationRotationModel::get_measurement_estimate(const Eigen::VectorXd &x,
														Eigen::VectorXd &out) {
	out << x(2) + x(16);
}

void TranslationRotationModel::get_H_matrix(Eigen::MatrixXd &out) {
	out << 0, 0, 1,  0, 0, 0,  0, 0, 0, 0,  0, 0, 0,  0, 0, 0,  0;
}

void TranslationRotationModel::get_R_matrix(Eigen::MatrixXd &out) {
	out << 0.1;
}

void TranslationRotationModel::get_noise_vect(Eigen::VectorXd &out) {
	double Fs = 1/_dt;
	out << Eigen::Vector3d::Constant(_accel_params.noise), Eigen::Vector3d::Constant(_gyro_params.noise), _baro_params.noise;
	out *= 0.5*(1./(Fs*Fs));
}

void TranslationRotationModel::get_Qs_matrix(Eigen::MatrixXd &out) {
	double Fs = 1 / _dt;
	double scale_var = 0.5 * (1 /(Fs*Fs));

	double vel = scale_var * _accel_params.drift;
	double ang = scale_var * _gyro_params.drift;
	double pos = scale_var * _baro_params.drift;
	out.setZero();
	out.diagonal() << Eigen::VectorXd::Ones(10), Eigen::Vector3d::Constant(vel), Eigen::Vector3d::Constant(ang), pos;
}