//
// Created by Samuel on 06/05/2023.
//

#include <cmath>

#include "TranslationRotationModelV2.hpp"
#include "Utils.hpp"

using namespace std;
using namespace shai::models;
using namespace Eigen;

VectorXd TranslationRotationModelV2::get_init_state() {
	Vector3d ba {_accel_params.bias[0], _accel_params.bias[1], _accel_params.bias[2]};
	Vector3d bg {_gyro_params.bias[0], _gyro_params.bias[1], _gyro_params.bias[2]};
	double bb = _baro_params.bias[0];

	VectorXd x_init(nx);
	x_init << VectorXd::Zero(10), ba, bg, bb;

	return x_init;
}

MatrixXd TranslationRotationModelV2::get_init_cov() {
	MatrixXd P_init(nx, nx);
	P_init.diagonal() << VectorXd::Constant(nx, 1e-9);
	return P_init;
}

void TranslationRotationModelV2::compute_x_new(const VectorXd &x, const VectorXd &u, VectorXd& out) {

	Vector3d pos = {x(0), x(1), x(2)};
	Vector3d vel = {x(3), x(4), x(5)};
	Quaternion<double> q = {x(6), x(7), x(8), x(9)};
	Vector3d ba {x(10), x(11), x(12)};
	Vector3d bg{x(13), x(14), x(15)};
	double bb = x(16);

	Vector3d a = {u(0), u(1), u(2)};
	Vector3d w = {u(3), u(4), u(5)};

	Vector3d delta_v{};

	body_to_inertial({q.w(), q.x(), q.y(), q.z()}, {_dt*a(0) - ba(0), _dt*a(1) - ba(1), _dt*a(2) - ba(2)}, delta_v);
	Quaternion<double> delta_q = {1, 0.5*(_dt*w(0) - bg(0)), 0.5*(_dt*w(1) - bg(1)), 0.5*(_dt*w(2) - bg(2))};

	q *= delta_q;

	out << pos + vel*_dt, vel + delta_v, q.coeffs(), ba, bg, bb;
}

void TranslationRotationModelV2::get_F_matrix(const VectorXd &x, const VectorXd &u, MatrixXd& out) {
	Quaternion<double> q = {x(6), x(7), x(8), x(9)};
	Vector3d ba {x(10), x(11), x(12)};
	Vector3d bg {x(13), x(14), x(15)};

	Vector3d a = {u(0) * _dt, u(1) * _dt, u(2) * _dt};
	Vector3d w = {u(3) * _dt, u(4) * _dt, u(5) * _dt};

	// dp^dot
	out(0, 3) = _dt;
	out(1, 4) = _dt;
	out(2, 5) = _dt;

	//dv^dot
	out.row(3) << 0, 0, 0, 0, 0, 0, 2*q.w()*(a(0) - ba(0)) - 2*q.z()*(a(1) - ba(1)) + 2*q.y()*(a(2) - ba(2)), 2*q.x()*(a(0) - ba(0)) + 2*q.y()*(a(1) - ba(1)) + 2*q.z()*(a(2) - ba(2)), 2*q.x()*(a(1) - ba(1)) - 2*q.y()*(a(0) - ba(0)) + 2*q.w()*(a(2) - ba(2)), 2*q.x()*(a(2) - ba(2)) - 2*q.w()*(a(1) - ba(1)) - 2*q.z()*(a(0) - ba(0)), 0, 0, 0, - q.w()*q.w()- q.x()*q.x()+ q.y()*q.y()+ q.z()*q.z(), 2*q.w()*q.z() - 2*q.x()*q.y(), - 2*q.w()*q.y() - 2*q.x()*q.z(), 0;
	out.row(4) << 0, 0, 0, 0, 0, 0, 2*q.z()*(a(0) - ba(0)) + 2*q.w()*(a(1) - ba(1)) - 2*q.x()*(a(2) - ba(2)), 2*q.y()*(a(0) - ba(0)) - 2*q.x()*(a(1) - ba(1)) - 2*q.w()*(a(2) - ba(2)), 2*q.x()*(a(0) - ba(0)) + 2*q.y()*(a(1) - ba(1)) + 2*q.z()*(a(2) - ba(2)), 2*q.w()*(a(0) - ba(0)) - 2*q.z()*(a(1) - ba(1)) + 2*q.y()*(a(2) - ba(2)), 0, 0, 0, - 2*q.w()*q.z() - 2*q.x()*q.y(), - q.w()*q.w()+ q.x()*q.x()- q.y()*q.y()+ q.z()*q.z(), 2*q.w()*q.x() - 2*q.y()*q.z(), 0;
	out.row(5) << 0, 0, 0, 0, 0, 0, 2*q.x()*(a(1) - ba(1)) - 2*q.y()*(a(0) - ba(0)) + 2*q.w()*(a(2) - ba(2)), 2*q.z()*(a(0) - ba(0)) + 2*q.w()*(a(1) - ba(1)) - 2*q.x()*(a(2) - ba(2)), 2*q.z()*(a(1) - ba(1)) - 2*q.w()*(a(0) - ba(0)) - 2*q.y()*(a(2) - ba(2)), 2*q.x()*(a(0) - ba(0)) + 2*q.y()*(a(1) - ba(1)) + 2*q.z()*(a(2) - ba(2)), 0, 0, 0, 2*q.w()*q.y() - 2*q.x()*q.z(), - 2*q.w()*q.x() - 2*q.y()*q.z(), - q.w()*q.w()+ q.x()*q.x()+ q.y()*q.y()- q.z()*q.z(), 0;

	//dq^dot
	out.row(6) << 0, 0, 0, 0, 0, 0, 0, bg(0)/2 - w(0)/2, bg(1)/2 - w(1)/2, bg(2)/2 - w(2)/2, 0, 0, 0, q.x()/2, q.y()/2, q.z()/2, 0;
	out.row(7) << 0, 0, 0, 0, 0, 0, w(0)/2 - bg(0)/2, 0, w(2)/2 - bg(2)/2, bg(1)/2 - w(1)/2, 0, 0, 0, -q.w()/2, q.z()/2, -q.y()/2, 0;
	out.row(8) << 0, 0, 0, 0, 0, 0, w(1)/2 - bg(1)/2, bg(2)/2 - w(2)/2, 0, w(0)/2 - bg(0)/2, 0, 0, 0, -q.z()/2, -q.w()/2, q.x()/2, 0;
	out.row(9) << 0, 0, 0, 0, 0, 0, w(2)/2 - bg(2)/2, w(1)/2 - bg(1)/2, bg(0)/2 - w(0)/2, 0, 0, 0, 0, q.y()/2, -q.x()/2, -q.w()/2, 0;
}

void TranslationRotationModelV2::get_Q_matrix(const VectorXd &x, const VectorXd &u, const VectorXd &w, MatrixXd& out) {
	Quaternion<double> q = {x(6), x(7), x(8), x(9)};
	Vector3d dv{w(0), w(1), w(2)};
	Vector3d da{w(3), w(4), w(5)};

	MatrixXd W = MatrixXd::Zero(nw + 1, nw + 1);
	W.diagonal() << dv, 0, da;

	MatrixXd G = MatrixXd::Zero(nx, nw + 1);

	MatrixXd temp = MatrixXd::Zero(3, 3);
	rotation_matrix(q, temp);
	G.block<3, 3>(3, 0) = temp;

	temp = MatrixXd::Zero(4, 4);
	hamilton_product_as_matrix(q, temp);
	G.block<4, 4>(6, 3) = temp;

	out = G * W * G.transpose();
}

void TranslationRotationModelV2::get_measurement_estimate(const VectorXd &x, VectorXd& out) {
	out << x(2) + x(16);
}

void TranslationRotationModelV2::get_H_matrix(MatrixXd& out) {
	out << 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
}

void TranslationRotationModelV2::get_R_matrix(MatrixXd& out) {
	out << 0.1;
}

void TranslationRotationModelV2::get_noise_vect(VectorXd& out) {
	double scale_var = 0.5*_dt*_dt;

	out << _accel_params.noise, _accel_params.noise, _accel_params.noise,
			_gyro_params.noise, _gyro_params.noise, _gyro_params.noise,
			_baro_params.noise;

	out *= scale_var;
}

void TranslationRotationModelV2::get_Qs_matrix(MatrixXd& out) {
	double scale_var = 0.5*_dt*_dt;
	double acc_drift = scale_var * _accel_params.drift;
	double gyro_drift = scale_var * _gyro_params.drift;
	double baro_drift = scale_var * _baro_params.drift;

	out.diagonal() <<  additive_noise, additive_noise, additive_noise,
			additive_noise, additive_noise, additive_noise,
			additive_noise, additive_noise, additive_noise, additive_noise,
			acc_drift, acc_drift, acc_drift,
			gyro_drift, gyro_drift, gyro_drift,
			baro_drift;
}