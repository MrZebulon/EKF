//
// Created by Samuel on 02/05/2023.
//

#include <cmath>

#include "../includes/TranslationRotationModel.hpp"
#include "../includes/Utils.hpp"

using namespace std;
using namespace shai::models;
using namespace Eigen;

Eigen::VectorXd TranslationRotationModel::get_init_state() {
	Eigen::Vector3d ba {_accel_params.bias[0], _accel_params.bias[1], _accel_params.bias[2]};
	Eigen::Vector3d bg {_gyro_params.bias[0], _gyro_params.bias[1], _gyro_params.bias[2]};
	double bb = _baro_params.bias[0];

	Eigen::VectorXd x_init(nx);
	x_init << Eigen::VectorXd::Zero(10), ba, bg, bb;

	return x_init;
}

Eigen::MatrixXd TranslationRotationModel::get_init_cov() {
	Eigen::MatrixXd P_init(nx, nx);
	P_init.diagonal() << Eigen::VectorXd::Constant(nx, 1e-9);
	return P_init;
}

void TranslationRotationModel::compute_x_new(const VectorXd &x, const VectorXd &u, VectorXd& out) {

	Eigen::Vector3d pos = {x(0), x(1), x(2)};
	Eigen::Vector3d vel = {x(3), x(4), x(5)};
	Eigen::Quaternion<double> q = {x(6), x(7), x(8), x(9)};
	Eigen::Vector3d ba {x(10), x(11), x(12)};
	Eigen::Vector3d bg{x(13), x(14), x(15)};
	double bb = x(16);

	Eigen::Vector3d a = {u(0), u(1), u(2)};
	Eigen::Vector3d w = {u(3), u(4), u(5)};

	Eigen::Vector3d delta_v{};

	body_to_inertial({q.w(), q.x(), q.y(), q.z()}, {_dt*(a(0) - ba(0)), _dt*(a(1) - ba(1)), _dt*(a(2) - ba(2))}, delta_v);
	Eigen::Quaternion<double> delta_q = {1,  0.5*_dt*(w(0) - bg(0)),  0.5*_dt*(w(1) - bg(1)),  0.5*_dt*(w(2) - bg(2))};

	q *= delta_q;

	out << pos + vel*_dt, vel + delta_v, q.coeffs(), ba, bg, bb;
}

void TranslationRotationModel::get_F_matrix(const VectorXd &x, const VectorXd  &u, MatrixXd& out) {
	Eigen::Quaternion<double> q = {x(6), x(7), x(8), x(9)};
	Eigen::Vector3d ba {x(10) * _dt, x(11) * _dt, x(12) * _dt};
	Eigen::Vector3d bg {x(13) * _dt, x(14) * _dt, x(15) * _dt};

	Eigen::Vector3d a = {u(0) * _dt, u(1) * _dt, u(2) * _dt};
	Eigen::Vector3d w = {u(3) * _dt, u(4) * _dt, u(5) * _dt};

	// dp^dot
	out(0, 3) = _dt;
	out(1, 4) = _dt;
	out(2, 5) = _dt;

	// dv^dot
	out.row(3) << 0, 0, 0, 0, 0, 0, 2*q.w()*(a(0) - ba(0)) - 2*q.z()*(a(1) - ba(1)) + 2*q.y()*(a(2) - ba(2)), 2*q.x()*(a(0) - ba(0)) + 2*q.y()*(a(1) - ba(1)) + 2*q.z()*(a(2) - ba(2)), 2*q.x()*(a(1) - ba(1)) - 2*q.y()*(a(0) - ba(0)) + 2*q.w()*(a(2) - ba(2)), 2*q.x()*(a(2) - ba(2)) - 2*q.w()*(a(1) - ba(1)) - 2*q.z()*(a(0) - ba(0)), 0, 0, 0, - q.w()*q.w() - q.x()*q.x() + q.y()*q.y() + q.z()*q.z(), 2*q.w()*q.z() - 2*q.x()*q.y(), - 2*q.w()*q.y() - 2*q.x()*q.z(), 0;
	out.row(4) << 0, 0, 0, 0, 0, 0, 2*q.z()*(a(0) - ba(0)) + 2*q.w()*(a(1) - ba(1)) - 2*q.x()*(a(2) - ba(2)), 2*q.y()*(a(0) - ba(0)) - 2*q.x()*(a(1) - ba(1)) - 2*q.w()*(a(2) - ba(2)), 2*q.x()*(a(0) - ba(0)) + 2*q.y()*(a(1) - ba(1)) + 2*q.z()*(a(2) - ba(2)), 2*q.w()*(a(0) - ba(0)) - 2*q.z()*(a(1) - ba(1)) + 2*q.y()*(a(2) - ba(2)), 0, 0, 0, - 2*q.w()*q.z() - 2*q.x()*q.y(), - q.w()*q.w() + q.x()*q.x() - q.y()*q.y() + q.z()*q.z(), 2*q.w()*q.x() - 2*q.y()*q.z(), 0;
	out.row(5) << 0, 0, 0, 0, 0, 0, 2*q.x()*(a(1) - ba(1)) - 2*q.y()*(a(0) - ba(0)) + 2*q.w()*(a(2) - ba(2)), 2*q.z()*(a(0) - ba(0)) + 2*q.w()*(a(1) - ba(1)) - 2*q.x()*(a(2) - ba(2)), 2*q.z()*(a(1) - ba(1)) - 2*q.w()*(a(0) - ba(0)) - 2*q.y()*(a(2) - ba(2)), 2*q.x()*(a(0) - ba(0)) + 2*q.y()*(a(1) - ba(1)) + 2*q.z()*(a(2) - ba(2)), 0, 0, 0, 2*q.w()*q.y() - 2*q.x()*q.z(), - 2*q.w()*q.x() - 2*q.y()*q.z(), - q.w()*q.w() + q.x()*q.x() + q.y()*q.y() - q.z()*q.z(), 0;

	//dq^dot
	out.row(6) << 0, 0, 0, 0, 0, 0, 0,  bg(0)/2 - w(0)/2,  bg(1)/2 - w(1)/2,  bg(2)/2 - w(2)/2, 0, 0, 0, q.x()/2, q.y()/2, q.z()/2, 0;
	out.row(7) << 0, 0, 0, 0, 0, 0, w(0)/2 - bg(0)/2, 0, w(2)/2 - bg(2)/2,  bg(1)/2 - w(1)/2, 0, 0, 0, -q.w()/2, q.z()/2, -q.y()/2, 0;
	out.row(8) << 0, 0, 0, 0, 0, 0, w(1)/2 - bg(1)/2,  bg(2)/2 - w(2)/2, 0, w(0)/2 - bg(0)/2, 0, 0, 0, -q.z()/2, -q.w()/2, q.x()/2, 0;
	out.row(9) << 0, 0, 0, 0, 0, 0, w(2)/2 - bg(2)/2, w(1)/2 - bg(1)/2,  bg(0)/2 - w(0)/2, 0, 0, 0, 0, q.y()/2, -q.x()/2, -q.w()/2, 0;
}

void TranslationRotationModel::get_Q_matrix(const VectorXd &x, const VectorXd &u, const VectorXd &w, MatrixXd& out) {
	Eigen::Vector3d dvCov = {w(0), w(1), w(2)};
	Eigen::Vector3d daCov = {w(3), w(4), w(5)};

	Eigen::Quaternion<double> q = {x(6), x(7), x(8), x(9)};

	out.row(3) << 0, 0, 0, dvCov(1)* pow(2*q.w()*q.z() - 2*q.x()*q.y(), 2) + dvCov(2)*pow(2*q.w()*q.y() + 2*q.x()*q.z(), 2) + dvCov(0)*pow(q.w()*q.w() + q.x()*q.x() - q.y()*q.y() - q.z()*q.z(), 2), dvCov(0)*(2*q.w()*q.z() + 2*q.x()*q.y())*(q.w()*q.w() + q.x()*q.x() - q.y()*q.y() - q.z()*q.z()) - dvCov(1)*(2*q.w()*q.z() - 2*q.x()*q.y())*(q.w()*q.w() - q.x()*q.x() + q.y()*q.y() - q.z()*q.z()) - dvCov(2)*(2*q.w()*q.x() - 2*q.y()*q.z())*(2*q.w()*q.y() + 2*q.x()*q.z()), dvCov(2)*(2*q.w()*q.y() + 2*q.x()*q.z())*(q.w()*q.w() - q.x()*q.x() - q.y()*q.y() + q.z()*q.z()) - dvCov(0)*(2*q.w()*q.y() - 2*q.x()*q.z())*(q.w()*q.w() + q.x()*q.x() - q.y()*q.y() - q.z()*q.z()) - dvCov(1)*(2*q.w()*q.x() + 2*q.y()*q.z())*(2*q.w()*q.z() - 2*q.x()*q.y()), 0, 0, 0, 0,  0, 0, 0,  0, 0, 0, 0;
	out.row(4) << 0, 0, 0, dvCov(0) * (2*q.w()*q.z() + 2*q.x()*q.y())*(q.w()*q.w() + q.x()*q.x() - q.y()*q.y() - q.z()*q.z()) - dvCov(1) *(2*q.w()*q.z() - 2*q.x()*q.y())*(q.w()*q.w() - q.x()*q.x() + q.y()*q.y() - q.z()*q.z()) - dvCov(2) *(2*q.w()*q.x() - 2*q.y()*q.z())*(2*q.w()*q.y() + 2*q.x()*q.z()), dvCov(0) * pow(2*q.w()*q.z() + 2*q.x()*q.y(), 2) + dvCov(2) * pow(2*q.w()*q.x() - 2*q.y()*q.z(), 2) + dvCov(1) * pow(q.w()*q.w() - q.x()*q.x() + q.y()*q.y() - q.z()*q.z(), 2), dvCov(1) *(2*q.w()*q.x() + 2*q.y()*q.z())*(q.w()*q.w() - q.x()*q.x() + q.y()*q.y() - q.z()*q.z()) - dvCov(2) *(2*q.w()*q.x() - 2*q.y()*q.z())*(q.w()*q.w() - q.x()*q.x() - q.y()*q.y() + q.z()*q.z()) - dvCov(0) *(2*q.w()*q.y() - 2*q.x()*q.z())*(2*q.w()*q.z() + 2*q.x()*q.y()), 0, 0, 0, 0,  0, 0, 0,  0, 0, 0, 0;
	out.row(5) << 0, 0, 0, dvCov(2) * (2*q.w()*q.y() + 2*q.x()*q.z())*(q.w()*q.w() - q.x()*q.x() - q.y()*q.y() + q.z()*q.z()) - dvCov(0) *(2*q.w()*q.y() - 2*q.x()*q.z())*(q.w()*q.w() + q.x()*q.x() - q.y()*q.y() - q.z()*q.z()) - dvCov(1) *(2*q.w()*q.x() + 2*q.y()*q.z())*(2*q.w()*q.z() - 2*q.x()*q.y()), dvCov(1) *(2*q.w()*q.x() + 2*q.y()*q.z())*(q.w()*q.w() - q.x()*q.x() + q.y()*q.y() - q.z()*q.z()) - dvCov(2) *(2*q.w()*q.x() - 2*q.y()*q.z())*(q.w()*q.w() - q.x()*q.x() - q.y()*q.y() + q.z()*q.z()) - dvCov(0) *(2*q.w()*q.y() - 2*q.x()*q.z())*(2*q.w()*q.z() + 2*q.x()*q.y()), dvCov(0) * pow(2*q.w()*q.y() - 2*q.x()*q.z(), 2) + dvCov(1) * pow(2*q.w()*q.x() + 2*q.y()*q.z(), 2) + dvCov(2) * pow(q.w()*q.w() - q.x()*q.x() - q.y()*q.y() + q.z()*q.z(), 2), 0, 0, 0, 0,  0, 0, 0,  0, 0, 0, 0;

	out.row(6) << 0, 0, 0,  0, 0, 0, (daCov(0) *q.x()*q.x())/4 + (daCov(1) *q.y()*q.y())/4 + (daCov(2) *q.z()*q.z())/4, (daCov(1) *q.y()*q.z())/4 - (daCov(0) *q.w()*q.x())/4 - (daCov(2) *q.y()*q.z())/4, (daCov(2) *q.x()*q.z())/4 - (daCov(1) *q.w()*q.y())/4 - (daCov(0) *q.x()*q.z())/4, (daCov(0) *q.x()*q.y())/4 - (daCov(1) *q.x()*q.y())/4 - (daCov(2) *q.w()*q.z())/4,  0, 0, 0,   0, 0, 0,  0;
	out.row(7) << 0, 0, 0,  0, 0, 0, (daCov(1) *q.y()*q.z())/4 - (daCov(0) *q.w()*q.x())/4 - (daCov(2) *q.y()*q.z())/4, (daCov(0) *q.w()*q.w())/4 + (daCov(2) *q.y()*q.y())/4 + (daCov(1) *q.z()*q.z())/4, (daCov(0) *q.w()*q.z())/4 - (daCov(1) *q.w()*q.z())/4 - (daCov(2) *q.x()*q.y())/4, (daCov(2) *q.w()*q.y())/4 - (daCov(1) *q.x()*q.z())/4 - (daCov(0) *q.w()*q.y())/4,  0, 0, 0,   0, 0, 0,  0;
	out.row(8) << 0, 0, 0,  0, 0, 0, (daCov(2) *q.x()*q.z())/4 - (daCov(1) *q.w()*q.y())/4 - (daCov(0) *q.x()*q.z())/4, (daCov(0) *q.w()*q.z())/4 - (daCov(1) *q.w()*q.z())/4 - (daCov(2) *q.x()*q.y())/4, (daCov(1) *q.w()*q.w())/4 + (daCov(2) *q.x()*q.x())/4 + (daCov(0) *q.z()*q.z())/4, (daCov(1) *q.w()*q.x())/4 - (daCov(0) *q.y()*q.z())/4 - (daCov(2) *q.w()*q.x())/4,  0, 0, 0,   0, 0, 0,  0;
	out.row(9) << 0, 0, 0, 0, 0, 0, (daCov(0) *q.x()*q.y())/4 - (daCov(1) *q.x()*q.y())/4 - (daCov(2) *q.w()*q.z())/4, (daCov(2) *q.w()*q.y())/4 - (daCov(1) *q.x()*q.z())/4 - (daCov(0) *q.w()*q.y())/4, (daCov(1) *q.w()*q.x())/4 - (daCov(0) *q.y()*q.z())/4 - (daCov(2) *q.w()*q.x())/4,(daCov(2) *q.w()*q.w())/4 + (daCov(1) *q.x()*q.x())/4 + (daCov(0) *q.y()*q.y())/4,   0, 0, 0,   0, 0, 0,  0;
}

void  TranslationRotationModel::get_measurement_estimate(const VectorXd &x, VectorXd& out) {
	out << x(2) + x(16);
}

void TranslationRotationModel::get_H_matrix(MatrixXd& out) {
	out << 0, 0, 1,  0, 0, 0,  0, 0, 0, 0,  0, 0, 0,  0, 0, 0,  0;
}

void TranslationRotationModel::get_R_matrix(MatrixXd& out) {
	out << 0.1;
}

void TranslationRotationModel::get_noise_vect(VectorXd& out) {
	double scale_var = 0.5*_dt*_dt;

	out <<  _accel_params.noise, _accel_params.noise, _accel_params.noise,
			_gyro_params.noise, _gyro_params.noise, _gyro_params.noise,
			_baro_params.noise;

	out *= scale_var;
}

void TranslationRotationModel::get_Qs_matrix(MatrixXd& out) {
	double scale_var = 0.5*_dt*_dt;
	double acc_drift = scale_var * _accel_params.drift;
	double gyro_drift = scale_var * _gyro_params.drift;
	double baro_drift = scale_var * _baro_params.drift;

	out.diagonal() <<   additive_noise, additive_noise, additive_noise,
						additive_noise, additive_noise, additive_noise,
						additive_noise, additive_noise, additive_noise, additive_noise,
						acc_drift, acc_drift, acc_drift,
						gyro_drift, gyro_drift, gyro_drift,
						baro_drift;
}