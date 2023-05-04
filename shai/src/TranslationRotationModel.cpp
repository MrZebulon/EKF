//
// Created by Samuel on 02/05/2023.
//

#include "../includes/TranslationRotationModel.hpp"

using namespace shai::models;
using namespace Eigen;

void TranslationRotationModel::compute_x_new(const VectorXd  &x, const VectorXd &u, VectorXd& out) {
	out = VectorXd::Zero(nx);
}

void TranslationRotationModel::get_F_matrix(const VectorXd &x, const VectorXd  &u, MatrixXd& out) {
	out = MatrixXd::Zero(nx, nx);
}

void TranslationRotationModel::get_Q_matrix(const VectorXd &x, const VectorXd &u, const VectorXd &w, MatrixXd& out) {
	out = MatrixXd::Zero(nx, nx);
}

void  TranslationRotationModel::get_measurement_estimate(const VectorXd &x, VectorXd& out) {
	out = VectorXd::Zero(nz);
}

void TranslationRotationModel::get_H_matrix(MatrixXd& out) {
	out = MatrixXd::Zero(nz, nx);
}

void TranslationRotationModel::get_R_matrix(MatrixXd& out) {
	out = MatrixXd::Zero(nz, nz);
}

void TranslationRotationModel::get_noise_vect(VectorXd& out) {
	out = VectorXd::Zero(nw);
}

void TranslationRotationModel::get_Qs_matrix(MatrixXd& out) {
	out = MatrixXd::Zero(nx, nx);
}
