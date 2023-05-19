//
// Created by Samuel on 02/05/2023.
//

#ifndef SHAI_BASEMODEL_H
#define SHAI_BASEMODEL_H

#include <eigen3/Eigen/Core>

using namespace Eigen;

namespace shai::models {

	class BaseModel{
	protected:
		double additive_noise = 1e-8;
		double _dt;

		explicit BaseModel(double dt) : _dt(dt) {}

	public:
		virtual VectorXd get_init_state() = 0;
		virtual MatrixXd get_init_cov() = 0;

		virtual VectorXd compute_x_new(const VectorXd &x, const VectorXd& u) = 0;
		virtual MatrixXd get_F_matrix(const VectorXd &x, const VectorXd& u) = 0;
		virtual MatrixXd get_Q_matrix(const VectorXd &x, const VectorXd& u, const VectorXd &w) = 0;

		virtual VectorXd get_measurement_estimate(const VectorXd &x) = 0;
		virtual MatrixXd get_H_matrix() = 0;
		virtual MatrixXd get_R_matrix() = 0;

		virtual VectorXd get_noise_vect() = 0;
		virtual MatrixXd get_Qs_matrix() = 0;
	};
}
#endif //SHAI_BASEMODEL_H
