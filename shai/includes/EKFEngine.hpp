//
// Created by Samuel on 15/05/2023.
//

#ifndef SHAI_EKFENGINE_HPP
#define SHAI_EKFENGINE_HPP

#include <eigen3/Eigen/Dense>
#include <iostream>
#include "models/BaseModel.hpp"
#include "data/EigenDebugger.hpp"

using namespace shai::models;

namespace shai {
	class EKFEngine {
	public:
		explicit EKFEngine(BaseModel* model) : model(model), debugger("../../data/debug_shai.txt") {}
	private:
		BaseModel* model;
		VectorXd x;
		MatrixXd P;
		EigenDebugger debugger;

	public:
		inline void init(){
			x = model->get_init_state();
			P = model->get_init_cov();
		}

		inline void predict(const VectorXd& u) {
			MatrixXd F = model->get_F_matrix(x, u);
			MatrixXd A = MatrixXd::Identity(F.rows(), F.cols()) + F;
			VectorXd w = model->get_noise_vect();
			MatrixXd Qs = model->get_Qs_matrix();

			MatrixXd Q = model->get_Q_matrix(x, u, w);

			x = model->compute_x_new(x, u);
			P = A * P * A.transpose() + Q + Qs;
			debugger.dump_matrix("P (@prior)", P);
		}

		inline void update(const VectorXd& z) {
			MatrixXd H = model->get_H_matrix();
			MatrixXd S = H * P * H.transpose() + model->get_R_matrix();
			MatrixXd K = P * H.transpose() * S.inverse();
			x = x + K*(z - model->get_measurement_estimate(x));
			P = (MatrixXd::Identity(P.rows(), P.cols()) - K*H)*P;
		}

		const VectorXd &get_x() {
			return x;
		}
	};

}

#endif //SHAI_EKFENGINE_HPP
