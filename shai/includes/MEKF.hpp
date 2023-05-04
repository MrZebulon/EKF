//
// Created by Samuel on 02/05/2023.
//

#ifndef SHAI_MEKF_HPP
#define SHAI_MEKF_HPP

#include "BaseModel.hpp"

namespace shai{

	struct MEKFParameters {
		Eigen::VectorXd x;
		Eigen::MatrixXd P;

		Eigen::VectorXd w;

		Eigen::MatrixXd F;
		Eigen::MatrixXd A;
		Eigen::MatrixXd Q;
		Eigen::MatrixXd Qs;

		Eigen::VectorXd h;
		Eigen::MatrixXd H;
		Eigen::MatrixXd R;

		Eigen::VectorXd inov;
		Eigen::MatrixXd S;
		Eigen::MatrixXd K;
	};

	template<class Model = models::BaseModel>
	class MEKF {
	private:

	public:
		explicit MEKF(const Model& model) : _model(model) {}

	private:
		Model _model;
		MEKFParameters _params;
	public:

		void init(){
			_params.x = _model.init_state();
			_params.P = _model.init_cov();
		}

		void predict(const Eigen::VectorXd& u) {
			_model.get_F_matrix(_params.x, u, _params.F);
			_params.A = Eigen::MatrixXd::Identity(_model._nx, _model._nx) + _params.F;
			_model.get_Qs_matrix(_params.Qs);
			_model.get_noise_vect(_params.w);
			_model.get_Q_matrix(_params.x, u, _params.w);
			_params.P = _params.A * _params.P * _params.A.transpose() + _params.Q + _params.Qs;

			_model.compute_x_new(_params.x, u, _params.x);
			_params.P = 0.5 * (_params.P + _params.P.transpose());
		}

		void update(const Eigen::VectorXd& z) {
			_model.get_H_matrix(_params.H);
			_params.inov = z - _model.get_measurement_estimate(_params.x);
			_model.get_R_matrix(_params.R);
			_params.S = _params.H * _params.P * _params.H.transpose() + _params.R;
			_params.K = _params.P * (_params.H.transpose()) * _params.S.inverse();
			_params.x = _params.x + _params.K * _params.inov;
			_params.P = (Eigen::MatrixXd::Identity(_model._nx, _model._nx) - _params.K * _params.H) * _params.P;
		}

		const MEKFParameters &get_state() const {
			return _params;
	};
};

#endif //SHAI_MEKF_HPP
