//
// Created by Samuel on 02/05/2023.
//

#ifndef SHAI_EKFENGINE_HPP
#define SHAI_EKFENGINE_HPP

#include "BaseModel.hpp"

namespace shai{

	struct EKFEngineParameters {
		Eigen::VectorXd x;
		Eigen::MatrixXd P;

		Eigen::MatrixXd F;
		Eigen::MatrixXd Q;

		Eigen::VectorXd w;
		Eigen::MatrixXd Qs;

		Eigen::VectorXd h;
		Eigen::MatrixXd H;
		Eigen::MatrixXd R;

		Eigen::VectorXd inov;
		Eigen::MatrixXd S;
		Eigen::MatrixXd K;
	};

	template<class Model = models::BaseModel>
	class EKFEngine {
	private:

	public:
		explicit EKFEngine(const Model& model) : _model(model) {}

		void init(){
			_params.x = _model.init_state();
			_params.P = _model.init_cov();
		}

		void run_once(const Eigen::VectorXd& u, const Eigen::VectorXd& z){
			predict(u);
			update(z);
			clear();
		}

	private:
		Model _model;
		EKFEngineParameters _params;
	protected:
		void predict(const Eigen::VectorXd& u) {
			_model.get_F_matrix(_params.x, u, _params.F);
			_params.F = Eigen::MatrixXd::Identity(_model._nx, _model._nx) + _params.F;
			_model.get_Qs_matrix(_params.Qs);
			_model.get_noise_vect(_params.w);
			_model.get_Q_matrix(_params.x, u, _params.w);
			_params.P = _params.F * _params.P * _params.F.transpose() + _params.Q + _params.Qs;

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

		void clear() {
			_params.F.setZero();
			_params.Q.setZero();

			_params.w.setZero();
			_params.Qs.setZero();

			_params.h.setZero();
			_params.H.setZero();
			_params.R.setZero();

			_params.inov.setZero();
			_params.S.setZero();
			_params.K.setZero();
		}

		const EKFEngineParameters& get_params() const {
			return _params;
		}
	};
};

#endif //SHAI_EKFENGINE_HPP
