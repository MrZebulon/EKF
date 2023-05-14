//
// Created by Samuel on 02/05/2023.
//

#ifndef SHAI_BASEMODEL_H
#define SHAI_BASEMODEL_H

#include <eigen3/Eigen/Eigen>

namespace shai {
	class EKFEngine;
}

namespace shai::models {
	class BaseModel{
		friend class shai::EKFEngine;

	protected:
		double additive_noise = 1e-8;
		double _dt;
		size_t nx;
		size_t nu;
		size_t nw;
		size_t nz;
	protected:
		explicit BaseModel(double dt, size_t nx, size_t nu, size_t nw, size_t nz)
			: _dt(dt), nx(nx), nu(nu), nw(nw), nz(nz) {}

		virtual Eigen::VectorXd get_init_state() = 0;
		virtual Eigen::MatrixXd get_init_cov() = 0;

		virtual void compute_x_new(const Eigen::VectorXd &x, const Eigen::VectorXd &u, Eigen::VectorXd &out) = 0;
		virtual void get_F_matrix(const Eigen::VectorXd &x, const Eigen::VectorXd &u, Eigen::MatrixXd &out) = 0;
		virtual void get_Q_matrix(const Eigen::VectorXd &x, const Eigen::VectorXd &u, const Eigen::VectorXd &w, Eigen::MatrixXd &out) = 0;

		virtual void get_measurement_estimate(const Eigen::VectorXd &x, Eigen::VectorXd &out) = 0;
		virtual void get_H_matrix(Eigen::MatrixXd &out) = 0;
		virtual void get_R_matrix(Eigen::MatrixXd &out) = 0;

		virtual void get_noise_vect(Eigen::VectorXd &out) = 0;
		virtual void get_Qs_matrix(Eigen::MatrixXd &out) = 0;
	public:
		double get_dt() const {
			return _dt;
		}

		std::size_t get_nx() const {
			return nx;
		}
		std::size_t get_nu() const {
			return nu;
		}
		std::size_t get_nw() const {
			return nw;
		}
		std::size_t get_nz() const {
			return nz;
		}
	};
}
#endif //SHAI_BASEMODEL_H