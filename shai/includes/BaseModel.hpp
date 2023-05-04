//
// Created by Samuel on 02/05/2023.
//

#ifndef SHAI_BASEMODEL_H
#define SHAI_BASEMODEL_H

#include "eigen3/Eigen/Eigen"

namespace shai::models {
	class BaseModel{
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

		virtual Eigen::VectorXd compute_x_new(const Eigen::VectorXd& x, const Eigen::VectorXd& u) = 0;
		virtual Eigen::MatrixXd get_F_matrix(const Eigen::VectorXd& x, const Eigen::VectorXd& u) = 0;
		virtual Eigen::MatrixXd get_Q_matrix(const Eigen::VectorXd& x, const Eigen::VectorXd& u, const Eigen::VectorXd & w) = 0;

		virtual Eigen::VectorXd get_measurement_estimate(const Eigen::VectorXd& x) = 0;
		virtual Eigen::MatrixXd get_H_matrix() = 0;
		virtual Eigen::MatrixXd get_R_matrix() = 0;

		virtual Eigen::VectorXd get_noise_vect() = 0;
		virtual Eigen::MatrixXd get_Qs_matrix() = 0;
	public:
		double get_dt() const {
			return _dt;
		}

		size_t get_nx() const {
			return nx;
		}
		size_t get_nu() const {
			return nu;
		}
		size_t get_nw() const {
			return nw;
		}
		size_t get_nz() const {
			return nz;
		}
	};
}
#endif //SHAI_BASEMODEL_H
