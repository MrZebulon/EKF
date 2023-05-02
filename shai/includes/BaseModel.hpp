//
// Created by Samuel on 02/05/2023.
//

#ifndef SHAI_BASEMODEL_H
#define SHAI_BASEMODEL_H

#include <eigen3/Eigen/Eigen>

namespace shai::models {
	template<std::size_t n_x, std::size_t n_u, std::size_t n_w, std::size_t n_z>
	class BaseModel{
	public:
		typedef Eigen::Vector<double, n_x> xVector;
		typedef Eigen::Vector<double, n_u> uVector;
		typedef Eigen::Vector<double, n_z> zVector;
		typedef Eigen::Vector<double, n_w> wVector;

		typedef Eigen::Matrix<double, n_x, n_x> xxMatrix;
		typedef Eigen::Matrix<double, n_x, n_u> xuMatrix;
		typedef Eigen::Matrix<double, n_u, n_u> uuMatrix;
		typedef Eigen::Matrix<double, n_z, n_x> zxMatrix;
		typedef Eigen::Matrix<double, n_z, n_z> zzMatrix;
	protected:
		double _dt;
	protected:
		explicit BaseModel(double dt) : _dt(dt) {}
		virtual xVector compute_x_new(const xVector& x, const uVector& u) = 0;
		virtual xxMatrix get_F_matrix(const xVector& x, const uVector& u) = 0;
		virtual xxMatrix get_Q_matrix(const xVector& x, const uVector& u, const wVector & w) = 0;

		virtual zVector get_measurement_estimate(const xVector& x) = 0;
		virtual zxMatrix get_H_matrix() = 0;
		virtual zzMatrix get_R_matrix() = 0;

		virtual wVector get_noise_vect() = 0;
		virtual xxMatrix get_Qs_matrix() = 0;
	};
}
#endif //SHAI_BASEMODEL_H
