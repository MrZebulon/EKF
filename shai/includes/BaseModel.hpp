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
		typedef Eigen::Vector<double, n_x> StateVector;
		typedef Eigen::Vector<double, n_u> InputVector;
		typedef Eigen::Vector<double, n_z> ObservationVector;
		typedef Eigen::Vector<double, n_w> NoiseVector;
	protected:
		double _dt;
	public:
		explicit BaseModel(double dt) : _dt(dt) {}
		virtual StateVector compute_x_new(const StateVector& x, const InputVector& u) = 0;
		virtual Eigen::Matrix<double, n_x, n_x> get_F_matrix(const StateVector& x, const InputVector& u) = 0;
		virtual Eigen::Matrix<double, n_x, n_x> get_Q_matrix(const StateVector& x, const InputVector& u, const NoiseVector & w) = 0;

		virtual ObservationVector get_measurement_estimate(const StateVector& x) = 0;
		virtual Eigen::Matrix<double, n_z, n_z> get_H_matrix() = 0;

		virtual NoiseVector get_noise_vect() = 0;
		virtual Eigen::Matrix<double, n_x, n_x> get_Qs_matrix() = 0;
	};
}
#endif //SHAI_BASEMODEL_H
