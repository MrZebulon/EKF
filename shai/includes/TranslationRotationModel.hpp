//
// Created by Samuel on 02/05/2023.
//

#ifndef SHAI_TRANSLATIONROTATIONMODEL_HPP
#define SHAI_TRANSLATIONROTATIONMODEL_HPP

#include "BaseModel.hpp"
#include "SensorParameters.hpp"
namespace shai::models {
	class TranslationRotationModel : public virtual BaseModel{
	public:
		explicit TranslationRotationModel(double dt) : BaseModel(dt, 17, 6, 14, 1) {}
	public:
		SensorParameters<3> _accel_params = {0.0029394957520758797, -0.0009818523846304172, 0.0028761351612281976, 687736578001504e-08, 1.30284870799534e-07};
		SensorParameters<3> _gyro_params = {0.003828461092875141, -0.001784311643174541, -0.002920534852288187, 5.047998984990613e-07, 1.0104519055951595e-06};
		SensorParameters<1> _baro_params = {399.23599189364484, 0.014769875002625985, 6.28236143573932e-05};
	protected:

		Eigen::VectorXd get_init_state() override;
		Eigen::MatrixXd get_init_cov() override;

		void compute_x_new(const Eigen::VectorXd &x, const Eigen::VectorXd &u, Eigen::VectorXd &out) override;
		void get_F_matrix(const Eigen::VectorXd &x, const Eigen::VectorXd &u, Eigen::MatrixXd &out) override;
		void get_Q_matrix(const Eigen::VectorXd &x, const Eigen::VectorXd &u, const Eigen::VectorXd &w, Eigen::MatrixXd &out) override;

		void get_measurement_estimate(const Eigen::VectorXd &x, Eigen::VectorXd &out) override;
		void get_H_matrix(Eigen::MatrixXd &out) override;
		void get_R_matrix(Eigen::MatrixXd &out) override;

		void get_noise_vect(Eigen::VectorXd &out) override;
		void get_Qs_matrix(Eigen::MatrixXd &out) override;
	};
}

#endif //SHAI_TRANSLATIONROTATIONMODEL_HPP
