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
		explicit TranslationRotationModel(double dt) : BaseModel(dt, 17, 6, 7, 1) {}
	public:
		SensorParameters<3> _accel_params = {0.0029394893123127377, -0.0009818539510592773, 0.0028762637247315066, 4.7358607479753485e-09, 3.314312818032142e-10};
		SensorParameters<3> _gyro_params = {0.0038284331173227743, -0.001784190927555866, -0.0028707138021085243, 2.5617707075843632e-08, 0};
		SensorParameters<1> _baro_params = {399.23657624056926, 0.0007454259701653068, 2.8486463440220755e-06};
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
