//
// Created by Samuel on 15/05/2023.
//

#ifndef SHAI_TRANSLATIONROTATIONMODEL_HPP
#define SHAI_TRANSLATIONROTATIONMODEL_HPP

#include "BaseModel.hpp"
#include "SensorParameters.hpp"

namespace shai::models {

	class TranslationRotationModel : public virtual BaseModel {
	private:
		SensorParameters<3> accel = {{0.0029394893123127377, -0.0009818539510592773, 0.0028762637247315066}, 4.7358607479753485e-09, 3.314312818032142e-10};
		SensorParameters<3> gyro = {{0.0038284331173227743, -0.001784190927555866, -0.002920534852288187}, 1.0102261028815603e-08, 3.9979150848914756e-10};
		SensorParameters<1> baro = {Array<double, 1, 1>(399.23657624056926), 0.014769875002697693, 6.282361435771973e-05};
	public:
		explicit TranslationRotationModel() : BaseModel(1./100) {}

		VectorXd get_init_state() override;
		MatrixXd get_init_cov() override;

		VectorXd compute_x_new(const VectorXd &x, const VectorXd &u) override;
		MatrixXd get_F_matrix(const VectorXd &x, const VectorXd &u) override;
		MatrixXd get_Q_matrix(const VectorXd &x, const VectorXd &u, const VectorXd &w) override;

		VectorXd get_measurement_estimate(const VectorXd &x) override;
		MatrixXd get_H_matrix() override;
		MatrixXd get_R_matrix() override;

		VectorXd get_noise_vect() override;
		MatrixXd get_Qs_matrix() override;
	};
}

#endif //SHAI_TRANSLATIONROTATIONMODEL_HPP
