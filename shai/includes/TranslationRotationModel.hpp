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
	private:
		SensorParameters<3> _accel_params = {-0.026842656242568, 0.033420780321046, -0.007947030636161,1, 2e-4};
		SensorParameters<3> _gyro_params = {5.735129093702078e-05,-5.361232250514282e-06,3.135626159376753e-05,1e-5, 0};
		SensorParameters<1> _baro_params = {361.3487972164834, 0.0007454259701653068, 2.8486463440220755e-06};
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
