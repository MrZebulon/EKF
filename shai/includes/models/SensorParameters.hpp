//
// Created by Samuel on 04/05/2023.
//

#ifndef SHAI_SENSORPARAMETERS_HPP
#define SHAI_SENSORPARAMETERS_HPP

#include <eigen3/Eigen/Core>

namespace shai::models{
	template<std::size_t dims>
	struct SensorParameters{
		Eigen::Array<double, dims, 1> bias;
		double noise;
		double drift;
	};
}

#endif //SHAI_SENSORPARAMETERS_HPP
