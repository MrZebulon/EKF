//
// Created by Samuel on 04/05/2023.
//

#ifndef SHAI_SENSORPARAMETERS_HPP
#define SHAI_SENSORPARAMETERS_HPP

#include <cmath>

namespace shai::models{
	template<std::size_t dims>
	struct SensorParameters{
		double bias[dims];
		double noise;
		double drift;
	};
}

#endif //SHAI_SENSORPARAMETERS_HPP
