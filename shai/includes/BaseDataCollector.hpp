//
// Created by Samuel on 10/05/2023.
//

#ifndef SHAI_BASEDATACOLLECTOR_HPP
#define SHAI_BASEDATACOLLECTOR_HPP

#include <eigen3/Eigen/Core>

namespace shai::data {
	class BaseDataCollector {
		virtual void collect_data(Eigen::VectorXd& out) = 0;
	};
}


#endif //SHAI_BASEDATACOLLECTOR_HPP
