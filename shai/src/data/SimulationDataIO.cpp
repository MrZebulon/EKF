//
// Created by Samuel on 13/05/2023.
//

#include <iostream>
#include <cmath>
#include "data/SimulationDataIO.hpp"

bool SimulationDataIO::read_next(Eigen::VectorXd &out) {
	double v = 0.;

	if(!std::getline(file_in, _in_str))
		return false;

	_in_stream.str(_in_str);

	int i = 0;
	while (std::getline(_in_stream, _in_str, ',')) {
		_in_stream >> v;
		out(i++) = v;
	}

	return true;
}

void SimulationDataIO::write_next(const Eigen::VectorXd &in) {

	_out_stream.str(std::string());

	for (const double &v: in)
		_out_stream << v << ',';

	_out_stream.seekp(-1, std::ios_base::end); // Remove last char : ','
	_out_stream << std::endl;

	file_out << _out_stream.str();
}
