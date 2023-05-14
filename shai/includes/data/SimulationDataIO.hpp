//
// Created by Samuel on 13/05/2023.
//

#ifndef SHAI_SIMULATIONDATAIO_HPP
#define SHAI_SIMULATIONDATAIO_HPP

#include <sstream>
#include <fstream>
#include <cstring>
#include <eigen3/Eigen/Core>

class SimulationDataIO {
private:
	std::ifstream file_in;
	std::stringstream _in_stream;
	std::string _in_str;

	std::ofstream file_out;
	std::stringstream _out_stream;
public:
	explicit SimulationDataIO(const char* in_path, const char* out_path){
		file_in = std::ifstream(in_path);
		file_out = std::ofstream(out_path, std::ofstream::out | std::ofstream::trunc);
		file_out.close();
		file_out = std::ofstream(out_path, std::ios_base::app);
	}

	bool read_next(Eigen::VectorXd& out);
	void write_next(const Eigen::VectorXd& in);
};


#endif //SHAI_SIMULATIONDATAIO_HPP
