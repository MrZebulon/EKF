//
// Created by Samuel on 17/05/2023.
//

#ifndef SHAI_EIGENDEBUGGER_HPP
#define SHAI_EIGENDEBUGGER_HPP
#include <fstream>
#include <eigen3/Eigen/Dense>

class EigenDebugger {
private:
	std::ofstream file_out;
public:
	explicit EigenDebugger(const char* debug_path){
		file_out = std::ofstream(debug_path, std::ofstream::out | std::ofstream::trunc);
		file_out.close();
		file_out = std::ofstream(debug_path, std::ios_base::app);
	}

	inline void dump_matrix(const char* name, Eigen::MatrixXd& m){
		file_out << "Matrix " << name << ": \n" << m << '\n';
	}

	inline void dump_vector(const char* name, Eigen::VectorXd& v){
		file_out << "Vector " << name << ": \n"<< v << '\n';
	}
};

#endif //SHAI_EIGENDEBUGGER_HPP
