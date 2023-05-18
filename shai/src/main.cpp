#include "data/SimulationDataIO.hpp"
#include "EKFEngine.hpp"
#include "models/TranslationRotationModel.hpp"

using namespace shai;

int main() {
	SimulationDataIO io("../../../Data/static/raw_shai.csv", "../../../Data/static/run_shai.csv");

	models::TranslationRotationModel model(1./100);
	EKFEngine engine(&model);
	engine.init();

	Eigen::VectorXd data_point = Eigen::VectorXd::Zero(7);

	while(io.read_next(data_point)){
		engine.predict(data_point.segment<6>(1));
		engine.update(data_point.segment<1>(0));
		io.write_next(engine.get_x());
	}
}
