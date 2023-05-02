//
// Created by Samuel on 02/05/2023.
//

#ifndef SHAI_TRANSLATIONROTATIONMODEL_HPP
#define SHAI_TRANSLATIONROTATIONMODEL_HPP

#include "BaseModel.hpp"

namespace shai::models::TRM {
	class TranslationRotationModel : public virtual BaseModel<17, 6, 14, 1>{
	protected:
		xVector compute_x_new(const xVector &x, const uVector &u) override;
		xxMatrix get_F_matrix(const xVector &x, const uVector &u) override;
		xxMatrix get_Q_matrix(const xVector &x, const uVector &u, const wVector &w) override;

		zVector get_measurement_estimate(const xVector &x) override;
		zxMatrix get_H_matrix() override;
		zzMatrix get_R_matrix() override;

		wVector get_noise_vect() override;
		xxMatrix get_Qs_matrix() override;
	};

	using xVector = TranslationRotationModel::xVector;
	using uVector = TranslationRotationModel::uVector;
	using zVector = TranslationRotationModel::zVector;
	using wVector = TranslationRotationModel::wVector;
	using xxMatrix = TranslationRotationModel::xxMatrix;
	using xuMatrix = TranslationRotationModel::xuMatrix;
	using uuMatrix = TranslationRotationModel::uuMatrix;
	using zxMatrix = TranslationRotationModel::zxMatrix;
	using zzMatrix = TranslationRotationModel::zzMatrix;
}

#endif //SHAI_TRANSLATIONROTATIONMODEL_HPP
