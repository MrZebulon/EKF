//
// Created by Samuel on 02/05/2023.
//

#include "../includes/TranslationRotationModel.hpp"

using namespace shai::models::TRM;

xVector TranslationRotationModel::compute_x_new(const xVector &x, const uVector &u) {

}

xxMatrix TranslationRotationModel::get_F_matrix(const xVector &x, const uVector &u) {

}

xxMatrix TranslationRotationModel::get_Q_matrix(const xVector &x, const uVector &u, const wVector &w) {

}

zVector TranslationRotationModel::get_measurement_estimate(const xVector &x) {

}

zxMatrix TranslationRotationModel::get_H_matrix() {

}

zzMatrix TranslationRotationModel::get_R_matrix() {

}

wVector TranslationRotationModel::get_noise_vect() {

}

xxMatrix TranslationRotationModel::get_Qs_matrix() {

}
