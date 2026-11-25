/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2026  Laurent Forthomme
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include <algorithm>  // std::clamp
#include <cmath>
#include <iostream>

#include "CepGen/HelicityAmplitudes/Kinematics.h"

using namespace cepgen::amplitudes;

Kinematics::Kinematics(double mv2, double s_hat, double t_hat, double u_hat)
    : Kinematics(mv2, mv2, s_hat, t_hat, u_hat) {}

Kinematics::Kinematics(double mv32, double mv42, double s_hat, double t_hat, double u_hat)
    : shat(s_hat),
      that(t_hat),
      uhat(u_hat),
      mv32_(mv32),
      mv42_(mv42),
      shat2(s_hat * s_hat),
      // sqrt(lambda(s,m3^2,m4^2))/s, reducing to 1-4m^2/s squared for identical masses
      beta2(((s_hat - mv32_ - mv42_) * (s_hat - mv32_ - mv42_) - 4. * mv32_ * mv42_) / (s_hat * s_hat)),
      beta(sqrt(beta2)),
      inv_gamma2(1. - beta2),
      gamma2(1. / inv_gamma2),
      gamma(sqrt(gamma2)),
      inv_gamma(1. / gamma) {
  setCosTheta((t_hat - u_hat) / s_hat / beta);
}

Kinematics Kinematics::fromSCosTheta(double shat, double cos_theta, double mv2) {
  Kinematics kin(mv2, shat, 0., 0.);
  kin.setCosTheta(cos_theta);
  return kin;
}

void Kinematics::setCosTheta(double cth) {
  // with off-shell photons, (t-u)/(s*beta) may exceed unity, possibly reaching the 1/A pole at beta*|cos(theta)|=1
  cos_theta = std::clamp(cth, -1., 1.);
  cos_theta2 = cos_theta * cos_theta;
  sin_theta2 = 1. - cos_theta2;
  sin_theta = sqrt(fabs(sin_theta2));
  invA = 1. / (1. - beta2 * cos_theta2);
}

bool Kinematics::operator!=(const Kinematics& oth) const {
  // checking only the base variables as all others are computed from these three
  return shat != oth.shat || that != oth.that || uhat != oth.uhat;
}

namespace cepgen::amplitudes {
  std::ostream& operator<<(std::ostream& os, const Kinematics& kin) {
    os << "Kin{mV2=" << kin.mv32_;
    if (kin.mv42_ != kin.mv32_)
      os << "/" << kin.mv42_;
    return os << ",shat=" << kin.shat << ",that=" << kin.that << ",uhat=" << kin.uhat << ",beta=" << kin.beta
              << ",gamma=" << kin.gamma << ",cos(theta)=" << kin.cos_theta << "->1/A=" << kin.invA << "}";
  }
}  // namespace cepgen::amplitudes
