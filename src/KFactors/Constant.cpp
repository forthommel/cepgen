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

#include "CepGen/Event/Event.h"
#include "CepGen/KFactors/KFactor.h"
#include "CepGen/Modules/KFactorFactory.h"

using namespace std::string_literals;

namespace cepgen::kfactor {
  /// Constant K-factor, independent of kinematics (and of the scale)
  class Constant final : public KFactor {
  public:
    explicit Constant(const ParametersList& params) : KFactor(params), value_(steer<double>("value")) {}

    static ParametersDescription description() {
      auto desc = KFactor::description();
      desc.setDescription("Constant K-factor");
      desc.add("value", 1.).setDescription("K-factor value");
      return desc;
    }

  protected:
    double compute(const Event&, double) const override { return value_; }

  private:
    const double value_;
  };
}  // namespace cepgen::kfactor
using ConstantKFactor = cepgen::kfactor::Constant;
REGISTER_KFACTOR("constant", ConstantKFactor);
