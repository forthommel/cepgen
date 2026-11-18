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

#include <CepGen/Modules/SurvivalFactorFactory.h>
#include <CepGen/SurvivalFactors/SurvivalFactor.h>

namespace cepgen::superchic {
  /// An interfacing for SuperChic's gap survival factor \f$S^2\f$ modelling
  /// \date Nov 2026
  class SurvivalFactor : public cepgen::SurvivalFactor {
  public:
    explicit SurvivalFactor(const ParametersList& params) : cepgen::SurvivalFactor(params) {}

    static ParametersDescription description() {
      auto desc = cepgen::SurvivalFactor::description();
      desc.setDescription("SuperChic survival factor modelling interface");
      return desc;
    }
    bool singleThreaded() const override { return true; }  ///< SuperChic common blocks are process-wide
    double operator()(const Event&) const {
      if (!initialised_)
        initialise();
      return 1.;
    }

  private:
    void initialise() const { initialised_ = true; }

    mutable bool initialised_{false};
  };
}  // namespace cepgen::superchic
using SuperChicSurvivalFactor = cepgen::superchic::SurvivalFactor;
REGISTER_SURVIVAL_FACTOR("superchic", SuperChicSurvivalFactor);
