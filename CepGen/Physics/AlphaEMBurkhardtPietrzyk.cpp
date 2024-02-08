/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2024  Laurent Forthomme
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

#include "CepGen/Core/SteeredObject.h"
#include "CepGen/Modules/CouplingFactory.h"
#include "CepGen/Physics/Constants.h"
#include "CepGen/Physics/Coupling.h"
#include "CepGen/Physics/PDG.h"

namespace cepgen {
  /// Electromagnetic alpha running calculator
  /// \note Implementation stolen from Grape MC
  /// \cite Burkhardt:1995tt
  class AlphaEMBurkhardtPietrzyk final : public Coupling {
  public:
    explicit AlphaEMBurkhardtPietrzyk(const ParametersList& params) : Coupling(params) {
      for (const auto& region : steer<std::vector<ParametersList> >("regions"))
        regions_.emplace_back(region);
      for (auto pdgid : {11, 13, 15})
        lepton_sq_masses_.emplace_back(std::pow(PDG::get().mass(pdgid), 2));
    }

    struct ParameterisationRegion : SteeredObject<ParameterisationRegion> {
      explicit ParameterisationRegion(const ParametersList& params) : SteeredObject(params) {
        (*this).add("q2max", q2max).add("a", a).add("b", b).add("c", c);
      }
      static ParametersDescription description() {
        auto desc = ParametersDescription();
        desc.setDescription("parameterisation of a Q^2 kinematic region");
        desc.add<double>("q2max", -1.).setDescription("maximum Q^2 range of validity, in GeV^2");
        desc.add<double>("a", -1.).setDescription("constant A parameter");
        desc.add<double>("b", -1.).setDescription("linear B parameter");
        desc.add<double>("c", -1.).setDescription("logarithmic C parameter");
        return desc;
      }
      double q2max, a, b, c;
    };

    static ParametersDescription description() {
      auto desc = Coupling::description();
      desc.setDescription("Burkhardt & Pietrzyk alpha(EM) evolution algorithm");
      desc.addParametersDescriptionVector("regions",
                                          ParameterisationRegion::description(),
                                          {ParametersList()
                                               .set<double>("q2max", 2. * 2.)
                                               .set<double>("a", 0.)
                                               .set<double>("b", 2.28770e-3)
                                               .set<double>("c", 4.08041425),
                                           ParametersList()
                                               .set<double>("q2max", 4. * 4.)
                                               .set<double>("a", 0.)
                                               .set<double>("b", 2.51507e-3)
                                               .set<double>("c", 3.09624477),
                                           ParametersList()
                                               .set<double>("q2max", 10. * 10.)
                                               .set<double>("a", 0.)
                                               .set<double>("b", 2.79328e-3)
                                               .set<double>("c", 2.07463133),
                                           ParametersList()
                                               .set<double>("q2max", 91.2 * 91.2)
                                               .set<double>("a", 1.22270e-3)
                                               .set<double>("b", 2.96694e-3)
                                               .set<double>("c", 1.),
                                           ParametersList()
                                               .set<double>("q2max", 1.e5 * 1.e5)
                                               .set<double>("a", 1.64178e-3)
                                               .set<double>("b", 2.92051e-3)
                                               .set<double>("c", 1.)})
          .setDescription("list of Q^2 regions parameterised for the hadronic contribution to the vacuum polarisation");
      return desc;
    }

    double operator()(double q) const override {
      const auto q2 = q * q;
      double cr = 0.;
      for (const auto& ml2 : lepton_sq_masses_) {
        if (q2 > ml2)
          cr += std::log(q2 / ml2);
      }
      return 1. / (1. - (prefac_ * cr + re_ph(q2)));
    }

  private:
    static constexpr double prefac_ = constants::ALPHA_EM * M_1_PI / 3.;

    /// Hadronic contribution to the vacuum polarisation
    double re_ph(double q2) const {
      for (const auto& reg : regions_) {
        if (q2 < reg.q2max)
          return reg.a + reg.b * std::log1p(reg.c * q2);
      }
      return 0.;
    }
    std::vector<ParameterisationRegion> regions_;
    std::vector<double> lepton_sq_masses_;
  };
}  // namespace cepgen

REGISTER_ALPHAEM_MODULE("burkhardtpietrzyk", AlphaEMBurkhardtPietrzyk);
