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

#include <CepGen/Event/Event.h>
#include <CepGen/Modules/SurvivalFactorFactory.h>
#include <CepGen/SurvivalFactors/SurvivalFactor.h>

#include <complex>

extern "C" {
extern struct {
  int sfac;
} survival_;
extern void schimcgam_(double&, double&, double&, double&, std::complex<double>*);
}

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
    double operator()(const Event& event) const {
      if (!initialised_)
        initialise(event);

      const auto &ob1 = event.oneWithRole(Particle::Role::OutgoingBeam1).momentum(),
                 &ob2 = event.oneWithRole(Particle::Role::OutgoingBeam2).momentum();
      double ob1_px = ob1.px(), ob1_py = ob1.py(), ob2_px = ob1.px(), ob2_py = ob2.py();
      std::array<std::complex<double>, 10> wts, wtb;
      survival_.sfac = 1;
      schimcgam_(ob1_px, ob1_py, ob2_px, ob2_py, wts.data());
      survival_.sfac = 0;
      schimcgam_(ob1_px, ob1_py, ob2_px, ob2_py, wtb.data());
      const auto num_s = std::pow(std::abs(wts.at(0)), 2), num_b = std::pow(std::abs(wtb.at(0)), 2);
      if (num_b > 0.)
        return num_s / num_b;
      return 0.;
    }

  private:
    void initialise(const Event& first_event) const {
      const auto &ib1 = first_event.oneWithRole(Particle::Role::IncomingBeam1).momentum(),
                 &ib2 = first_event.oneWithRole(Particle::Role::IncomingBeam2).momentum();
      initialised_ = true;
    }

    mutable bool initialised_{false};
  };
}  // namespace cepgen::superchic
using SuperChicSurvivalFactor = cepgen::superchic::SurvivalFactor;
REGISTER_SURVIVAL_FACTOR("superchic", SuperChicSurvivalFactor);
