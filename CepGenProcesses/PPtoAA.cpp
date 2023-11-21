/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2018-2023  Laurent Forthomme
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

#include "CepGen/Core/Exception.h"
#include "CepGen/Event/Event.h"
#include "CepGen/Modules/ProcessFactory.h"
#include "CepGen/Modules/RandomGeneratorFactory.h"
#include "CepGen/Physics/Constants.h"
#include "CepGen/Physics/PDG.h"
#include "CepGen/Process/Process2to4.h"
#include "CepGen/Utils/RandomGenerator.h"
#include "CepGen/Utils/String.h"

using namespace cepgen;

/// \brief Compute the matrix element for a CE \f$\gamma\gamma\rightarrow\gamma\gamma\f$ process using \f$k_{\rm T}\f$-factorization approach
class PPtoAA : public cepgen::proc::Process2to4 {
public:
  explicit PPtoAA(const ParametersList& params)
      : Process2to4(params, PDG::photon),
        method_((Method)params.get<int>("method", 1)),
        pol_state_((Polarisation)params.get<int>("polarisationStates", 0)),
        par_alp_(params.get<ParametersList>("alpParameters")),
        rnd_gen_(RandomGeneratorFactory::get().build(steer<ParametersList>("randomGenerator"))) {
    switch (method_) {
      case Method::sm:
      default: {
        //FIXME
        switch (pol_state_) {
          case Polarisation::LL:
            pol_gam1_ = pol_gam2_ = {0};
            break;
          case Polarisation::LT:
            pol_gam1_ = {0};
            pol_gam2_ = {-1, 1};
            break;
          case Polarisation::TL:
            pol_gam1_ = {-1, 1};
            pol_gam2_ = {0};
            break;
          case Polarisation::TT:
            pol_gam1_ = pol_gam2_ = {-1, 1};
            break;
          default:
          case Polarisation::full:
            pol_gam1_ = pol_gam2_ = {-1, 1};
            break;
        }
      } break;
      case Method::alp:
        pol_gam1_ = pol_gam2_ = {0};
        break;
    }
    CG_DEBUG("PPtoAA:mode") << "matrix element computation method: " << (int)method_ << ".";
  }

  proc::ProcessPtr clone() const override { return proc::ProcessPtr(new PPtoAA(parameters())); }

  static ParametersDescription description() {
    auto desc = Process2to4::description();
    desc.setDescription("ɣɣ → ɣɣ light-by-light process");
    desc.add<ParametersDescription>("alpParameters", ALPParameters::description());
    desc.add<ParametersDescription>("randomGenerator", ParametersDescription().setName<std::string>("stl"))
        .setDescription("random number generator engine");
    return desc;
  }

  enum class Polarisation { full = 0, LL = 1, LT = 2, TL = 3, TT = 4 };
  enum class Method { sm = 0, alp = 1 };

private:
  void prepareProcessKinematics() override {}
  double computeCentralMatrixElement() const override {
    switch (method_) {
      case Method::sm:
        throw CG_FATAL("PPtoAA") << "Process not yet handled!";
        return smOnShellME();
      case Method::alp:
        return alpOnShellME();
      default:
        throw CG_FATAL("PPtoAA") << "Process not yet handled!";
    }
  }

  double smOnShellME() const { return 0.; }
  double alpOnShellME() const;

  const Method method_;
  const Polarisation pol_state_;
  /// List of parameters for axion-like particles exchanges
  struct ALPParameters : SteeredObject<ALPParameters> {
    explicit ALPParameters(const ParametersList& params)
        : SteeredObject(params),
          mass(steer<double>("mass")),
          mass2(mass * mass),
          coupling(steer<double>("coupling")),
          width(coupling * coupling * mass * mass * mass / 64. * M_1_PI),
          width2(width * width),
          pseudo_scalar(steer<double>("pseudoScalar")) {}
    static ParametersDescription description() {
      auto desc = ParametersDescription();
      desc.setDescription("collection of parameters for the axion-like particle description");
      desc.add<double>("mass", 0.).setDescription("intermediate particle mass, in GeV/c^2");
      desc.add<double>("coupling", 0.).setDescription("intermediate particle coupling value");
      desc.add<bool>("pseudoScalar", false).setDescription("intermediate particle spin");
      return desc;
    }
    const double mass, mass2, coupling;
    const double width, width2;
    const bool pseudo_scalar;
  } par_alp_;
  const std::unique_ptr<utils::RandomGenerator> rnd_gen_;

  std::vector<short> pol_gam1_, pol_gam2_;
};

double PPtoAA::alpOnShellME() const {
  const auto s_hat = shat();

  const auto norma = (0.5 * par_alp_.coupling) * (s_hat * 0.5);
  const auto ampli_pp = norma, ampli_mm = (par_alp_.pseudo_scalar ? +1. : -1.) * norma;

  const Limits mass_range{par_alp_.mass - 10. * par_alp_.width, par_alp_.mass + 10. * par_alp_.width},
      al_range = mass_range.compute(
          [this](double ext) { return std::atan((ext * ext - par_alp_.mass2) / par_alp_.width / par_alp_.mass); });
  const auto norm = al_range.range();
  const auto al1 = rnd_gen_->uniform(al_range.min(), al_range.max());
  const auto mout = std::sqrt(std::tan(al1) * par_alp_.mass * par_alp_.width + par_alp_.mass2);

  //--- compute Jacobian
  auto jalp = norm * par_alp_.width * par_alp_.mass * (1. + std::pow(std::tan(al1), 2));
  //--- up to this point just change of variables
  //--- now BW
  jalp *= mout * par_alp_.width / norm;
  jalp /= (std::pow(par_alp_.mass2 - mout * mout, 2) + mout * mout * par_alp_.width2);

  return (ampli_pp * ampli_pp + ampli_mm * ampli_mm) * jalp / s_hat / s_hat;
}
// register process
REGISTER_PROCESS("pptoaa", PPtoAA);
