/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2013-2021  Laurent Forthomme
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
#include "CepGen/Physics/Constants.h"
#include "CepGen/Physics/PDG.h"
#include "CepGen/Process/Process2to4.h"
#include "CepGen/Utils/String.h"

namespace cepgen {
  namespace proc {
    /// \brief Compute the matrix element for a CE \f$\gamma\gamma\rightarrow\gamma\gamma\f$ process using \f$k_{\rm T}\f$-factorization approach
    class PPtoAA : public Process2to4 {
    public:
      PPtoAA(const ParametersList& params = ParametersList());

      ProcessPtr clone() const override { return ProcessPtr(new PPtoAA(*this)); }

      static ParametersDescription description() {
        auto desc = Process2to4::description();
        desc.setDescription("ɣɣ → ɣɣ light-by-light process");
        return desc;
      }

      enum class Polarisation { full = 0, LL = 1, LT = 2, TL = 3, TT = 4 };
      enum class Method { sm = 0, alp = 1 };

    private:
      void prepareProcessKinematics() override {}
      double computeCentralMatrixElement() const override;

      double smOnShellME() const;
      double alpOnShellME() const;

      Method method_;
      Polarisation pol_state_;
      std::vector<short> pol_gam1_, pol_gam2_;

      /// List of parameters for axion-like particles exchanges
      ParametersList par_alp_;
    };

    PPtoAA::PPtoAA(const ParametersList& params)
        : Process2to4(params, {PDG::photon, PDG::photon}, PDG::photon),
          method_((Method)params.get<int>("method", 1)),
          pol_state_((Polarisation)params.get<int>("polarisationStates", 0)),
          par_alp_(params.get<ParametersList>("alpParameters")) {
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

    double PPtoAA::computeCentralMatrixElement() const {
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

    double PPtoAA::smOnShellME() const { return 0.; }

    double PPtoAA::alpOnShellME() const {
      const double s_hat = shat();

      const double mA = par_alp_.get<double>("mass"), gA = par_alp_.get<double>("coupling");
      const bool pseudoscalar = par_alp_.get<bool>("pseudoScalar", true);

      const double norma = (0.5 * gA) * (s_hat * 0.5);
      double ampli_pp = norma, ampli_mm = (pseudoscalar ? +1. : -1.) * norma;

      const double width = gA * gA * pow(mA, 3) / 64. * M_1_PI;

      const double mminr = mA - 10. * width, mmaxr = mA + 10. * width;
      const double almin = atan(-(mA * mA - mminr * mminr) / width / mA);
      const double almax = atan(-(mA * mA - mmaxr * mmaxr) / width / mA);

      const double norm = almax - almin;

      const double al1 = almin + (almax - almin) * drand();
      const double mout = sqrt(tan(al1) * mA * width + mA * mA);

      //--- compute Jacobian
      double jalp = (almax - almin) * width * mA * (1. + pow(tan(al1), 2));
      //--- up to this point just change of variables
      //--- now BW
      jalp *= mout * width / norm;
      jalp /= (pow(mA * mA - mout * mout, 2) + mout * mout * width * width);

      return (ampli_pp * ampli_pp + ampli_mm * ampli_mm) * jalp;
    }
  }  // namespace proc
}  // namespace cepgen
// register process
typedef cepgen::proc::PPtoAA PPtoAA;
REGISTER_PROCESS("pptoaa", PPtoAA);
