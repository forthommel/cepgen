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

#include <CepGen/Core/Exception.h>
#include <CepGen/Event/Event.h>
#include <CepGen/Modules/CouplingFactory.h>
#include <CepGen/Modules/SurvivalFactorFactory.h>
#include <CepGen/Physics/Coupling.h>
#include <CepGen/Physics/PDG.h>
#include <CepGen/SurvivalFactors/SurvivalFactor.h>

#include <array>
#include <complex>
#include <cstring>
#include <filesystem>
#include <fstream>

using namespace std::complex_literals;
using namespace std::string_literals;

extern "C" {
extern struct {
  char beam[10];
} btype_;
extern struct {
  char diff[3], diffsd[3];
} diff_com_;
extern struct {
  int diss1, diss2, offshell, elcoll;
} diss_com_;
extern struct {
  double gf, v, mt, mb, mw, me, mtau, alpha;
} ewpars_;
extern struct {
  char intag[100];
} in_;
extern struct {
  int pol;
} polarization_;
extern struct {
  int proc, dps, hel;
} process_;
extern struct {
  int s2int;
} nsurv_;
extern struct {
  double pi;
} pi_com_;
extern struct {
  double mp;
} pmass_;
extern struct {
  std::array<std::complex<double>, 10> ppa, mma, pma, mpa;
} ppamps_;
extern struct {
  int sfac;
} survival_;
extern struct {
  double s, rts, yx, mx, rtsnn;
} vars_;
extern struct {
  double x1, x2;
} x_com_;
extern void apfelinit_();
extern void calcop_();
extern void calcs2diss_();
extern void calcscreen_();
extern void dd_();
extern void gaussinit_();
extern void inpdf_();
extern void initparsr_(int&);
extern void readcoh_();
extern void readscreen_();
extern void sdcoh_();
extern void sdincoh_();
extern void schimcgam_(double&, double&, double&, double&, std::complex<double>*);
}

namespace cepgen::superchic {
  namespace {
    /// Fill a Fortran CHARACTER*N common block member (blank-padded, no NUL terminator, as expected for comparisons)
    template <size_t N>
    void setFortranString(char (&dest)[N], const std::string& src) {
      std::memset(dest, ' ', N);
      std::memcpy(dest, src.data(), std::min(src.size(), N));
    }
  }  // namespace
  /// An interfacing for SuperChic's gap survival factor \f$S^2\f$ modelling
  /// \date Nov 2026
  class SurvivalFactor : public cepgen::SurvivalFactor {
  public:
    explicit SurvivalFactor(const ParametersList& params)
        : cepgen::SurvivalFactor(params),
          alpha_em_(AlphaEMFactory::get().build(steer<ParametersList>("alphaEM"))),
          survival_factor_model_(steer<int>("survivalFactorModel")),
          integration_parameter_(steer<int>("integrationParameter")),
          tag_(steer<std::string>("tag")),
          regenerate_(steer<bool>("regenerate")),
          proc_(steer<int>("proc")) {
      if (const auto model = steer<std::string>("model"); model == "scalar")
        amplitudes_ = std::array<std::complex<double>, 4>{1., 1., 0., 0.};
      else if (model == "pseudoscalar")
        amplitudes_ = std::array<std::complex<double>, 4>{1., -1., 0., 0.};
      else if (model == "tensor")
        amplitudes_ = std::array<std::complex<double>, 4>{0., 0., 1., 1.};
    }

    static ParametersDescription description() {
      auto desc = cepgen::SurvivalFactor::description();
      desc.setDescription("SuperChic survival factor modelling interface");
      desc.add("alphaEM", AlphaEMFactory::get().describeParameters("fixed"))
          .setDescription("alpha(EM) modelling for scale evolution");
      desc.add("survivalFactorModel", 4)
          .allow(1, ""s)  //TODO
          .allow(2, ""s)
          .allow(3, ""s)
          .allow(4, ""s)
          .setDescription("SuperChic soft survival model");
      desc.add("integrationParameter", 16).setDescription("SuperChic survival integration parameter");
      desc.add("tag", "user"s).setDescription("SuperChic 'intag' used when running its 'init' program");
      desc.add("regenerate", false)
          .setDescription("force the (re)computation of SuperChic's screening amplitude, even if its file exists");
      desc.add("proc", 0).setDescription(
          "SuperChic process number, used for (single/double) dissociative events only (requires its off-shell "
          "amplitude)");
      return desc;
    }
    bool singleThreaded() const override { return true; }  ///< SuperChic common blocks are process-wide
    double operator()(const Event& event) const override {
      if (!initialised_)
        initialise(event);

      const auto &ib1 = event.oneWithRole(Particle::Role::IncomingBeam1).momentum(),
                 &ib2 = event.oneWithRole(Particle::Role::IncomingBeam2).momentum();
      const auto &ob1 = event.oneWithRole(Particle::Role::OutgoingBeam1).momentum(),
                 &ob2 = event.oneWithRole(Particle::Role::OutgoingBeam2).momentum();
      double ob1_px = ob1.px(), ob1_py = ob1.py(), ob2_px = ob2.px(), ob2_py = ob2.py();
      process_.proc = 0;
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
      pmass_.mp = PDG::get().mass(PDG::proton);
      pi_com_.pi = M_PI;
      ewpars_.alpha = (*alpha_em_)(0. /*FIXME*/);

      nsurv_.s2int = integration_parameter_;

      diss_com_.diss1 = 0;
      diss_com_.diss2 = 0;
      setFortranString(btype_.beam, "prot");
      setFortranString(diff_com_.diff, "el");
      polarization_.pol = 1;  // only the first polarisation state amplitudes are filled below

      // fill the common blocks first: both the screening amplitude generation and its reading rely on them
      setFortranString(in_.intag, tag_);
      vars_.rts = first_event.cmEnergy();
      vars_.s = vars_.rts * vars_.rts;

      gaussinit_();
      int surv = survival_factor_model_;
      initparsr_(surv);

      // equivalent of SuperChic's 'init' program (restricted to what the survival factor requires):
      // the screening amplitude is written to 'inputs/screening<tag>.dat' in the current working directory
      const auto screening_file = std::filesystem::path("inputs") / ("screening"s + tag_ + ".dat"s);
      if (regenerate_ || !std::filesystem::exists(screening_file)) {
        CG_INFO("SuperChic:initialise") << "Computing the screening amplitude for tag '" << tag_
                                        << "' (sqrt(s) = " << vars_.rts << " GeV, model " << surv << ") into "
                                        << screening_file << ".";
        std::filesystem::create_directories("inputs");
        calcop_();  // proton opacity (needed for the screening amplitude computation)
        calcscreen_();
      } else
        calcop_();  // opacity table: assumed needed by the survival amplitudes as well
      readscreen_();

      ppamps_.ppa[0] = amplitudes_.at(0);
      ppamps_.mma[0] = amplitudes_.at(1);
      ppamps_.pma[0] = amplitudes_.at(2);
      ppamps_.mpa[0] = amplitudes_.at(3);

      initialised_ = true;
    }

    const std::unique_ptr<Coupling> alpha_em_;
    const int survival_factor_model_;
    const int integration_parameter_;
    const std::string tag_;
    const bool regenerate_;
    const int proc_;

    std::array<std::complex<double>, 4> amplitudes_;

    mutable bool initialised_{false};
  };
}  // namespace cepgen::superchic
using SuperChicSurvivalFactor = cepgen::superchic::SurvivalFactor;
REGISTER_SURVIVAL_FACTOR("superchic", SuperChicSurvivalFactor);
