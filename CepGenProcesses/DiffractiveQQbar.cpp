/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2023-2024  Barbara Linek, Wolfgang Schaefer, Marta Luszczak
 *                2024       Laurent Forthomme
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

#include <fstream>

#include "CepGen/Core/Exception.h"
#include "CepGen/Event/Event.h"
#include "CepGen/Modules/ProcessFactory.h"
#include "CepGen/Physics/Constants.h"
#include "CepGen/Physics/PDG.h"
#include "CepGen/Physics/Utils.h"
#include "CepGen/Process/Process.h"
#include "CepGen/Utils/Filesystem.h"
#include "CepGen/Utils/GridDrawer.h"
#include "CepGen/Utils/GridHandler.h"
#include "CepGen/Utils/Math.h"
#include "CepGen/Utils/String.h"

using namespace cepgen;

class DiffractiveQQbar final : public cepgen::proc::Process {
public:
  explicit DiffractiveQQbar(const ParametersList& params)
      : Process(params),
        pair_(steer<ParticleProperties>("pair")),
        fy_hatta_grid_path_(steerPath("gridPathHatta")),
        mq2_(pair_.mass * pair_.mass) {}

  proc::ProcessPtr clone() const override { return proc::ProcessPtr(new DiffractiveQQbar(parameters())); }

  static ParametersDescription description() {
    auto desc = Process::description();
    desc.setDescription("diffractive q-qbar");
    desc.add<int>("pair", PDG::up);
    desc.add<std::string>("gridPathHatta", "F0_KTv3_logk.dat").setDescription("path to the Hatta et al. F_Y grid");
    return desc;
  }

  void addEventContent() override {
    proc::Process::setEventContent(
        {{Particle::Role::IncomingBeam1, {PDG::electron}},
         {Particle::Role::IncomingBeam2, {PDG::proton}},
         {Particle::Role::Parton1, {PDG::photon}},
         {Particle::Role::Parton2, {PDG::photon}},
         {Particle::Role::OutgoingBeam1, {PDG::electron}},
         {Particle::Role::OutgoingBeam2, {PDG::proton}},
         {Particle::Role::CentralSystem, {(spdgid_t)(+pair_.pdgid), (spdgid_t)(-pair_.pdgid)}}});
  }

private:
  void prepareKinematics() override {
    defineVariable(m_qt1_, Mapping::linear, {0.05, 50.}, "qt1", "first incoming parton transverse virtuality");
    defineVariable(m_qt2_, Mapping::linear, {0.05, 50.}, "qt2", "second incoming parton transverse virtuality");
    defineVariable(m_phi1_, Mapping::linear, {0., 2. * M_PI}, "phi1", "first incoming parton azimuthal angle");
    defineVariable(m_phi2_, Mapping::linear, {0., 2. * M_PI}, "phi2", "second incoming parton azimuthal angle");
    defineVariable(m_pt_, Mapping::linear, {0., 100.}, "pt", "jets transverse momentum");
    defineVariable(
        m_delta_, Mapping::linear, {0., 3.}, "delta", "jets transverse momentum difference in gamma-p frame");
    defineVariable(
        m_phi_, Mapping::linear, {0., 2. * M_PI}, "phi", "jets transverse momentum azimuthal angle in gamma-p frame");
    defineVariable(m_z_, Mapping::linear, {1.e-5, 1.}, "z");
    defineVariable(m_y_in_, Mapping::linear, {0.05, 0.7}, "y_in");
    defineVariable(m_q2_, Mapping::linear, {1.e-5, 1.e3}, "q2", "photon virtuality");

    {  // initialise the F_Y interpolation grid
      if (!utils::fileExists(fy_hatta_grid_path_))
        throw CG_FATAL("DiffractiveQQbar")
            << "Failed to load the Hatta et al. F_Y interpolation grid from '" << fy_hatta_grid_path_ << "'!";
      CG_INFO("DiffractiveQQbar") << "Loading Hatta et al. F_Y values from '" << fy_hatta_grid_path_ << "' file.";
      std::ifstream grid_file(fy_hatta_grid_path_);
      std::string tmp;
      while (std::getline(grid_file, tmp)) {
        const auto toks = utils::split(tmp, ' ', true);
        if (toks.size() != 4)
          continue;
        fy_hatta_.insert({std::stod(toks.at(0)), std::stod(toks.at(1)), std::stod(toks.at(2))},
                         {std::stod(toks.at(3))});
      }
      fy_hatta_.initialise();
      utils::GridDrawer::draw(fy_hatta_);
    }
  }

  double computeWeight() override {
    q1() = Momentum::fromPtEtaPhiE(m_qt1_, 0., m_phi1_);
    q2() = Momentum::fromPtEtaPhiE(m_qt2_, 0., m_phi2_);
    const auto diff_pt = Momentum::fromPtEtaPhiE(m_pt_, 0., 0.);
    const auto pair_pt = Momentum::fromPtEtaPhiE(m_delta_, 0., m_phi_);
    //const auto pt1 = diff_pt + 0.5 * pair_pt, pt2 = diff_pt - 0.5 * pair_pt;
    const auto z1 = m_z_, z2 = 1. - z1;

    //const auto mt2 = (pt1.p2() + mq2_) / z1 + (pt2.p2() + mq2_) / z2;
    //const auto m2_qqbar = mt2 - m_delta_ * m_delta_, m_qqbar = std::sqrt(m2_qqbar);

    const auto w2_gammap = m_y_in_ * s();
    const auto xbj = utils::xBj(m_q2_, mp2_, w2_gammap);
    //const auto x_pom = xbj / utils::xBj(m_q2_, 0., m2_qqbar);
    //const auto eta1 = std::log((w2_gammap * z1) / (mp_ * std::sqrt(pt1.p2() + mq2_) * std::exp(std::acosh(ep / mp_))));
    //const auto eta2 = std::log((w2_gammap * z2) / (mp_ * std::sqrt(pt2.p2() + mq2_) * std::exp(std::acosh(ep / mp_))));

    const auto y_evol = std::log(1.e-2 / xbj);
    const auto tt = fy_hatta_.eval({y_evol, m_delta_, std::log10(m_qt1_)}).at(0),
               tt_prime = fy_hatta_.eval({y_evol, m_delta_, std::log10(m_qt2_)}).at(0);

    const auto alpha_em = alphaEM((q1() + q2()).mass());

    const auto factor = 0.5 * alpha_em / num_colours_ * ef_sq_;
    const auto num_t = (pair_pt - q1()) * (pair_pt - q2());
    const auto den = (z1 * z2 * m_q2_ + (diff_pt - q1()).p2()) * (z1 * z2 * m_q2_ + (diff_pt - q2()).p2());
    const auto mat2_t = 0.25 * (z1 * z1 + z2 * z2) * tt * tt_prime * num_t,
               mat2_l = (z1 * z1 * z2 * z2 * m_q2_) * tt * tt_prime;

    const auto mat2 = factor / den * ((1. - m_y_in_ + 0.5 * m_y_in_ * m_y_in_) * mat2_t + (1. - m_y_in_) * mat2_l);

    const auto factor2 = alpha_em * M_1_PI / m_y_in_ / m_q2_;
    const auto weight = factor2 * mat2 * (m_qt1_ * m_qt2_ * m_pt_ * m_delta_);
    if (utils::positive(weight))
      return weight;
    return 0.;
  }

  void fillKinematics() override {
    //event().oneWithRole(Particle::Role::OutgoingBeam1).momentum().setMass(m_delta_);
  }

  const ParticleProperties pair_;
  const std::string fy_hatta_grid_path_;

  static constexpr int num_colours_ = 3;
  static constexpr double ef_sq_ = std::pow(2. / 3, 2) + std::pow(1. / 3, 2) + std::pow(1. / 3, 2);
  const double mq2_;
  GridHandler<3, 1> fy_hatta_{GridType::linear};

  double m_qt1_, m_qt2_, m_phi1_, m_phi2_;
  double m_pt_, m_delta_, m_phi_;
  double m_z_;
  double m_y_in_;
  double m_q2_;
};
// register process
REGISTER_PROCESS("diffqqbar", DiffractiveQQbar);
