/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2026 Laurent Forthomme
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program. If not, see <http://www.gnu.org/licenses/>.
 */

#include <algorithm>
#include <cmath>
#include <memory>
#include <numeric>
#include <vector>

#include "CepGen/Core/Exception.h"
#include "CepGen/Event/Event.h"
#include "CepGen/Modules/PartonsPhaseSpaceGeneratorFactory.h"
#include "CepGen/Modules/PhaseSpaceGeneratorFactory.h"
#include "CepGen/Physics/Cuts.h"
#include "CepGen/Physics/PDG.h"
#include "CepGen/Process/FactorisedProcess.h"
#include "CepGen/Process/PartonsPhaseSpaceGenerator.h"
#include "CepGen/Process/PhaseSpaceGenerator.h"
#include "CepGen/Utils/Math.h"
#include "CepGen/Utils/Message.h"

using namespace cepgen;
using namespace std::string_literals;

/// A 2-to-(2+N) phase space generator for a central system of N>=2 particles, based on the RAMBO algorithm
/// (Kleiss, Stirling, Ellis, Comput. Phys. Commun. 40 (1986) 359), with its massive extension.
/// The central system four-momentum P is built from three mapped quantities (transverse momentum fixed by the partons,
/// rapidity y_cs and invariant mass M_cs), then decayed isotropically (flat in phase space) in its rest frame.
class PhaseSpaceGeneratorRAMBO : public PhaseSpaceGenerator {
public:
  explicit PhaseSpaceGeneratorRAMBO(const ParametersList& params)
      : PhaseSpaceGenerator(params),
        part_psgen_(PartonsPhaseSpaceGeneratorFactory::get().build(steer<std::string>("partonsGenerator"), params_)) {}

  static ParametersDescription description() {
    auto desc = PhaseSpaceGenerator::description();
    desc.setDescription("RAMBO N-body phase space mapper");
    desc.add("partonsGenerator", ""s).setDescription("type of partons generator algorithm to use");
    return desc;
  }

  bool ktFactorised() const override {
    CG_ASSERT(part_psgen_);
    return part_psgen_->ktFactorised();
  }

  void setCentralCuts(const cuts::Central& single) override { single_limits_ = single; }

  void initialise(proc::FactorisedProcess* process) override {
    proc_ = process;
    CG_ASSERT(part_psgen_);
    if (particles_.size() < 2)
      throw CG_FATAL("PhaseSpaceGeneratorRAMBO:initialise")
          << "This phase space mapper requires at least 2 central particles (use 2to3 for a single one).";
    part_psgen_->initialise(process);
    masses_.clear();
    for (const auto& pdgid : particles_)
      masses_.emplace_back(PDG::get().mass(pdgid));
    sum_masses_ = std::accumulate(masses_.begin(), masses_.end(), 0.);

    const auto& kin_cuts = proc_->kinematics().cuts().central;
    // the system rapidity is bounded by the single particles' rapidities
    const auto lim_rap = kin_cuts.rapidity_single.truncate(Limits{-6., 6.});
    const auto lim_mass = kin_cuts.mass_sum.truncate(
        Limits{std::max(sum_masses_, NUM_LIMITS), proc_->kinematics().incomingBeams().sqrtS()});
    // the mass range must be strictly positive and above the kinematic threshold for the log mapping to be defined
    lim_mass_ = Limits{std::max({lim_mass.min(), sum_masses_, NUM_LIMITS}), lim_mass.max()};
    if (lim_mass_.max() <= lim_mass_.min())
      throw CG_FATAL("PhaseSpaceGeneratorRAMBO:initialise")
          << "Invalid central system mass range: " << lim_mass_ << ".";
    (*proc_)
        .defineVariable(m_y_cs_, proc::Process::Mapping::linear, lim_rap, "y_cs", "Central system rapidity")
        // power-law mapping: M = M_min*(M_max/M_min)^x, with the Jacobian (M ln(M_max/M_min)) handled by Process
        .defineVariable(m_m_cs_, proc::Process::Mapping::power_law, lim_mass_, "m_cs", "Central system invariant mass");
    // 4 uniform random numbers per particle (RAMBO)
    m_rho_.assign(4 * particles_.size(), 0.5);
    for (size_t i = 0; i < m_rho_.size(); ++i)
      (*proc_).defineVariable(m_rho_[i],
                              proc::Process::Mapping::linear,
                              Limits{0., 1.},
                              "rambo"s + std::to_string(i),
                              "RAMBO random number #"s + std::to_string(i));
  }

  bool generate() override {
    CG_ASSERT(part_psgen_);
    if (!part_psgen_->generatePartonKinematics())
      return false;
    central_weight_ = generateCentralKinematics();
    return utils::positive(central_weight_);
  }

  double weight() const override {
    const auto fluxes_weight = part_psgen_->fluxes();
    if (!utils::positive(fluxes_weight))
      return 0.;
    return fluxes_weight * central_weight_;
  }

  spdgids_t partons() const override {
    CG_ASSERT(part_psgen_);
    return {part_psgen_->positiveFlux().partonPdgId(), part_psgen_->negativeFlux().partonPdgId()};
  }

  spdgids_t central() const override { return particles_; }

  void setCentral(const spdgids_t& central_pdgids) override {
    if (central_pdgids.size() < 2)
      throw CG_FATAL("PhaseSpaceGeneratorRAMBO:setCentral")
          << "Invalid central particles multiplicity: expecting at least 2, got " << central_pdgids << ".";
    particles_ = central_pdgids;
  }

  // Mandelstam-like invariants are ill-defined for a generic N-body central system: processes relying on them
  // should use the 2to4 mapper instead.
  double that() const override { return 0.; }
  double uhat() const override { return 0.; }

private:
  struct Vec4 {
    double px{0.}, py{0.}, pz{0.}, e{0.};
    double p() const { return std::hypot(px, px, py); }
  };

  double generateCentralKinematics() const {
    const auto n = particles_.size();
    const auto& kin = proc_->kinematics();
    if (m_m_cs_ <= sum_masses_)
      return 0.;

    //--- four-momentum of the whole central system, from the partons' kT and the mapped (y, M)
    const auto qt_sum = (proc_->q1() + proc_->q2()).transverse();
    const auto p_cs = Momentum::fromPtYPhiM(qt_sum.pt(), m_y_cs_, qt_sum.phi(), m_m_cs_);
    const Vec4 total{p_cs.px(), p_cs.py(), p_cs.pz(), p_cs.energy()};

    //--- RAMBO decay in the system's rest frame, then boost to the laboratory frame
    std::vector<Vec4> p(n);
    const auto weight_rambo = rambo(m_m_cs_, p);
    if (!utils::positive(weight_rambo))
      return 0.;
    const auto gamma = total.e / m_m_cs_;
    const auto beta2 = (total.px * total.px + total.py * total.py + total.pz * total.pz) / (total.e * total.e);
    const auto gamma2 = beta2 > 0. ? (gamma - 1.) / beta2 : 0.;
    const Vec4 beta{total.px / total.e, total.py / total.e, total.pz / total.e, 0.};
    for (auto& k : p) {
      const auto bp = beta.px * k.px + beta.py * k.py + beta.pz * k.pz;
      const auto e = gamma * (k.e + bp);
      k.px += (gamma2 * bp + gamma * k.e) * beta.px;
      k.py += (gamma2 * bp + gamma * k.e) * beta.py;
      k.pz += (gamma2 * bp + gamma * k.e) * beta.pz;
      k.e = e;
    }

    //--- user cuts on each central particle
    std::vector<Momentum> moms;
    for (const auto& k : p) {
      moms.emplace_back(k.px, k.py, k.pz, k.e);
      const auto& mom = moms.back();
      if (!kin.cuts().central.pt_single.contains(mom.pt()) || !single_limits_.pt_single.contains(mom.pt()))
        return 0.;
      if (!kin.cuts().central.rapidity_single.contains(mom.rapidity()))
        return 0.;
      if (!kin.cuts().central.phi_single.contains(mom.phi()) || !single_limits_.phi_single.contains(mom.phi()))
        return 0.;
    }
    for (size_t i = 0; i < n; ++i) {
      proc_->pc(i) = moms.at(i);
      proc_->event()[Particle::Role::CentralSystem][i].get().setStatus(Particle::Status::FinalState);
    }

    // d^4P = d^2P_T dy M dM, with d^2P_T absorbed by the partons' kT variables (same convention as 2to4/2to3 mappers)
    // dM Jacobian (M ln(M_max/M_min)) is applied by Process::generateVariables for the power-law mapping
    return m_m_cs_ * weight_rambo;
  }

  /// Generate a flat, massive N-body phase space point in the rest frame of a system of invariant mass M.
  /// \return RAMBO weight, i.e. the (2 pi)^(4-3N)-normalised total phase space volume at this point
  double rambo(double mass_tot, std::vector<Vec4>& p) const {
    const auto n = particles_.size();
    //--- massless momenta, isotropic, with energy distributed as q0 * exp(-q0)
    std::vector<Vec4> q(n);
    Vec4 total;
    for (size_t i = 0; i < n; ++i) {
      const auto* r = &m_rho_[4 * i];
      const auto c = 2. * r[0] - 1., s = std::sqrt(std::max(0., 1. - c * c)), phi = 2. * M_PI * r[1];
      const auto e = -std::log(std::max(r[2] * r[3], 1.e-300));
      q[i] = Vec4{e * s * std::cos(phi), e * s * std::sin(phi), e * c, e};
      total.px += q[i].px, total.py += q[i].py, total.pz += q[i].pz, total.e += q[i].e;
    }
    const auto mass_q = std::sqrt(total.e * total.e - total.px * total.px - total.py * total.py - total.pz * total.pz);
    if (!utils::positive(mass_q))
      return 0.;
    //--- conformal transformation to a system of total invariant mass M, at rest
    const auto bx = -total.px / mass_q, by = -total.py / mass_q, bz = -total.pz / mass_q;
    const auto gamma = total.e / mass_q, a = 1. / (1. + gamma), x = mass_tot / mass_q;
    for (size_t i = 0; i < n; ++i) {
      const auto bq = bx * q[i].px + by * q[i].py + bz * q[i].pz;
      p[i].e = x * (gamma * q[i].e + bq);
      p[i].px = x * (q[i].px + bx * q[i].e + a * bq * bx);
      p[i].py = x * (q[i].py + by * q[i].e + a * bq * by);
      p[i].pz = x * (q[i].pz + bz * q[i].e + a * bq * bz);
    }
    //--- massive extension: rescale the 3-momenta by xi, solving sum_i sqrt(m_i^2 + xi^2 |p_i|^2) = M by Newton's method
    //    (f is increasing and convex, f(1) >= 0: monotonic convergence from xi=1)
    std::vector<double> p0(n);
    for (size_t i = 0; i < n; ++i)
      p0[i] = p[i].p();
    auto xi = 1.;
    for (size_t iter = 0; iter < MAX_NEWTON_ITER; ++iter) {
      double f = -mass_tot, df = 0.;
      for (size_t i = 0; i < n; ++i) {
        const auto e = std::sqrt(masses_[i] * masses_[i] + xi * xi * p0[i] * p0[i]);
        f += e;
        df += xi * p0[i] * p0[i] / e;
      }
      if (std::fabs(f) < NEWTON_PRECISION * mass_tot)
        break;
      if (df <= 0.)
        return 0.;
      xi -= f / df;
    }
    double prod = 1., sum_k2_over_e = 0.;
    for (size_t i = 0; i < n; ++i) {
      p[i].px *= xi, p[i].py *= xi, p[i].pz *= xi;
      p[i].e = std::sqrt(masses_[i] * masses_[i] + p[i].p() * p[i].p());
      prod *= p[i].p() / p[i].e;
      sum_k2_over_e += p[i].p() * p[i].p() / p[i].e;
    }
    if (!utils::positive(sum_k2_over_e))
      return 0.;
    //--- weights: massless total volume W0 times the massive correction (reduces to W0 for xi=1)
    const auto w0 = std::pow(2. * M_PI, 4. - 3. * n) * std::pow(0.5 * M_PI, n - 1.) * std::pow(mass_tot, 2. * n - 4.) /
                    (std::tgamma(n) * std::tgamma(n - 1.));
    return w0 * std::pow(xi, 2. * n - 3.) * prod * mass_tot / sum_k2_over_e;
  }

  static constexpr double NUM_LIMITS = 1.e-3;  ///< Numerical limits for sanity comparisons (MeV/mm-level)
  static constexpr size_t MAX_NEWTON_ITER = 50;
  static constexpr double NEWTON_PRECISION = 1.e-12;

  const std::unique_ptr<PartonsPhaseSpaceGenerator> part_psgen_;
  spdgids_t particles_;                     ///< Type of particles produced in the final state (integer values)
  std::vector<double> masses_;              ///< Central particles' masses (GeV/c^2)
  double sum_masses_{0.};                   ///< Sum of the central particles' masses (GeV/c^2)
  proc::FactorisedProcess* proc_{nullptr};  // NOT owning
  cuts::Central single_limits_;             ///< Limits to be applied on single central system's particles

  // mapped variables
  double m_y_cs_{0.};          ///< Rapidity of the central system
  Limits lim_mass_;            ///< Central system invariant mass range (GeV/c^2)
  double m_m_cs_{0.};          ///< Invariant mass of the central system
  std::vector<double> m_rho_;  ///< RAMBO uniform random numbers (4 per central particle)
  double central_weight_{0.};
};

REGISTER_PHASE_SPACE_GENERATOR("rambo", PhaseSpaceGeneratorRAMBO);
