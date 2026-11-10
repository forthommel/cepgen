/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2016-2025  Laurent Forthomme
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

#include <algorithm>
#include <cmath>

#include "CepGen/Core/Exception.h"
#include "CepGen/Core/RunParameters.h"
#include "CepGen/Event/Event.h"
#include "CepGen/Event/Particle.h"
#include "CepGen/Physics/PDG.h"
#include "CepGenPythia8/CepGenEvent.h"

using namespace Pythia8;

/// Convert a CepGen particle momentum into its Pythia8 counterpart
Vec4 momToVec4(const cepgen::Momentum& mom) { return {mom.px(), mom.py(), mom.pz(), mom.energy()}; }

CepGenEvent::CepGenEvent() : LHAup(3), mp_(cepgen::PDG::get().mass(cepgen::PDG::proton)), mp2_(mp_ * mp_) {}

void CepGenEvent::initialise(const cepgen::RunParameters& params) {
  params_ = &params;
  inel1_ = !params_->kinematics().incomingBeams().positive().elastic();
  inel2_ = !params_->kinematics().incomingBeams().negative().elastic();

  setBeamA(params_->kinematics().incomingBeams().positive().integerPdgId(),
           params_->kinematics().incomingBeams().positive().momentum().pz());
  setBeamB(params_->kinematics().incomingBeams().negative().integerPdgId(),
           params_->kinematics().incomingBeams().negative().momentum().pz());
  //addProcess( 0, params_->integration().result, params_->integration().err_result, 100. );
}

void CepGenEvent::addComments(const std::string& comments) {
#if PYTHIA_VERSION_INTEGER >= 8200
  osLHEF << comments;
#else
  CG_WARNING("CepGenEvent:addComments") << "Pythia 8 is too outdated... Unused comments: " << comments;
#endif
}

void CepGenEvent::setCrossSection(int id, double cross_section, double cross_section_err) {
  addProcess(0, cross_section, cross_section_err, 100.);
  setXSec(id, cross_section);
  setXErr(id, cross_section_err);
  //listInit();
}

namespace {
  /// Light-cone component of a four-momentum along a beam direction (sign = +1 for a +z beam, -1 for a -z beam)
  double lightCone(const Vec4& mom, double sign) { return mom.e() + sign * mom.pz(); }
  /// Massless four-momentum collinear to a beam, carrying a fraction x of its light-cone momentum
  Vec4 collinear(double x, double beam_light_cone, double sign) {
    return Vec4(0., 0., sign * 0.5 * x * beam_light_cone, 0.5 * x * beam_light_cone);
  }
}  // namespace

bool CepGenEvent::feedEvent(const cepgen::Event& ev, const Type& type, double event_weight) {
  const double scale = ev(cepgen::Particle::Role::Intermediate)[0].momentum().mass();
  setProcess(0, event_weight, scale, ev.metadata("alphaEM"), ev.metadata("alphaS"));
  boost_lha_to_cg_.reset();
  boost_cg_to_lha_.reset();
  collinearised_ = false;

  unsigned short colour_index = MIN_COLOUR_INDEX;

  if (type == Type::centralAndBeamRemnants) {  // full event content, with a collinear initial state
    if (!feedCollinearInitialState(ev, scale, colour_index))
      return false;
  } else {
    //===========================================================================================
    // incoming partons
    //===========================================================================================

    addCepGenParticle(ev(cepgen::Particle::Role::Parton1)[0], -2);
    addCepGenParticle(ev(cepgen::Particle::Role::Parton2)[0], -2);

    if (type == Type::centralAndFullBeamRemnants)  // full beam remnants content
      for (const auto& forward_system : {cepgen::Particle::Role::OutgoingBeam1, cepgen::Particle::Role::OutgoingBeam2})
        for (const auto& p : ev(forward_system))
          addCepGenParticle(p, INVALID_ID, findMothers(ev, p));
  }

  //=============================================================================================
  // central system
  //=============================================================================================

  const unsigned short central_colour = colour_index++;
  for (const auto& p : ev(cepgen::Particle::Role::CentralSystem)) {
    std::pair<int, int> colours = {0, 0}, mothers = {1, 2};
    if (type != Type::centralAndBeamRemnants)
      mothers = findMothers(ev, p);
    try {
      if (cepgen::PDG::get().colours(p.pdgId()) > 1) {
        if (p.integerPdgId() > 0)  // particle
          colours.first = central_colour;
        else  // anti-particle
          colours.second = central_colour;
      }
    } catch (const cepgen::Exception&) {
    }
    int status = 1;
    if (type == Type::centralAndFullBeamRemnants && p.status() == cepgen::Particle::Status::Resonance)
      status = 2;
    addCepGenParticle(p, status, mothers, colours);
  }
  return true;
}

bool CepGenEvent::feedCollinearInitialState(const cepgen::Event& ev, double scale, unsigned short& colour_index) {
  using cepgen::Particle;
  const Particle* beams[2] = {&ev(Particle::Role::IncomingBeam1)[0], &ev(Particle::Role::IncomingBeam2)[0]};
  const Particle* partons[2] = {&ev(Particle::Role::Parton1)[0], &ev(Particle::Role::Parton2)[0]};
  const Particle* outgoing[2] = {&ev(Particle::Role::OutgoingBeam1)[0], &ev(Particle::Role::OutgoingBeam2)[0]};
  const bool inel[2] = {inel1_, inel2_};

  Vec4 mom_beam[2], mom_parton[2];
  double sign[2], beam_lc[2];
  for (size_t i = 0; i < 2; ++i) {
    mom_beam[i] = momToVec4(beams[i]->momentum());
    mom_parton[i] = momToVec4(partons[i]->momentum());
    sign[i] = mom_beam[i].pz() > 0. ? 1. : -1.;
    beam_lc[i] = lightCone(mom_beam[i], sign[i]);
  }
  Vec4 mom_central;
  for (const auto& p : ev(Particle::Role::CentralSystem))
    mom_central += momToVec4(p.momentum());
  const double m2_central = mom_central.m2Calc();
  if (m2_central <= 0.)
    return false;

  //--- momentum fractions of the collinear initiators
  double x[2] = {0., 0.};
  for (size_t i = 0; i < 2; ++i) {
    if (!inel[i])
      continue;
    const double q_lc_opposite = lightCone(mom_parton[i], -sign[i]);  // q^- for a +z beam
    x[i] = std::fabs(q_lc_opposite) > 1.e-12 ? mom_parton[i].m2Calc() / (beam_lc[i] * q_lc_opposite)
                                             : lightCone(mom_parton[i], sign[i]) / beam_lc[i];
  }
  if (!inel[0] && !inel[1]) {  // both elastic: keep the central-system mass and rapidity
    const double c_1 = lightCone(mom_central, sign[0]), c_2 = lightCone(mom_central, sign[1]);
    if (c_1 <= 0. || c_2 <= 0.)
      return false;
    x[0] = std::sqrt(m2_central * c_1 / c_2) / beam_lc[0];
    x[1] = std::sqrt(m2_central * c_2 / c_1) / beam_lc[1];
  } else if (inel[0] != inel[1]) {  // single dissociation: (x_e P_e + q_r)^2 = M^2
    const size_t el = inel[0] ? 1 : 0, rec = 1 - el;
    const double denominator = 2. * (collinear(1., beam_lc[el], sign[el]) * mom_parton[rec]);
    if (denominator <= 0.)
      return false;
    x[el] = (m2_central - mom_parton[rec].m2Calc()) / denominator;
  }
  for (const auto& xi : x)
    if (xi <= 0. || xi >= 1.)
      return false;

  //--- collinear initiators, and transformation between the collinear hard system and the CepGen central system
  Vec4 mom_init[2], mom_hard;
  for (size_t i = 0; i < 2; ++i) {
    mom_init[i] = collinear(x[i], beam_lc[i], sign[i]);
    mom_hard += inel[i] ? mom_parton[i] : mom_init[i];
  }
  boost_lha_to_cg_.reset();
  boost_lha_to_cg_.bstback(mom_hard);  // collinear hard system -> its rest frame
  boost_lha_to_cg_.bst(mom_central);   // rest frame -> CepGen central system
  boost_cg_to_lha_ = boost_lha_to_cg_;
  boost_cg_to_lha_.invert();
  collinearised_ = !inel[0] || !inel[1];  // for double dissociation the transformation is the identity

  //--- flavour / x value of hard-process initiators
  int pdg_id[2], colour[2] = {0, 0};
  for (size_t i = 0; i < 2; ++i) {
    pdg_id[i] = inel[i] ? 2 : partons[i]->integerPdgId();  //FIXME select quark flavours accordingly
    if (inel[i])
      colour[i] = colour_index++;
  }
  setIdX(pdg_id[0], pdg_id[1], x[0], x[1]);
  setPdf(pdg_id[0], pdg_id[1], x[0], x[1], scale, 0., 0., false);

  //===========================================================================================
  // incoming (massless, collinear) initiators
  //===========================================================================================

  unsigned short init_id[2];
  for (size_t i = 0; i < 2; ++i) {
    init_id[i] = sizePart();
    // valence quark from a dissociated beam is attached to the outgoing CepGen remnant; an elastic parton to itself
    addCorresp(init_id[i], inel[i] ? outgoing[i]->id() : partons[i]->id());
    addParticle(pdg_id[i],
                -1,
                0,
                0,
                colour[i],
                0,
                mom_init[i].px(),
                mom_init[i].py(),
                mom_init[i].pz(),
                mom_init[i].e(),
                0.,
                0.,
                1.);
  }

  //===========================================================================================
  // outgoing valence quarks, absorbing the kT recoil of the emitted partons
  //===========================================================================================

  for (size_t i = 0; i < 2; ++i) {
    if (!inel[i])
      continue;
    const Vec4 mom_out = mom_init[i] - mom_parton[i];
    addCorresp(sizePart(), outgoing[i]->id());
    addParticle(pdg_id[i],
                1,
                init_id[i],
                0,
                colour[i],
                0,
                mom_out.px(),
                mom_out.py(),
                mom_out.pz(),
                mom_out.e(),
                std::sqrt(std::max(0., mom_out.m2Calc())),
                0.,
                1.);
  }
  return true;
}

void CepGenEvent::setProcess(int id, double cross_section, double q2_scale, double alpha_qed, double alpha_qcd) {
  LHAup::setProcess(id, cross_section, q2_scale, alpha_qed, alpha_qcd);
  py_cg_corresp_.clear();
}

unsigned short CepGenEvent::cepgenId(unsigned short pythia_id) const {
  if (py_cg_corresp_.count(pythia_id) == 0)
    return INVALID_ID;
  return py_cg_corresp_.at(pythia_id);
}

unsigned short CepGenEvent::pythiaId(unsigned short cepgen_id) const {
  const auto it = std::find_if(py_cg_corresp_.begin(), py_cg_corresp_.end(), [&cepgen_id](const auto& py_cg) {
    return py_cg.second == cepgen_id;
  });
  if (it != py_cg_corresp_.end())
    return it->first;
  return INVALID_ID;
}

void CepGenEvent::addCepGenParticle(const cepgen::Particle& part,
                                    int status,
                                    const std::pair<int, int>& mothers,
                                    const std::pair<int, int>& colours) {
  Vec4 mom_part(momToVec4(part.momentum()));
  if (collinearised_)
    mom_part.rotbst(boost_cg_to_lha_);
  int pdg_id = part.integerPdgId();
  if (status == INVALID_ID)
    switch (part.status()) {
      case cepgen::Particle::Status::Resonance:
      case cepgen::Particle::Status::Fragmented:
        status = 2;
        break;
      default: {
        if (part.pdgId() == 21 && static_cast<int>(part.status()) == 12)
          pdg_id = -21;  // workaround for HepMC2 interface
        else
          status = 1;
      } break;
    }
  addCorresp(sizePart(), part.id());
  addParticle(pdg_id,
              status,
              mothers.first,
              mothers.second,
              colours.first,
              colours.second,
              mom_part.px(),
              mom_part.py(),
              mom_part.pz(),
              mom_part.e(),
              mom_part.mCalc(),
              0.,
              0.,
              0.);
}

void CepGenEvent::addCorresp(unsigned short pythia_id, unsigned short cepgen_id) {
  py_cg_corresp_[pythia_id] = cepgen_id;
}

void CepGenEvent::dumpCorresp() const {
  CG_INFO("CepGenEvent:dump").log([&](auto& msg) {
    msg << "List of Pythia ←|→ CepGen particle ids correspondence";
    for (const auto& py_cg : py_cg_corresp_)
      msg << "\n\t" << py_cg.first << " <-> " << py_cg.second;
  });
}

std::pair<int, int> CepGenEvent::findMothers(const cepgen::Event& cepgen_event,
                                             const cepgen::Particle& cepgen_particle) const {
  std::pair out = {0, 0};

  const auto& mothers = cepgen_particle.mothers();
  if (mothers.empty())
    return out;
  const unsigned short moth1_cg_id = *mothers.begin();
  if (out.first = pythiaId(moth1_cg_id); out.first == INVALID_ID) {
    const auto& moth = cepgen_event(moth1_cg_id);
    out = {(moth.mothers().size() > 0) ? pythiaId(*moth.mothers().begin()) : 0,
           (moth.mothers().size() > 1) ? pythiaId(*moth.mothers().rbegin()) : 0};
  }
  if (mothers.size() > 1) {
    const unsigned short moth2_cg_id = *mothers.rbegin();
    out.second = pythiaId(moth2_cg_id);
    if (out.second == INVALID_ID)
      out.second = 0;
  }
  return out;
}
