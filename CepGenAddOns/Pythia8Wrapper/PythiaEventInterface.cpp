/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2020-2023  Laurent Forthomme
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
#include "CepGen/Event/Particle.h"
#include "CepGen/Physics/Beam.h"
#include "CepGen/Physics/PDG.h"
#include "CepGen/Physics/Utils.h"
#include "CepGenAddOns/Pythia8Wrapper/PythiaEventInterface.h"

namespace Pythia8 {
  /// Convert a CepGen particle momentum into its Pythia8 counterpart
  Vec4 momToVec4(const cepgen::Momentum& mom) { return Vec4(mom.px(), mom.py(), mom.pz(), mom.energy()); }

  CepGenEvent::CepGenEvent() : LHAup(3), mp_(cepgen::PDG::get().mass(cepgen::PDG::proton)), mp2_(mp_ * mp_) {}

  void CepGenEvent::initialise(const cepgen::Beam& pos_beam, const cepgen::Beam& neg_beam) {
    setBeamA((short)pos_beam.pdgId(), pos_beam.momentum().pz());
    setBeamB((short)neg_beam.pdgId(), neg_beam.momentum().pz());
    inel1_ = !pos_beam.elastic();
    inel2_ = !neg_beam.elastic();
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
  }

  void CepGenEvent::feedEvent(const cepgen::Event& ev, unsigned int type) {
    const auto &part1 = ev.oneWithRole(cepgen::Particle::Parton1), &part2 = ev.oneWithRole(cepgen::Particle::Parton2);
    const auto p4_part1 = momToVec4(part1.momentum()), p4_part2 = momToVec4(part2.momentum());
    const double scale = (part1.momentum() + part2.momentum()).mass();
    setProcess(0, ev.metadata.at("weight"), scale, ev.metadata.at("alphaEM"), ev.metadata.at("alphaS"));

    unsigned short colour_index = MIN_COLOUR_INDEX;

    if (type & Type::partonsCollinear) {  // full event content (with collinear partons)
      const auto &op1 = ev(cepgen::Particle::OutgoingBeam1)[0], &op2 = ev(cepgen::Particle::OutgoingBeam2)[0];
      const double x1 = cepgen::utils::xBj(std::fabs(part1.momentum().mass2()), mp2_, op1.momentum().mass2()),
                   x2 = cepgen::utils::xBj(std::fabs(part2.momentum().mass2()), mp2_, op2.momentum().mass2());

      //===========================================================================================
      // incoming valence quarks
      //FIXME select quark flavours accordingly
      auto p4_iq1 = p4_part1;
      unsigned short q1_pdg = part1.integerPdgId();
      unsigned short q1_col = 0;
      if (inel1_) {
        q1_pdg = 2;
        q1_col = colour_index++;
        p4_iq1 = momToVec4(x1 * ev.oneWithRole(cepgen::Particle::IncomingBeam1).momentum());
      }
      const auto q1_id = sizePart();
      addCorresp(q1_id, op1.id());
      addParticle(q1_pdg, -1, 0, 0, q1_col, 0, p4_iq1.px(), p4_iq1.py(), p4_iq1.pz(), p4_iq1.e(), p4_iq1.mCalc(), 0, 1);

      auto p4_iq2 = p4_part2;
      unsigned short q2_pdg = part2.integerPdgId();
      unsigned short q2_col = 0;
      if (inel2_) {
        q2_pdg = 2;
        q2_col = colour_index++;
        p4_iq2 = momToVec4(x2 * ev.oneWithRole(cepgen::Particle::IncomingBeam2).momentum());
      }
      const auto q2_id = sizePart();
      addCorresp(q2_id, op2.id());
      addParticle(q2_pdg, -1, 0, 0, q2_col, 0, p4_iq2.px(), p4_iq2.py(), p4_iq2.pz(), p4_iq2.e(), p4_iq2.mCalc(), 0, 1);
      //===========================================================================================

      setIdX(q1_pdg, q2_pdg, x1, x2);
      //setIdX(part1.integerPdgId(), part2.integerPdgId(), x1, x2);
      setPdf(q1_pdg, q2_pdg, x1, x2, scale, 0., 0., false);  // flavour/x value of hard-process initiators

      //===========================================================================================
      // outgoing valence quarks
      if (inel1_) {
        const Vec4 p4_oq1 = p4_iq1 - p4_part1;
        addParticle(
            q1_pdg, 1, q1_id, q2_id, q1_col, 0, p4_oq1.px(), p4_oq1.py(), p4_oq1.pz(), p4_oq1.e(), p4_oq1.mCalc(), 0, 1);
      }
      if (inel2_) {
        const Vec4 p4_oq2 = p4_iq2 - p4_part2;
        addParticle(
            q2_pdg, 1, q1_id, q2_id, q2_col, 0, p4_oq2.px(), p4_oq2.py(), p4_oq2.pz(), p4_oq2.e(), p4_oq2.mCalc(), 0, 1);
      }
      //===========================================================================================
    } else if (type & Type::partonsKT) {  // addition of incoming partons (possibly with kT)
      addCepGenParticle(part1, -2);
      addCepGenParticle(part2, -2);
    }
    if (type & Type::beamRemnants)  // full beam remnants content
      for (const auto& syst : {cepgen::Particle::OutgoingBeam1, cepgen::Particle::OutgoingBeam2})
        for (const auto& p : ev(syst))
          addCepGenParticle(p, INVALID_ID, findMothers(ev, p));

    //=============================================================================================
    // central system
    if (type & Type::central) {
      const unsigned short central_colour = colour_index++;
      for (const auto& p : ev(cepgen::Particle::CentralSystem)) {
        std::pair<int, int> colours = {0, 0}, mothers = {1, 2};
        if (!(type & Type::partonsCollinear))
          mothers = findMothers(ev, p);
        try {
          if (cepgen::PDG::get().colours(p.pdgId()) > 1) {
            if (p.integerPdgId() > 0)  //--- particle
              colours.first = central_colour;
            else  //--- anti-particle
              colours.second = central_colour;
          }
        } catch (const cepgen::Exception&) {
        }
        int status = 1;
        if (type & Type::beamRemnants && p.status() == cepgen::Particle::Status::Resonance)
          status = 2;
        addCepGenParticle(p, status, mothers, colours);
      }
    }
    //=============================================================================================
  }

  void CepGenEvent::setProcess(int id, double cross_section, double q2_scale, double alpha_qed, double alpha_qcd) {
    LHAup::setProcess(id, cross_section, q2_scale, alpha_qed, alpha_qcd);
    py_cg_corresp_.clear();
  }

  unsigned short CepGenEvent::cepgenId(unsigned short py_id) const {
    if (py_cg_corresp_.count(py_id) > 0)
      return py_cg_corresp_.at(py_id);
    return INVALID_ID;
  }

  unsigned short CepGenEvent::pythiaId(unsigned short cg_id) const {
    auto it = std::find_if(
        py_cg_corresp_.begin(), py_cg_corresp_.end(), [&cg_id](const auto& py_cg) { return py_cg.second == cg_id; });
    if (it != py_cg_corresp_.end())
      return it->first;
    return INVALID_ID;
  }

  void CepGenEvent::addCepGenParticle(const cepgen::Particle& part,
                                      int status,
                                      const std::pair<int, int>& mothers,
                                      const std::pair<int, int>& colours) {
    const auto p4 = momToVec4(part.momentum());
    int pdgid = part.integerPdgId();
    if (status == INVALID_ID) {
      if (part.status() == cepgen::Particle::Status::Resonance || part.status() == cepgen::Particle::Status::Fragmented)
        status = 2;
      else if (part.pdgId() == cepgen::PDG::gluon && (int)part.status() == 12)
        pdgid = -21;  // workaround for HepMC2 interface
      else
        status = 1;
    }
    const auto &moth1 = mothers.first, &moth2 = mothers.second;
    const auto &col1 = colours.first, &col2 = colours.second;
    const auto tau = 0., spin = 0., scale = 0.;
    addCorresp(sizePart(), part.id());
    addParticle(
        pdgid, status, moth1, moth2, col1, col2, p4.px(), p4.py(), p4.pz(), p4.e(), p4.mCalc(), tau, spin, scale);
  }

  void CepGenEvent::addCorresp(unsigned short py_id, unsigned short cg_id) { py_cg_corresp_[py_id] = cg_id; }

  void CepGenEvent::dumpCorresp(std::ostream& os) const {
    os << "List of Pythia ←|→ CepGen particle ids correspondence";
    for (const auto& py_cg : py_cg_corresp_)
      os << "\n\t" << py_cg.first << " <-> " << py_cg.second;
  }

  std::pair<int, int> CepGenEvent::findMothers(const cepgen::Event& ev, const cepgen::Particle& p) const {
    std::pair<int, int> out = {0, 0};

    const auto& mothers = p.mothers();
    if (mothers.empty())
      return out;
    const unsigned short moth1_cg_id = *mothers.begin();
    out.first = pythiaId(moth1_cg_id);
    if (out.first == INVALID_ID) {
      const auto& moth = ev[moth1_cg_id];
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
}  // namespace Pythia8
