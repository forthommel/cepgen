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
#include "CepGenAddOns/Pythia8Wrapper/CepGenEventInterface.h"

namespace Pythia8 {
  /// Convert a CepGen particle momentum into its Pythia8 counterpart
  Vec4 momToVec4(const cepgen::Momentum& cg_mom) {
    return Vec4(cg_mom.px(), cg_mom.py(), cg_mom.pz(), cg_mom.energy());
  }

  CepGenEventInterface::CepGenEventInterface()
      : LHAup(3), mp_(cepgen::PDG::get().mass(cepgen::PDG::proton)), mp2_(mp_ * mp_) {}

  void CepGenEventInterface::initialise(const cepgen::Beam& pos_beam, const cepgen::Beam& neg_beam) {
    setBeamA((short)pos_beam.pdgId(), pos_beam.momentum().pz());
    setBeamB((short)neg_beam.pdgId(), neg_beam.momentum().pz());
    inel1_ = !pos_beam.elastic();
    inel2_ = !neg_beam.elastic();
  }

  void CepGenEventInterface::addComments(const std::string& comments) {
#if PYTHIA_VERSION_INTEGER >= 8200
    osLHEF << comments;
#else
    CG_WARNING("CepGenEventInterface:addComments") << "Pythia 8 is too outdated... Unused comments: " << comments;
#endif
  }

  void CepGenEventInterface::setCrossSection(int id, double cross_section, double cross_section_err) {
    addProcess(0, cross_section, cross_section_err, 100.);
    setXSec(id, cross_section);
    setXErr(id, cross_section_err);
  }

  void CepGenEventInterface::feedEvent(const cepgen::Event& ev, unsigned int type) {
    const auto &parton1 = ev.oneWithRole(cepgen::Particle::Parton1),
               &parton2 = ev.oneWithRole(cepgen::Particle::Parton2);
    const auto p4_parton1 = momToVec4(parton1.momentum()), p4_parton2 = momToVec4(parton2.momentum());
    const auto scale = (parton1.momentum() + parton2.momentum()).mass();
    setProcess(0, ev.metadata.at("weight"), scale, ev.metadata.at("alphaEM"), ev.metadata.at("alphaS"));

    auto colour_index = MIN_COLOUR_INDEX;

    if (type & Type::incomingBeams) {
      const auto &ip1 = ev.oneWithRole(cepgen::Particle::IncomingBeam1),
                 &ip2 = ev.oneWithRole(cepgen::Particle::IncomingBeam2);
      const auto p4_ip1 = momToVec4(ip1.momentum()), p4_ip2 = momToVec4(ip2.momentum());
      addParticle(
          ip1.integerPdgId(), -9, 0, 0, 0, 0, p4_ip1.px(), p4_ip1.py(), p4_ip1.pz(), p4_ip1.e(), p4_ip1.mCalc(), 0, 0);
      addParticle(
          ip2.integerPdgId(), -9, 0, 0, 0, 0, p4_ip2.px(), p4_ip2.py(), p4_ip2.pz(), p4_ip2.e(), p4_ip2.mCalc(), 0, 0);
    }
    if (type & Type::partonsCollinear) {  // full event content (with collinear partons)
      const auto &op1 = ev(cepgen::Particle::Role::OutgoingBeam1)[0],
                 &op2 = ev(cepgen::Particle::Role::OutgoingBeam2)[0];
      const auto x1 = cepgen::utils::xBj(std::fabs(parton1.momentum().mass2()), mp2_, op1.momentum().mass2()),
                 x2 = cepgen::utils::xBj(std::fabs(parton2.momentum().mass2()), mp2_, op2.momentum().mass2());

      //===========================================================================================
      // incoming valence quarks
      //FIXME select quark flavours accordingly
      auto p4_iq1 = p4_parton1;
      unsigned short q1_pdg = parton1.integerPdgId();
      unsigned short q1_col = 0;
      if (inel1_) {
        q1_pdg = 2;
        q1_col = colour_index++;
        p4_iq1 = momToVec4(x1 * ev.oneWithRole(cepgen::Particle::IncomingBeam1).momentum());
      }
      const auto q1_id = sizePart();
      addCorresp(q1_id, op1.id());
      addParticle(q1_pdg, -1, 0, 0, q1_col, 0, p4_iq1.px(), p4_iq1.py(), p4_iq1.pz(), p4_iq1.e(), p4_iq1.mCalc(), 0, 1);

      auto p4_iq2 = p4_parton2;
      unsigned short q2_pdg = parton2.integerPdgId();
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
      //setIdX(parton1.integerPdgId(), parton2.integerPdgId(), x1, x2);
      setPdf(q1_pdg, q2_pdg, x1, x2, scale, 0., 0., false);  // flavour/x value of hard-process initiators

      //===========================================================================================
      // outgoing valence quarks
      if (inel1_) {
        const Vec4 p4_oq1 = p4_iq1 - p4_parton1;
        addParticle(
            q1_pdg, 1, q1_id, q2_id, q1_col, 0, p4_oq1.px(), p4_oq1.py(), p4_oq1.pz(), p4_oq1.e(), p4_oq1.mCalc(), 0, 1);
      }
      if (inel2_) {
        const Vec4 p4_oq2 = p4_iq2 - p4_parton2;
        addParticle(
            q2_pdg, 1, q1_id, q2_id, q2_col, 0, p4_oq2.px(), p4_oq2.py(), p4_oq2.pz(), p4_oq2.e(), p4_oq2.mCalc(), 0, 1);
      }
      //===========================================================================================
    } else if (type & Type::partonsKT) {  // addition of incoming partons (possibly with kT)
      addCepGenParticle(parton1, -2);
      addCepGenParticle(parton2, -2);
    }
    if (type & Type::beamRemnants) {  // full beam remnants content
      if (!inel1_)
        for (const auto& p : ev(cepgen::Particle::Role::OutgoingBeam1))
          addCepGenParticle(p, INVALID_ID, findMothers(ev, p));
      if (!inel2_)
        for (const auto& p : ev(cepgen::Particle::Role::OutgoingBeam2))
          addCepGenParticle(p, INVALID_ID, findMothers(ev, p));
    }

    //=============================================================================================
    // central system
    if (type & Type::central) {
      const unsigned short central_colour = colour_index++;
      for (const auto& p : ev(cepgen::Particle::Role::CentralSystem)) {
        std::pair<int, int> colours = {0, 0}, mothers = {1, 2};
        if (!(type & Type::partonsCollinear))
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
        if (type & Type::beamRemnants && p.status() == cepgen::Particle::Status::Resonance)
          status = 2;
        addCepGenParticle(p, status, mothers, colours);
      }
    }
    //=============================================================================================
  }

  void CepGenEventInterface::setProcess(
      int id, double cross_section, double q2_scale, double alpha_qed, double alpha_qcd) {
    LHAup::setProcess(id, cross_section, q2_scale, alpha_qed, alpha_qcd);
    py_cg_corresp_.clear();
  }

  unsigned short CepGenEventInterface::cepgenId(unsigned short py_id) const {
    if (py_cg_corresp_.count(py_id) > 0)
      return py_cg_corresp_.at(py_id);
    return INVALID_ID;
  }

  unsigned short CepGenEventInterface::pythiaId(unsigned short cg_id) const {
    auto it = std::find_if(
        py_cg_corresp_.begin(), py_cg_corresp_.end(), [&cg_id](const auto& py_cg) { return py_cg.second == cg_id; });
    if (it != py_cg_corresp_.end())
      return it->first;
    return INVALID_ID;
  }

  void CepGenEventInterface::addCepGenParticle(const cepgen::Particle& cg_part,
                                               int status,
                                               const std::pair<int, int>& mothers,
                                               const std::pair<int, int>& colours) {
    const auto p4 = momToVec4(cg_part.momentum());
    int pdgid = cg_part.integerPdgId();
    if (status == INVALID_ID) {
      if (cg_part.status() == cepgen::Particle::Status::Resonance ||
          cg_part.status() == cepgen::Particle::Status::Fragmented)
        status = 2;
      else if (cg_part.pdgId() == cepgen::PDG::gluon && (int)cg_part.status() == 12)
        pdgid = -21;  // workaround for HepMC2 interface
      else
        status = 1;
    }
    const auto &moth1 = mothers.first, &moth2 = mothers.second;
    const auto &col1 = colours.first, &col2 = colours.second;
    const auto tau = 0., spin = 0., scale = 0.;
    addCorresp(sizePart(), cg_part.id());
    addParticle(
        pdgid, status, moth1, moth2, col1, col2, p4.px(), p4.py(), p4.pz(), p4.e(), p4.mCalc(), tau, spin, scale);
  }

  void CepGenEventInterface::addCorresp(unsigned short py_id, unsigned short cg_id) { py_cg_corresp_[py_id] = cg_id; }

  void CepGenEventInterface::dumpCorresp(std::ostream& os) const {
    os << "List of Pythia ←|→ CepGen particle ids correspondence";
    for (const auto& py_cg : py_cg_corresp_)
      os << "\n\t" << py_cg.first << " <-> " << py_cg.second;
  }

  cepgen::Particle& CepGenEventInterface::addParticleToEvent(cepgen::Event& cg_event,
                                                             cepgen::Particle::Role cg_role,
                                                             const Particle& py_part,
                                                             const Vec4& py_mom) const {
    checkPDGid(py_part);
    auto& op = cg_event.addParticle(cg_role).get();  // add the particle to the event content
    op.setPdgId((long)py_part.id());
    op.setStatus(py_part.isFinal()                                  ? cepgen::Particle::Status::FinalState
                 : cg_role == cepgen::Particle::Role::CentralSystem ? cepgen::Particle::Status::Propagator
                                                                    : cepgen::Particle::Status::Fragmented);
    op.setMomentum(cepgen::Momentum(py_mom.px(), py_mom.py(), py_mom.pz(), py_mom.e()).setMass(py_mom.mCalc()));
    //FIXME addCorresp(py_part.index() - offset_, op.id());
    return op;
  }

  void CepGenEventInterface::updateEvent(const Pythia8::Event& py_event,
                                         cepgen::Event& cg_event,
                                         double& weight,
                                         bool correct_central) const {
    std::vector<unsigned short> central_parts;
    py_event.list();
    if (first_evt_) {
      offset_ = 0;
      for (unsigned short i = 1; i < py_event.size(); ++i)
        if (py_event[i].status() == -PYTHIA_STATUS_IN_BEAM)  // no incoming particles in further stages
          offset_++;
      first_evt_ = false;
    }

    for (unsigned short i = 1 + offset_; i < py_event.size(); ++i) {
      const auto& p = py_event[i];
      const unsigned short cg_id = cepgenId(i - offset_);
      if (cg_id != INVALID_ID) {  // particle already in the event
        auto& cg_part = cg_event[cg_id];
        //--- fragmentation result
        if (cg_part.role() == cepgen::Particle::Role::OutgoingBeam1 ||
            cg_part.role() == cepgen::Particle::Role::OutgoingBeam2) {
          cg_part.setStatus(cepgen::Particle::Status::Fragmented);
          continue;
        }
        if (cg_part.role() == cepgen::Particle::Role::CentralSystem && p.status() < 0) {  // resonance decayed
          weight *= p.particleDataEntry().pickChannel().bRatio();  // apply branching ratio for this decay
          cg_part.setStatus(cepgen::Particle::Status::Resonance);
          central_parts.emplace_back(i);
        }
        if (p.idAbs() != abs(cg_part.integerPdgId())) {  // particle is not what we expect
          CG_INFO("CepGenEventInterface:update") << "LHAEVT event content:";
          const_cast<CepGenEventInterface&>(*this).listEvent();
          CG_INFO("CepGenEventInterface:update") << "Pythia event content:";
          py_event.list();
          CG_INFO("CepGenEventInterface:update") << "CepGen event content:";
          cg_event.dump();
          CG_INFO("CepGenEventInterface:update").log([this](auto& log) {
            log << "Correspondence:";
            dumpCorresp(log.stream());
          });

          throw CG_FATAL("CepGenEventInterface:update")
              << "Event list corruption detected for (Pythia/CepGen) particle " << i << "/" << cg_id << ":\n\t"
              << "should be " << abs(p.id()) << ", "
              << "got " << cg_part.integerPdgId() << "!";
        }
      }
      //--- check for messed up particles parentage and discard incoming beam particles
      /*else if ( p.mother1() > i || p.mother1() <= offset_ )
           continue;
        else if ( p.mother2() > i || p.mother2() <= offset_ )
         continue;*/
      else {
        //----- new particle to be added
        const auto role = findRole(cg_event, py_event, p);
        switch (role) {
          case cepgen::Particle::Role::OutgoingBeam1:
            cg_event[cepgen::Particle::Role::OutgoingBeam1][0].get().setStatus(cepgen::Particle::Status::Fragmented);
            break;
          case cepgen::Particle::Role::OutgoingBeam2:
            cg_event[cepgen::Particle::Role::OutgoingBeam2][0].get().setStatus(cepgen::Particle::Status::Fragmented);
            break;
          default:
            break;
        }
        // found the role ; now we can add the particle
        auto& cg_part = addParticleToEvent(cg_event, role, p, p.p());
        if (correct_central && role == cepgen::Particle::Role::CentralSystem) {
          const auto& ip = std::find(central_parts.begin(), central_parts.end(), p.mother1());
          if (ip != central_parts.end())
            cg_part.setMomentum(cg_event[cepgenId(*ip - offset_)].momentum());
        }
        for (const auto& moth_id : p.motherList()) {
          if (moth_id <= offset_)
            continue;
          const unsigned short moth_cg_id = cepgenId(moth_id - offset_);
          if (moth_cg_id != INVALID_ID)
            cg_part.addMother(cg_event[moth_cg_id]);
          else
            cg_part.addMother(addParticleToEvent(cg_event, role, py_event[moth_id], p.p()));
          if (!p.isFinal()) {
            if (p.isResonance() || !p.daughterList().empty())
              cg_part.setStatus(cepgen::Particle::Status::Resonance);
            else
              cg_part.setStatus(cepgen::Particle::Status::Undefined);
          }
        }
      }
    }
  }

  void CepGenEventInterface::checkPDGid(const Particle& py_part) const {
    if (cepgen::PDG::get().has(py_part.id()))
      return;
    cepgen::ParticleProperties prop;
    prop.pdgid = py_part.idAbs();
    prop.name = py_part.name();
    prop.descr = py_part.name();
    prop.colours = py_part.col();  // colour factor
    prop.mass = py_part.m0();
    prop.width = py_part.mWidth();
    prop.charge = py_part.charge();  // charge
    prop.fermion = py_part.isLepton();
    cepgen::PDG::get().define(prop);
  }

  cepgen::Particle::Role CepGenEventInterface::findRole(const cepgen::Event& cg_event,
                                                        const Event& py_event,
                                                        const Particle& py_part) const {
    for (const auto& par_id : py_part.motherList()) {
      if (par_id == 1 && offset_ > 0)
        return cepgen::Particle::Role::OutgoingBeam1;
      if (par_id == 2 && offset_ > 0)
        return cepgen::Particle::Role::OutgoingBeam2;
      const auto par_cg_id = cepgenId(par_id - offset_);
      if (par_cg_id != INVALID_ID)
        return cg_event[par_cg_id].role();
      return findRole(cg_event, py_event, py_event[par_id]);
    }
    return cepgen::Particle::Role::UnknownRole;
  }

  std::pair<int, int> CepGenEventInterface::findMothers(const cepgen::Event& cg_event,
                                                        const cepgen::Particle& cg_part) const {
    std::pair<int, int> out = {0, 0};

    const auto& mothers = cg_part.mothers();
    if (mothers.empty())
      return out;
    const unsigned short moth1_cg_id = *mothers.begin();
    out.first = pythiaId(moth1_cg_id);
    if (out.first == INVALID_ID) {
      const auto& moth = cg_event[moth1_cg_id];
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
