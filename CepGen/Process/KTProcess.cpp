/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2016-2023  Laurent Forthomme
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
#include "CepGen/KTFluxes/KTFlux.h"
#include "CepGen/Modules/PartonFluxFactory.h"
#include "CepGen/Physics/Beam.h"
#include "CepGen/Process/KTProcess.h"

namespace cepgen {
  namespace proc {
    KTProcess::KTProcess(const ParametersList& params, const pdgids_t& central)
        : Process(params), produced_parts_(central) {
      event().map()[Particle::CentralSystem].resize(central.size());
    }

    void KTProcess::addEventContent() {
      Process::setEventContent(
          {// incoming state
           {Particle::IncomingBeam1, kinematics().incomingBeams().positive().pdgId()},
           {Particle::IncomingBeam2, kinematics().incomingBeams().negative().pdgId()},
           {Particle::Parton1, kinematics().incomingBeams().positive().daughterId()},
           {Particle::Parton2, kinematics().incomingBeams().negative().daughterId()}},
          {// outgoing state
           {Particle::OutgoingBeam1, {kinematics().incomingBeams().positive().pdgId()}},
           {Particle::OutgoingBeam2, {kinematics().incomingBeams().negative().pdgId()}},
           {Particle::CentralSystem, produced_parts_}});
      setExtraContent();
      CG_DEBUG("KTProcess:addEventContent") << "Addition of:\n\t"
                                            << "Intermediate partons: "
                                            << pdgids_t{kinematics().incomingBeams().positive().daughterId(),
                                                        kinematics().incomingBeams().negative().daughterId()}
                                            << "\n\t"
                                            << "Produced system: " << produced_parts_ << ".\n\t" << event();
    }

    void KTProcess::prepareBeams() {
      auto set_beam_properties = [](Beam& beam) {
        ParametersList params;
        switch (beam.mode()) {
          case Beam::Mode::HIElastic:
            params.set<ParametersList>("partonFlux",
                                       PartonFluxFactory::get().describeParameters("ElasticHeavyIonKT").parameters());
            break;
          case Beam::Mode::ProtonInelastic:
            params.set<ParametersList>("partonFlux",
                                       PartonFluxFactory::get().describeParameters("BudnevInelasticKT").parameters());
            break;
          case Beam::Mode::ProtonElastic:
            params.set<ParametersList>("partonFlux",
                                       PartonFluxFactory::get().describeParameters("BudnevElasticKT").parameters());
            break;
          default:
            throw CG_FATAL("KTProcess:prepareBeams")
                << "No fallback strategy defined for beam mode '" << beam.mode() << "'.";
        }
        beam.setParameters(params);
      };
      set_beam_properties(kinematics().incomingBeams().positive());
      set_beam_properties(kinematics().incomingBeams().negative());
    }

    void KTProcess::prepareKinematics() {
      if (!kinematics().incomingBeams().positive().flux().ktFactorised() ||
          !kinematics().incomingBeams().negative().flux().ktFactorised())
        throw CG_FATAL("KTProcess:prepareKinematics") << "Invalid incoming parton fluxes.";

      const auto& ib1 = event().oneWithRole(Particle::IncomingBeam1);
      const auto& ib2 = event().oneWithRole(Particle::IncomingBeam2);

      //============================================================================================
      // register the incoming partons' variables
      //============================================================================================

      const auto log_lim_kt = utils::log(kinematics().cuts().initial.qt).truncate(Limits{-10., 10.});
      defineVariable(qt1_, Mapping::exponential, log_lim_kt, "First incoming parton virtuality");
      defineVariable(qt2_, Mapping::exponential, log_lim_kt, "Second incoming parton virtuality");

      const auto lim_phi_kt = utils::log(kinematics().cuts().initial.phi_qt).truncate(Limits{0., 2. * M_PI});
      defineVariable(phi_qt1_, Mapping::linear, lim_phi_kt, "First incoming parton azimuthal angle");
      defineVariable(phi_qt2_, Mapping::linear, lim_phi_kt, "Second incoming parton azimuthal angle");

      //============================================================================================
      // register all process-dependent variables
      //============================================================================================

      preparePhaseSpace();

      //============================================================================================
      // register the outgoing remnants' variables
      //============================================================================================

      mX2() = ib1.mass2();
      mY2() = ib2.mass2();
      if (kinematics().incomingBeams().positive().fragmented())
        defineVariable(
            mX2(), Mapping::square, kinematics().cuts().remnants.mx, "Positive z proton remnant squared mass");
      if (kinematics().incomingBeams().negative().fragmented())
        defineVariable(
            mY2(), Mapping::square, kinematics().cuts().remnants.mx, "Negative z proton remnant squared mass");
    }

    double KTProcess::computeWeight() {
      const auto cent_me = computeKTFactorisedMatrixElement();
      if (cent_me <= 0)
        return 0.;  // avoid computing the fluxes if the matrix element is already null

      // convolute with fluxes according to modelling specified in parameters card
      const auto& flux1 = dynamic_cast<const KTFlux&>(kinematics().incomingBeams().positive().flux());
      const auto& flux2 = dynamic_cast<const KTFlux&>(kinematics().incomingBeams().negative().flux());

      return flux1.fluxMX2(x1_, qt1_ * qt1_, mX2()) * M_1_PI * flux2.fluxMX2(x2_, qt2_ * qt2_, mY2()) * M_1_PI *
             cent_me;
    }

    void KTProcess::fillKinematics(bool) {
      fillCentralParticlesKinematics();  // process-dependent!

      // set the kinematics of the incoming and outgoing beams (or remnants)
      const auto& ib1 = event().oneWithRole(Particle::IncomingBeam1);
      const auto& ib2 = event().oneWithRole(Particle::IncomingBeam2);
      auto& ob1 = event().oneWithRole(Particle::OutgoingBeam1);
      auto& ob2 = event().oneWithRole(Particle::OutgoingBeam2);
      auto& p1 = event().oneWithRole(Particle::Parton1);
      auto& p2 = event().oneWithRole(Particle::Parton2);
      auto& cm = event().oneWithRole(Particle::Intermediate);

      // beam systems
      if (kinematics().incomingBeams().positive().fragmented())
        ob1.setMass(mX());
      if (kinematics().incomingBeams().negative().fragmented())
        ob2.setMass(mY());

      // parton systems
      p1.setMomentum(ib1.momentum() - pX(), true);
      p2.setMomentum(ib2.momentum() - pY(), true);

      // two-parton system
      cm.setMomentum(p1.momentum() + p2.momentum());
    }

    ParametersDescription KTProcess::description() {
      auto desc = Process::description();
      desc.setDescription("Unnamed kT-factorised process");
      return desc;
    }
  }  // namespace proc
}  // namespace cepgen
