/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2023  Laurent Forthomme
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

#include "CepGen/Event/Event.h"
#include "CepGen/Modules/ProcessFactory.h"
#include "CepGen/Physics/Constants.h"
#include "CepGen/Physics/HeavyIon.h"
#include "CepGen/Physics/PDG.h"
#include "CepGen/Process/KTProcess.h"
#include "CepGen/Utils/Message.h"
#include "CepGenAddOns/PythonWrapper/Environment.h"
#include "CepGenAddOns/PythonWrapper/Error.h"
#include "CepGenAddOns/PythonWrapper/PythonUtils.h"

namespace cepgen {
  namespace proc {
    /// Compute the matrix element for a generic \f$k_{\rm T}\f$-factorised process defined in a Python weighting function
    class PythonKTProcess final : public KTProcess {
    public:
      explicit PythonKTProcess(const ParametersList& params) : KTProcess(params, {PDG::muon, PDG::muon}) { init(); }

      void init() {
        const auto mod_name = steer<std::string>("name");
        const auto class_name = steer<std::string>("class").empty() ? mod_name : steer<std::string>("class");
        const auto mod = python::importModule("Process." + mod_name);  // new
        if (!mod)
          throw PY_ERROR << "Failed to load module '" << mod_name << "'.";
        const auto obj = python::getAttribute(mod, class_name);
        if (!obj)
          throw PY_ERROR << "Failed to retrieve class '" << class_name << "'.";
        obj_ = python::call(obj);
        if (!obj_)
          throw PY_ERROR << "Failed to instantiate new class '" << class_name << "'.";
        python::setAttribute(obj_, "mp", Process::mp_);
        python::setAttribute(obj_, "units", constants::GEVM2_TO_PB);
        python::setAttribute(obj_, "processParameters", params_);
      }
      ProcessPtr clone() const override {
        auto proc = new PythonKTProcess(params_);
        proc->init();
        return ProcessPtr(proc);
      }

    private:
      void preparePhaseSpace() override;
      double computeKTFactorisedMatrixElement() override {
        return python::get<double>(python::callMethod(
            obj_, "__eval__", qt1_, qt2_, phi_qt1_, phi_qt2_, y1_, y2_, pt_diff_, phi_pt_diff_, mX(), mY()));
      }
      void fillCentralParticlesKinematics() override;

      python::Environment env_;
      python::ObjectPtr obj_;
      double y1_;           ///< First outgoing particle rapidity
      double y2_;           ///< Second outgoing particle rapidity
      double pt_diff_;      ///< Transverse momentum balance between outgoing particles
      double phi_pt_diff_;  ///< Azimuthal angle difference between outgoing particles
    };

    void PythonKTProcess::preparePhaseSpace() {
      defineVariable(y1_,
                     Mapping::linear,
                     kinematics().cuts().central.rapidity_single,
                     {-6., 6.},
                     "First central particle rapidity");
      defineVariable(y2_,
                     Mapping::linear,
                     kinematics().cuts().central.rapidity_single,
                     {-6., 6.},
                     "Second central particle rapidity");
      defineVariable(pt_diff_,
                     Mapping::linear,
                     kinematics().cuts().central.pt_diff,
                     {0., 50.},
                     "Transverse momentum difference between central particles");
      defineVariable(phi_pt_diff_,
                     Mapping::linear,
                     kinematics().cuts().central.phi_diff,
                     {0., 2. * M_PI},
                     "Central particles azimuthal angle difference");

      //===========================================================================================
      // feed phase space cuts to the common block
      //===========================================================================================

      // export the limits into external variables
      /*auto save_lim = [](const Limits& lim, int& on, double& min, double& max) {
        on = lim.valid();
        min = max = 0.;
        if (lim.hasMin())
          min = lim.min();
        max = lim.hasMax() ? lim.max() : 9999.999;
      };

      save_lim(kinematics().cuts().central.pt_single, kincuts_.ipt, kincuts_.pt_min, kincuts_.pt_max);
      save_lim(kinematics().cuts().central.energy_single, kincuts_.iene, kincuts_.ene_min, kincuts_.ene_max);
      save_lim(kinematics().cuts().central.eta_single, kincuts_.ieta, kincuts_.eta_min, kincuts_.eta_max);
      save_lim(kinematics().cuts().central.mass_sum, kincuts_.iinvm, kincuts_.invm_min, kincuts_.invm_max);
      save_lim(kinematics().cuts().central.pt_sum, kincuts_.iptsum, kincuts_.ptsum_min, kincuts_.ptsum_max);
      save_lim(kinematics().cuts().central.rapidity_diff, kincuts_.idely, kincuts_.dely_min, kincuts_.dely_max);*/

      //===========================================================================================
      // feed run parameters to the common block
      //===========================================================================================

      //genparams_.icontri = (int)kinematics().incomingBeams().mode();

      //-------------------------------------------------------------------------------------------
      // incoming beams information
      //-------------------------------------------------------------------------------------------

      //--- positive-z incoming beam
      /*genparams_.inp1 = kinematics().incomingBeams().positive().momentum().pz();
      //--- check if first incoming beam is a heavy ion
      if (HeavyIon::isHI(kinematics().incomingBeams().positive().pdgId())) {
        const auto in1 = HeavyIon::fromPdgId(kinematics().incomingBeams().positive().pdgId());
        genparams_.a_nuc1 = in1.A;
        genparams_.z_nuc1 = (unsigned short)in1.Z;
        if (genparams_.z_nuc1 > 1) {
          event().oneWithRole(Particle::IncomingBeam1).setPdgId((pdgid_t)in1);
          event().oneWithRole(Particle::OutgoingBeam1).setPdgId((pdgid_t)in1);
        }
      } else
        genparams_.a_nuc1 = genparams_.z_nuc1 = 1;

      //--- negative-z incoming beam
      genparams_.inp2 = kinematics().incomingBeams().negative().momentum().pz();
      //--- check if second incoming beam is a heavy ion
      if (HeavyIon::isHI(kinematics().incomingBeams().negative().pdgId())) {
        const auto in2 = HeavyIon::fromPdgId(kinematics().incomingBeams().negative().pdgId());
        genparams_.a_nuc2 = in2.A;
        genparams_.z_nuc2 = (unsigned short)in2.Z;
        if (genparams_.z_nuc2 > 1) {
          event().oneWithRole(Particle::IncomingBeam2).setPdgId((pdgid_t)in2);
          event().oneWithRole(Particle::OutgoingBeam2).setPdgId((pdgid_t)in2);
        }
      } else
        genparams_.a_nuc2 = genparams_.z_nuc2 = 1;*/

      //-------------------------------------------------------------------------------------------
      // intermediate partons information
      //-------------------------------------------------------------------------------------------

      //FIXME
      //genparams_.iflux1 = (int)kinematics().incomingBeams().positive().ktFlux();
      //genparams_.iflux2 = (int)kinematics().incomingBeams().negative().ktFlux();
    }

    void PythonKTProcess::fillCentralParticlesKinematics() {
      //===========================================================================================
      // outgoing beam remnants
      //===========================================================================================

      /*pX() = Momentum(evtkin_.px);
      pY() = Momentum(evtkin_.py);
      // express these momenta per nucleon
      pX() *= 1. / genparams_.a_nuc1;
      pY() *= 1. / genparams_.a_nuc2;

      //===========================================================================================
      // intermediate partons
      //===========================================================================================

      q1() = pA() - pX();
      q2() = pB() - pY();
      event().oneWithRole(Particle::Intermediate).setMomentum(q1() + q2());

      //===========================================================================================
      // central system
      //===========================================================================================

      auto oc = event()[Particle::CentralSystem];  // retrieve all references
                                                   // to central system particles
      for (int i = 0; i < evtkin_.nout; ++i) {
        auto& p = oc[i].get();  // retrieve a reference to the specific particle
        p.setPdgId((long)evtkin_.pdg[i]);
        p.setStatus(Particle::Status::FinalState);
        p.setMomentum(Momentum(evtkin_.pc[i]));
      }*/
    }
  }  // namespace proc
}  // namespace cepgen

// register process
REGISTER_PROCESS("python", PythonKTProcess);