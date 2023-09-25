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
#include "CepGen/Physics/PDG.h"
#include "CepGen/Process/FactorisedProcess.h"
#include "CepGen/Utils/Message.h"
#include "CepGenAddOns/PythonWrapper/Environment.h"
#include "CepGenAddOns/PythonWrapper/Error.h"
#include "CepGenAddOns/PythonWrapper/PythonUtils.h"

namespace cepgen {
  namespace proc {
    /// Compute the matrix element for a generic factorised process defined in a Python weighting function
    class PythonProcess final : public FactorisedProcess {
    public:
      explicit PythonProcess(const ParametersList& params) : FactorisedProcess(params, {PDG::muon, PDG::muon}) {
        initialise();
      }

      void initialise() {
        const auto mod_name = steer<std::string>("object");
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
        auto proc = new PythonProcess(params_);
        proc->initialise();
        return ProcessPtr(proc);
      }

      static ParametersDescription description() {
        auto desc = FactorisedProcess::description();
        desc.setDescription("Python generic process");
        return desc;
      }

    private:
      void prepareFactorisedPhaseSpace() override;
      double computeFactorisedMatrixElement() override {
        ParametersList args;
        for (size_t i = 0; i < vars_.size(); ++i)
          args.set<double>(vars_names_.at(i), vars_.at(i));
        python::setAttribute(obj_, "variable", args);
        return python::get<double>(python::callMethod(obj_, "__eval__"));
      }
      void fillCentralParticlesKinematics() override;

      python::Environment env_;
      python::ObjectPtr obj_;
      std::vector<double> vars_;
      std::vector<std::string> vars_names_;
    };

    void PythonProcess::prepareFactorisedPhaseSpace() {
      auto vars_def = python::getAttribute(obj_, "_variables");
      auto vars_coll = python::get<ParametersList>(vars_def);
      for (auto& var : vars_coll.keys()) {
        const auto& var_attrs = vars_coll.get<ParametersList>(var);
        const auto& type = var_attrs.get<std::string>("type");
        Mapping mapping;
        if (type == "linear")
          mapping = Mapping::linear;
        else if (type == "square")
          mapping = Mapping::square;
        else if (type == "exponential")
          mapping = Mapping::exponential;
        else if (type == "power_law")
          mapping = Mapping::power_law;
        else
          throw CG_FATAL("PythonProcess")
              << "Invalid mapping type for '" << var << "': '" << type << "' is not supported.";
        auto& var_ref = vars_.emplace_back(0.);
        vars_names_.emplace_back(var);
        defineVariable(var_ref, mapping, var_attrs.get<Limits>("range"), var);
      }
    }

    void PythonProcess::fillCentralParticlesKinematics() {
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
REGISTER_PROCESS("python", PythonProcess);
