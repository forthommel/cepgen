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

#include <marty.h>

#include "CepGen/Core/Exception.h"
#include "CepGen/Event/Event.h"
#include "CepGen/Modules/ProcessFactory.h"
#include "CepGen/Physics/PDG.h"
#include "CepGen/Process/Process2to4.h"
#include "CepGenAddOns/MartyWrapper/MartyModelFactory.h"

namespace cepgen {
  namespace proc {
    /// Let Marty compute the matrix element for a generic 2-to-4 kT-factorised process
    class MartyProcess final : public Process2to4 {
    public:
      explicit MartyProcess(const ParametersList& params)
          : Process2to4(params,
                        {0, 0},
                        0),  //FIXME
            model_name_(steer<std::string>("model")),
            reg_prop_(steer<double>("regProp")) {
        CG_INFO("MartyProcess") << "List of models registered: " << MartyModelFactory::get().models() << ".";
        // define the incoming-outgoing central system
        setIntermediatePartons({steer<ParticleProperties>("ip1").pdgid, steer<ParticleProperties>("ip2").pdgid});
        setProducedParticles(std::vector<pdgid_t>(2, steer<ParticleProperties>("output").pdgid));
        // start computing the Feynman rules for the given model
        auto model = MartyModelFactory::get().build(model_name_);
        try {
          auto ip1 = mty::GetParticle(*model, toMartyFields(intermediate_parts_.at(0)));
          auto ip2 = mty::GetParticle(*model, toMartyFields(intermediate_parts_.at(1)));
          auto central = mty::GetParticle(*model, toMartyFields(produced_parts_.at(0)));
          auto proc = model->computeAmplitude(
              mty::Order::TreeLevel,
              {mty::Incoming(ip1), mty::Incoming(ip2), mty::Outgoing(central), mty::Outgoing(mty::AntiPart(central))});
          auto sq_ampl = model->computeSquaredAmplitude(proc);
          auto extract_variables = [](auto& expr) -> std::set<csl::Expr> {
            std::set<csl::Expr> vars;
            csl::VisitEachLeaf(expr, [&](const csl::Expr& sub) {
              if (csl::IsConstant(sub))
                vars.insert(sub);
            });
            return vars;
          };
          me_expr_ = csl::Evaluated(csl::DeepExpanded(sq_ampl), csl::eval::all);
          me_params_ = extract_variables(me_expr_);
          CG_LOG << "Extracted model parameters for process: " << me_params_;
          //csl::LibFunction fct("me", me_expr_, nullptr, true);
        } catch (const mty::error::Type& err) {
          throw CG_FATAL("MartyProcess") << "Error encountered: " << err;
        }
      }
      ProcessPtr clone() const override { return ProcessPtr(new MartyProcess(*this)); }

      static std::string toMartyFields(pdgid_t pdgid) {
        switch (pdgid) {
          case PDG::gluon:
            return "G";
          case PDG::photon:
            return "A";
          case 12:
            return "nu_e";
          case 14:
            return "nu_mu";
          case 16:
            return "nu_tau";
          default:
            return PDG::get()(pdgid).name;
        }
      }

      static ParametersDescription description() {
        auto desc = Process2to4::description();
        desc.setDescription("Generic Marty process");
        desc.add<std::string>("model", "sm").setDescription("Marty model considered");
        desc.add<pdgid_t>("ip1", PDG::photon).setDescription("positive-z axis incoming parton");
        desc.add<pdgid_t>("ip2", PDG::photon).setDescription("negative-z axis incoming parton");
        desc.add<pdgid_t>("output", PDG::invalid).setDescription("outgoing central system's PDG identifier");
        desc.add<double>("regProp", 0.).setDescription("divergence suppression parameter");
        return desc;
      }

    private:
      void prepareProcessKinematics() override {}
      double computeCentralMatrixElement() const override {
        //throw CG_FATAL("MartyProcess:computeCentralMatrixElement") << "Not yet implemented!";
        static auto to_value = [&](auto& to_replace) {
          const auto& var_name = to_replace->getName();
          if (var_name == "s_12")
            return (q1() + q2()).mass2();
          if (var_name == "s_13")
            return (q1() + pc(0)).mass2();
          if (var_name == "s_14")
            return (q1() + pc(1)).mass2();
          if (var_name == "s_23")
            return (q2() + pc(0)).mass2();
          if (var_name == "s_24")
            return (q2() + pc(1)).mass2();
          if (var_name == "s_34")
            return (pc(0) + pc(1)).mass2();
          if (var_name == "reg_prop")
            return reg_prop_;
          throw CG_FATAL("MartyProcess") << "Failed to recognise the parameter name '" << var_name << "'.";
        };
        auto me_expr = me_expr_;
        for (const auto& param : me_params_)
          csl::Replace(me_expr, param, csl::float_s(to_value(param)));
        auto eval = csl::GetComplexModulus(csl::Evaluated(me_expr, csl::eval::all));
        if (csl::GetType(eval) != csl::Type::Float && csl::GetType(eval) != csl::Type::Integer)
          throw CG_FATAL("MartyProcess") << "Invalid type retrieved after evaluation: '" << csl::GetType(eval) << "'.";
        return eval->evaluateScalar();
      }

      const std::string model_name_;
      const double reg_prop_;
      csl::Expr me_expr_;
      std::set<csl::Expr> me_params_;
    };
  }  // namespace proc
}  // namespace cepgen
// register process
typedef cepgen::proc::MartyProcess MartyProcess;
REGISTER_PROCESS("marty", MartyProcess);
