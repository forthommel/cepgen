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

#include <HepLib.h>

#include "CepGen/Event/Event.h"
#include "CepGen/Modules/ProcessFactory.h"
#include "CepGen/Process/Process.h"

using namespace GiNaC;
using namespace HepLib;

namespace cepgen {
  class EEtoLL : public proc::Process {
  public:
    explicit EEtoLL(const ParametersList& params) : proc::Process(params) {
      proc_.Model = R"EOF(
[e, ebar, -]
[mu, mubar, -]
[A, A, +]
[ebar, e, A]
[mubar, mu, A]
)EOF";
      proc_.In = "e[p], ebar[P]";
      proc_.Out = " mubar[k], mu[K]";
      proc_.Options = "onshell";
      proc_.Loops = 0;
      Vector p("p"), P("P"), k("k"), K("K");
      symtab st{{"p", p}, {"P", P}, {"k", k}, {"K", K}};
      auto amps = proc_.Amplitudes(st);

      /*auto amps_feyn_rule = MapFunction([&](const ex& e, MapFunction& self) -> ex {
        if (isFunction(e, "OutField") || isFunction(e, "InField"))
          return 1;
        else if (isFunction(e, "Propagator")) {
          auto fi1 = e.op(0).op(1);
          auto fi2 = e.op(1).op(1);
          auto mom = e.op(2);
          if (e.op(0).op(0) == A)
            return (-I) * SP(LI(fi1), LI(fi2)) / SP(mom);  // Feynman Gauge
          else if (e.op(0).op(0) == ebar)
            return I * Matrix(GAS(mom) + GAS(1) * me, DI(fi1), DI(fi2)) / (SP(mom) - me * me);
          else if (e.op(0).op(0) == mubar)
            return I * Matrix(GAS(mom) + GAS(1) * mm, DI(fi1), DI(fi2)) / (SP(mom) - mm * mm);
        } else if (isFunction(e, "Vertex")) {
          auto fi1 = e.op(0).op(1);
          auto fi2 = e.op(1).op(1);
          auto fi3 = e.op(2).op(1);
          if (e.op(0).op(0) == ebar)
            return I * Symbol("e") * Matrix(GAS(LI(fi3)), DI(fi1), DI(fi2));
          else if (e.op(0).op(0) == mubar)
            return I * Symbol("e") * Matrix(GAS(LI(fi3)), DI(fi1), DI(fi2));
        }
        return e.map(self);
      })(amps);*/
    }

    proc::ProcessPtr clone() const override { return proc::ProcessPtr(new EEtoLL(*this)); }

    static ParametersDescription description() {
      auto desc = proc::Process::description();
      desc.setDescription("e+e- -> Z -> l+l- process");
      return desc;
    }

    void addEventContent() override {}

    void prepareKinematics() override {
      //defineVariable();
    }
    void fillKinematics(bool) override {}
    double computeWeight() override {
      Vector p("p"), P("P"), k("k"), K("K");
      symtab st{{"p", p}, {"P", P}, {"k", k}, {"K", K}};
      auto amps = proc_.Amplitudes(st);
      return 0.;
    }

  private:
    QGRAF::Process proc_;
  };
}  // namespace cepgen

REGISTER_PROCESS("eetoll", EEtoLL);
