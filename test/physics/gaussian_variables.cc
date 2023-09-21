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

#include "CepGen/Core/ParametersList.h"
#include "CepGen/Generator.h"
#include "CepGen/Parameters.h"
#include "CepGen/Process/Process.h"
#include "CepGen/Utils/ArgumentsParser.h"
#include "CepGen/Utils/Histogram.h"
#include "CepGen/Utils/Test.h"

using namespace std;

struct GaussianVarProcess : cepgen::proc::Process {
  explicit GaussianVarProcess()
      : cepgen::proc::Process(
            cepgen::proc::Process::description().validate(cepgen::ParametersList()).set<bool>("hasEvent", false)) {}
  cepgen::proc::ProcessPtr clone() const override { return cepgen::proc::ProcessPtr(new GaussianVarProcess); }
  double computeWeight() override { return 1.; }
  void prepareKinematics() override {
    defineVariable(gaussian_var, cepgen::proc::Process::Mapping::gaussian, {0., 1.});
  }

  double gaussian_var{0.};
};

int main(int argc, char* argv[]) {
  int num_events;
  string worker;
  cepgen::ArgumentsParser(argc, argv)
      .addOptionalArgument("num-events,n", "number of events to generate", &num_events, 100)
      .addOptionalArgument("worker,w", "generator worker algorithm", &worker, "grid_optimised")
      .parse();

  cepgen::Generator gen;
  gen.parametersRef().setProcess(new GaussianVarProcess);
  gen.parametersRef().generation().setParameters(
      cepgen::ParametersList().set("worker", cepgen::ParametersList().setName(worker)));
  cepgen::utils::Hist1D hist(100, {-1., 1.});

  gen.computeXsection();
  gen.generate(num_events, [&hist](const cepgen::proc::Process& proc) {
    hist.fill(dynamic_cast<const GaussianVarProcess&>(proc).gaussian_var);
    CG_LOG << dynamic_cast<const GaussianVarProcess&>(proc).gaussian_var;
  });

  CG_TEST_SUMMARY;
}
