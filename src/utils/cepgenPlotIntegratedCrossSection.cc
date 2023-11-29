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

#include "CepGen/Cards/Handler.h"
#include "CepGen/Core/Exception.h"
#include "CepGen/Generator.h"
#include "CepGen/Modules/DrawerFactory.h"
#include "CepGen/Parameters.h"
#include "CepGen/Process/Process.h"
#include "CepGen/Utils/ArgumentsParser.h"
#include "CepGen/Utils/Drawer.h"
#include "CepGen/Utils/Graph.h"
#include "CepGen/Utils/ProcessFunctional.h"
#include "CepGen/Utils/String.h"

using namespace std;

int main(int argc, char* argv[]) {
  string input_card, plotter;
  int npoints;
  bool draw_grid, log;
  cepgen::Limits s_range;

  cepgen::ArgumentsParser(argc, argv)
      .addArgument("input,i", "input card", &input_card)
      .addOptionalArgument("s-range,s", "range of s", &s_range, cepgen::Limits{1., 10.})
      .addOptionalArgument("num-points,n", "number of points to probe", &npoints, 100)
      .addOptionalArgument("draw-grid,g", "draw the x/y grid", &draw_grid, false)
      .addOptionalArgument("log,l", "logarithmic axis", &log, false)
      .addOptionalArgument("plotter,p", "type of plotter to user", &plotter, "")
      .parse();

  cepgen::Generator gen;
  gen.setParameters(cepgen::card::Handler::parseFile(input_card));

  auto pf = cepgen::proc::ProcessFunctional(
      cepgen::ParametersList{}, dynamic_cast<cepgen::proc::FactorisedProcess&>(gen.parametersRef().process()));

  const auto sqrts = gen.parametersRef().process().kinematics().incomingBeams().sqrtS();

  cepgen::utils::Graph1D gr_s_scan("test_scan");
  for (const auto& s : s_range.generate(npoints)) {
    const auto x = std::sqrt(s) / sqrts;
    gr_s_scan.addPoint(s, pf(x, x));
  }

  if (!plotter.empty()) {
    auto plt = cepgen::DrawerFactory::get().build(plotter);
    cepgen::utils::Drawer::Mode dm;
    if (draw_grid)
      dm |= cepgen::utils::Drawer::Mode::grid;
    if (log)
      dm |= cepgen::utils::Drawer::Mode::logy;
    gr_s_scan.xAxis().setLabel("s (GeV^{2})");
    gr_s_scan.yAxis().setLabel("d$\\sigma$/d$\\hat{s}$ (pb/GeV^{2})");
    plt->draw(gr_s_scan, dm);
  }

  return 0;
}
