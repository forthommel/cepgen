/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2022-2023  Laurent Forthomme
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

#include <fstream>

#include "CepGen/CollinearFluxes/IntegratedPartonFlux.h"
#include "CepGen/Core/Exception.h"
#include "CepGen/Generator.h"
#include "CepGen/Integration/AnalyticIntegrator.h"
#include "CepGen/Modules/AnalyticIntegratorFactory.h"
#include "CepGen/Modules/DrawerFactory.h"
#include "CepGen/Modules/PartonFluxFactory.h"
#include "CepGen/Physics/PDG.h"
#include "CepGen/Utils/ArgumentsParser.h"
#include "CepGen/Utils/Drawer.h"
#include "CepGen/Utils/Graph.h"
#include "CepGen/Utils/String.h"

using namespace std;

int main(int argc, char* argv[]) {
  vector<string> fluxes;
  int num_points;
  string integrator, output_file, plotter, beams;
  bool logx, logy, draw_grid;
  vector<double> rescl, sqrts, accept;
  cepgen::Limits mx_range, y_range;
  vector<cepgen::Limits> xi_ranges(1);

  cepgen::initialise();

  cepgen::ArgumentsParser(argc, argv)
      .addOptionalArgument(
          "flux,f", "(collinear) flux modelling(s)", &fluxes, cepgen::IntegratedPartonFluxFactory::get().modules())
      .addOptionalArgument("beams,b", "beams info (PDGid1:pz1[,PDGid2:pz2])", &beams, "")
      .addOptionalArgument("rescaling,r", "luminosity rescaling", &rescl, vector<double>{1.})
      .addOptionalArgument("integrator,i", "type of integration algorithm", &integrator, "gsl")
      .addOptionalArgument("sqrts,s", "two-proton centre of mass energy (GeV)", &sqrts, vector<double>{13.e3})
      .addOptionalArgument("mrange,m", "two-photon mass range", &mx_range, cepgen::Limits{0., 100.})
      .addOptionalArgument("npoints,n", "number of x-points to scan", &num_points, 50)
      .addOptionalArgument("output,o", "output file name", &output_file, "collflux.int.scan.output.txt")
      .addOptionalArgument("plotter,p", "type of plotter to user", &plotter, "")
      .addOptionalArgument("logx", "logarithmic x-scale", &logx, false)
      .addOptionalArgument("logy,l", "logarithmic y-scale", &logy, false)
      .addOptionalArgument("xi-range,x", "acceptance range for proton momentum loss", &xi_ranges[0])
      .addOptionalArgument("accept", "pairs of min/max acceptance ranges for proton momentum loss", &accept)
      .addOptionalArgument("yrange,y", "y plot range", &y_range)
      .addOptionalArgument("draw-grid,g", "draw the x/y grid", &draw_grid, false)
      .parse();

  ofstream out(output_file);
  if (logx && mx_range.min() == 0.)
    mx_range.min() = 1.e-3;

  //----- two possible strategies handled for the centre of mass energy computation

  if (beams.empty()) {  // centre of mass energy is directly specified _for each flux_.
    if (sqrts.size() != fluxes.size())
      sqrts = vector<double>(fluxes.size(), sqrts.at(0));
  } else {  // one single centre of mass energy for all fluxes ; unpack the beams info, and compute the centre of mass energy
    vector<int> pdgids;
    vector<double> pzs;
    size_t i = 0;
    for (const auto& beam_info : cepgen::utils::split(beams, ',')) {
      ++i;
      const auto info = cepgen::utils::split(beam_info, ':');
      if (info.size() != 2)
        throw CG_FATAL("main") << "Invalid format for beam " << i << ": should be 'id:pz(GeV)'.";
      pdgids.emplace_back(std::stoi(info.at(0)));
      pzs.emplace_back(std::stod(info.at(1)));
    }
    cepgen::Momentum mom1, mom2;
    if (i == 1) {
      mom1 = cepgen::Momentum::fromPxPyPzM(0., 0., pzs.at(0), cepgen::PDG::get().mass(pdgids.at(0)));
      mom2 = cepgen::Momentum::fromPxPyPzM(0., 0., -pzs.at(0), cepgen::PDG::get().mass(pdgids.at(0)));
    } else {
      mom1 = cepgen::Momentum::fromPxPyPzM(0., 0., pzs.at(0), cepgen::PDG::get().mass(pdgids.at(0)));
      mom2 = cepgen::Momentum::fromPxPyPzM(0., 0., -pzs.at(1), cepgen::PDG::get().mass(pdgids.at(1)));
    }
    sqrts = vector<double>(fluxes.size(), (mom1 + mom2).mass());
  }
  CG_LOG << sqrts;
  if (rescl.size() != fluxes.size())
    rescl = vector<double>(fluxes.size(), rescl.at(0));
  if (!accept.empty()) {
    xi_ranges.clear();
    if (accept.size() % 2 != 0)
      throw CG_FATAL("main")
          << "Invalid acceptance(s) list specified! Supported format is an array of (min1, max2, min2, max2, ...)";
    for (size_t i = 0; i < accept.size() / 2; ++i) {
      cepgen::Limits rng{accept.at(2 * i), accept.at(2 * i + 1)};
      if (rng == cepgen::Limits{0., 1.})
        rng = cepgen::Limits();
      xi_ranges.emplace_back(rng);
    }
    if (!xi_ranges.empty())
      CG_LOG << "x (xi) acceptance cuts defined: " << xi_ranges << ".";
  }

  //----- prepare the output file, start the computation of points

  out << "# fluxes: " << cepgen::utils::merge(fluxes, ",") << "\n"
      << "# two-photon mass range: " << mx_range;
  map<string, vector<cepgen::utils::Graph1D> > m_gr_fluxes;  // {collinear flux -> graph}
  const auto mxvals = mx_range.generate(num_points, logx);
  vector<vector<double> > values(num_points);
  auto integr = cepgen::AnalyticIntegratorFactory::get().build(
      cepgen::ParametersList().setName<string>(integrator).set<int>("mode", 0).set<int>("nodes", 2000));
  for (size_t i = 0; i < fluxes.size(); ++i) {
    const auto flux_names = cepgen::utils::split(fluxes.at(i), '+');
    auto flux1 = cepgen::IntegratedPartonFluxFactory::get().build(flux_names.at(0));
    auto flux2 =
        cepgen::IntegratedPartonFluxFactory::get().build(flux_names.size() > 1 ? flux_names.at(1) : flux_names.at(0));
    ostringstream foss;
    foss << flux1->name() << "/" << flux2->name();
    const auto s = sqrts.at(i) * sqrts.at(i);
    for (const auto& xi_range : xi_ranges) {
      ostringstream oss(foss.str());
      if (xi_range.valid())
        oss << " (" << xi_range.min() << " < \\xi < " << xi_range.max() << ")";
      m_gr_fluxes[fluxes.at(i)].emplace_back("", oss.str());
    }

    for (int j = 0; j < num_points; ++j) {
      const auto& mx = mxvals.at(j);
      for (size_t k = 0; k < xi_ranges.size(); ++k) {
        const auto& xi_range = xi_ranges.at(k);
        auto lumi_wgg = 2. * mx / s *
                        integr->integrate(
                            [&xi_range, &mx, &s, &flux1, &flux2](double x) {
                              if (xi_range.valid() && (!xi_range.contains(x) || !xi_range.contains(mx * mx / x / s)))
                                return 0.;
                              return flux1->flux(x) * flux2->flux(mx * mx / x / s) / x;
                            },
                            cepgen::Limits(mx * mx / s, 1.));
        lumi_wgg *= rescl.at(i);
        values.at(j).emplace_back(lumi_wgg);
        m_gr_fluxes[fluxes.at(i)][k].addPoint(mx, lumi_wgg);
      }
    }
  }
  for (int i = 0; i < num_points; ++i)
    out << "\n" << mxvals.at(i) << "\t" << cepgen::utils::merge(values.at(i), "\t");
  out.close();

  //----- plot things if necessary

  if (!plotter.empty()) {
    auto plt = cepgen::DrawerFactory::get().build(plotter);
    cepgen::utils::Drawer::Mode dm;
    if (logx)
      dm |= cepgen::utils::Drawer::Mode::logx;
    if (logy)
      dm |= cepgen::utils::Drawer::Mode::logy;
    if (draw_grid)
      dm |= cepgen::utils::Drawer::Mode::grid;
    cepgen::utils::DrawableColl coll;
    size_t i = 0;
    for (auto& cf_gr : m_gr_fluxes) {
      for (auto& gr_xi : cf_gr.second) {
        string units;
        if (rescl.at(i) != 1.)
          units = "cm$^{-2}$ s$^{-1}$ ";
        gr_xi.xAxis().setLabel("$\\omega_{\\gamma\\gamma}$ (GeV)");
        gr_xi.yAxis().setLabel("d$L_{\\gamma\\gamma}$/d$\\omega_{\\gamma\\gamma}$ (" + units + "GeV$^{-1}$)");
        if (y_range.valid())
          gr_xi.yAxis().setRange(y_range);
        coll.emplace_back(&gr_xi);
      }
      ++i;
    }
    plt->draw(coll, "comp_photonlumi", "", dm);
  }
  return 0;
}
