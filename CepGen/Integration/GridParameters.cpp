/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2018-2023  Laurent Forthomme
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

#include <cmath>  // pow

#include "CepGen/Core/Exception.h"
#include "CepGen/Integration/GridParameters.h"
#include "CepGen/Integration/Integrator.h"
#include "CepGen/Utils/String.h"

namespace cepgen {
  GridParameters::GridParameters(size_t ndim) : ndim_(ndim) {
    //--- build and populate the grid
    coord_t coord(ndim, 0);
    for (size_t i = 0; i < (size_t)pow(M_BIN, ndim_); ++i) {
      generateCoordinates(coord, i);
      coords_.emplace_back(coord);
      bin_correc_.emplace_back();
    }
  }

  size_t GridParameters::size() const { return coords_.size(); }

  const GridParameters::coord_t& GridParameters::n(size_t coord) const { return coords_.at(coord); }

  void GridParameters::setValue(size_t coord, float val) {
    //--- update function local and global maxima if needed
    bin_correc_.at(coord).f_max = std::max(bin_correc_.at(coord).f_max, val);
    f_max_global_ = std::max(f_max_global_, val);
  }

  float GridParameters::maxValue(size_t coord) const { return bin_correc_.at(coord).f_max; }

  size_t GridParameters::numPoints(size_t coord) const { return bin_correc_.at(coord).num_points; }

  void GridParameters::increment(size_t coord) { bin_correc_.at(coord).num_points++; }

  void GridParameters::shoot(const Integrator* integr, size_t coord, std::vector<double>& out) const {
    const auto& nv = coords_.at(coord);
    for (size_t i = 0; i < nv.size(); ++i)
      out[i] = (integr->uniform() + nv.at(i)) * INV_M_BIN;
  }

  void GridParameters::dump() const {
    CG_INFO("GridParameters:dump").log([&](auto& info) {
      for (size_t i = 0; i < coords_.size(); ++i) {
        const auto& correc = bin_correc_.at(i);
        info << "\nn[" << i << "]: "
             << "coord=" << coords_.at(i) << ", "
             << "num points: " << correc.num_points << ", "
             << "max=" << correc.f_max << ".";
      }
    });
  }

  void GridParameters::generateCoordinates(coord_t& coord, size_t i) const {
    size_t jj = i;
    for (size_t j = 0; j < ndim_; ++j) {
      size_t tmp = jj * INV_M_BIN;
      coord[j] = jj - tmp * M_BIN;
      jj = tmp;
    }
  }

  bool GridParameters::correct(size_t bin) {
    if (bin >= bin_correc_.size())
      throw CG_FATAL("GridParameters:correct") << "No correction parameter found for bin " << bin << ".";
    auto& correc = bin_correc_.at(bin);
    if (correc.f_max2 <= correc.f_max)
      return true;
    correc.f_max_old = correc.f_max;
    correc.f_max_diff = correc.f_max2 - correc.f_max_old;
    correc.correc = (correc.num_points - 1) * correc.f_max_diff / f_max_global_;
    if (correc.f_max2 >= f_max_global_)
      correc.correc *= correc.f_max2 / f_max_global_;
    setValue(bin, correc.f_max2);
    correc.correc -= correc.correc2;
    correc.correc2 = 0.;
    correc.f_max2 = 0.;
    return false;
  }

  void GridParameters::rescale(size_t bin, float weight) {
    if (bin >= bin_correc_.size())
      throw CG_FATAL("GridParameters:correct") << "No correction parameter found for bin " << bin << ".";
    auto& correc = bin_correc_.at(bin);
    if (weight <= correc.f_max)
      return;
    correc.correc += 1.;
    correc.correc2 -= 1.;
    correc.f_max2 = std::max(correc.f_max2, weight);
  }

  void GridParameters::initCorrectionCycle(size_t bin, float weight) {
    if (bin >= bin_correc_.size())
      throw CG_FATAL("GridParameters:correct") << "No correction parameter found for bin " << bin << ".";
    setValue(bin, weight);
    auto& correc = bin_correc_.at(bin);
    correc.f_max_old = correc.f_max;
    correc.f_max_diff = weight - correc.f_max_old;
    correc.correc = (correc.num_points - 1) * correc.f_max_diff / f_max_global_ - 1.;

    CG_DEBUG("GridParameters:initCorrectionCycle")
        << "Correction " << bin_correc_.at(bin).correc << " will be applied "
        << "for phase space bin " << bin << " (" << utils::s("point", correc.num_points, true) << "). "
        << "Maxima ratio: " << (correc.f_max_diff / f_max_global_) << ".";
  }

  double GridParameters::maxHist(const Integrator* integr, size_t bin) const {
    if (bin >= bin_correc_.size())
      throw CG_FATAL("GridParameters:correct") << "No correction parameter found for bin " << bin << ".";
    const auto& correc = bin_correc_.at(bin);
    return correc.f_max_old + integr->uniform(0., correc.f_max_diff);
  }
}  // namespace cepgen
