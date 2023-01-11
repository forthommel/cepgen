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

#ifndef CepGen_Integration_GridParameters_h
#define CepGen_Integration_GridParameters_h

#include <cstddef>
#include <unordered_map>
#include <vector>

namespace cepgen {
  class Integrator;
  /// A parameters placeholder for the grid integration helper
  class GridParameters {
  public:
    /// Build a generation grid for a ndim-dimensional phase space
    explicit GridParameters(size_t ndim);

    /// Integration grid size parameter
    static constexpr unsigned short M_BIN = 3;
    /// Weight of each grid coordinate
    static constexpr double INV_M_BIN = 1. / M_BIN;
    /// Coordinates definition
    typedef std::vector<unsigned short> coord_t;

    /// Dump the grid coordinates
    void dump() const;

    /// Grid multiplicity
    size_t size() const;
    /// Number of times a phase space point has been randomly selected
    const coord_t& n(size_t) const;
    /// Global function maximum
    float globalMax() const { return f_max_global_; }
    /// Maximal function value for a given grid coordinate
    float maxValue(size_t) const;
    /// Set the function value for a given grid coordinate
    void setValue(size_t, float);
    /// Shoot a phase space point for a grid coordinate
    void shoot(const Integrator*, size_t, std::vector<double>& out) const;
    /// Specify a new trial has been attempted for bin
    void increment(size_t coord);
    /// Number of points already shot for a given grid coordinate
    size_t numPoints(size_t coord) const;
    /// Has the grid been prepared
    bool prepared() const { return gen_prepared_; }
    /// Mark the grid as prepared
    void setPrepared(bool prepared = true) { gen_prepared_ = prepared; }
    /// Correction to apply on the next phase space point generation
    float correctionValue(size_t bin) const { return bin_correc_.at(bin).correc; }
    /// Set the correction to apply on the next phase space point generation
    void setCorrectionValue(size_t bin, float correc) { bin_correc_[bin].correc = correc; }
    /// Apply the correction requested at the previous generation
    bool correct(size_t);
    void rescale(size_t, float);
    void initCorrectionCycle(size_t, float);
    double maxHist(const Integrator*, size_t bin) const;

  private:
    void generateCoordinates(coord_t&, size_t) const;
    /// Phase space multiplicity
    size_t ndim_{0};
    /// Has the grid been already prepared?
    bool gen_prepared_{false};
    /// Point coordinates in grid
    std::vector<coord_t> coords_;
    /// Maximal value of the function in the considered integration range
    float f_max_global_{0.};
    /// Binned correction
    struct BinCorrection {
      size_t num_points{0ul};  ///< Number of functions values evaluated for this point
      float correc{0.};        ///< Correction to apply on the next phase space point generation
      float correc2{0.};
      float f_max{0.};  ///< Maximal value of the function at one given point
      float f_max2{0.};
      float f_max_diff{0.};
      float f_max_old{0.};
    };
    /// Correction to apply on the next phase space point generation
    std::vector<BinCorrection> bin_correc_;
  };
}  // namespace cepgen

#endif
