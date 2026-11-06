/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2026 Laurent Forthomme
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef CepGenMadGraph_OffShell_h
#define CepGenMadGraph_OffShell_h

// Header-only (no CepGen dependency): it is included both by the process template and, through the hook inserted by
// mg5amc::patchHelAmps, by the model's HelAmps_*.cc. Both are compiled in the same shared library, hence share the
// thread-local state below.

#include <algorithm>
#include <array>
#include <complex>
#include <cstddef>

namespace cepgen::mg5amc::offshell {
  /// State of one overridden incoming leg
  struct Leg {
    bool active{false};
    std::array<double, 4> q{};    ///< off-shell four-momentum (E, px, py, pz), as passed to the matrix element
    std::array<double, 4> eps{};  ///< polarisation four-vector (real) to use in place of the HELAS wavefunction
    std::size_t hits{0};          ///< number of times this override replaced a stock wavefunction
  };
  inline std::array<Leg, 2>& legs() {
    thread_local std::array<Leg, 2> legs_;
    return legs_;
  }
  /// Override the wavefunction of an incoming leg
  inline void set(std::size_t leg, const double* q, const std::array<double, 4>& eps) {
    auto& l = legs().at(leg);
    l.active = true;
    std::copy(q, q + 4, l.q.begin());
    l.eps = eps;
    l.hits = 0;
  }
  /// Restore the stock behaviour for all legs
  inline void reset() {
    for (auto& l : legs())
      l.active = false;
  }
  /// Number of times a leg override was applied since it was last set
  inline std::size_t hits(std::size_t leg) { return legs().at(leg).hits; }
  /// Fill a HELAS vector wavefunction if the incoming momentum matches an overridden leg
  /// \param[in] p four-momentum passed to HELAS vxxxxx
  /// \param[in] nsv incoming (-1) or outgoing (+1) flag
  /// \param[out] vc HELAS vector wavefunction
  /// \return true if the wavefunction was filled here, false if the stock routine must run
  inline bool fill(const double* p, int nsv, std::complex<double>* vc) {
    if (nsv != -1)  // incoming legs only
      return false;
    for (auto& l : legs()) {
      if (!l.active || !std::equal(p, p + 4, l.q.begin()))
        continue;
      vc[0] = std::complex<double>(p[0] * nsv, p[3] * nsv);
      vc[1] = std::complex<double>(p[1] * nsv, p[2] * nsv);
      for (std::size_t mu = 0; mu < 4; ++mu)
        vc[2 + mu] = l.eps[mu];
      ++l.hits;
      return true;
    }
    return false;
  }
}  // namespace cepgen::mg5amc::offshell

#endif
