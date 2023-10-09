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

#include <cmath>

#include "CepGen/Utils/Math.h"

namespace cepgen {
  namespace utils {
    template <typename T>
    bool positive(const T& val) {
      return val > T{} && std::isfinite(val);
    }
    template bool positive<double>(const double&);
    template bool positive<float>(const float&);
    template bool positive<int>(const int&);

    double delta(double x, double epsilon) {
      if (epsilon < 0.)
        epsilon = 1.e-5;
      if (fabs(x) < epsilon)
        return 2. * M_1_PI / epsilon / epsilon * std::hypot(epsilon, x);
      return 0.;
    }
  }  // namespace utils
}  // namespace cepgen
