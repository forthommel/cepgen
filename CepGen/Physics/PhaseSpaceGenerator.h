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

#ifndef CepGen_Physics_PhaseSpaceGenerator_h
#define CepGen_Physics_PhaseSpaceGenerator_h

#include <array>
#include <vector>

#include "CepGen/Physics/Momentum.h"
#include "CepGen/Physics/PDG.h"

namespace cepgen {
  namespace proc {
    class KTProcess;
  }
  /// A phase space mapping utility
  class PhaseSpaceGenerator {
  public:
    explicit PhaseSpaceGenerator(proc::KTProcess&);

    static const Limits x_limits_;  ///< Standard [0,1] limits for input variables
    /// Initialise the phase space generator
    virtual void initialise() = 0;
    /// Generate the output particles' momenta
    virtual bool generate() = 0;

  protected:
    proc::KTProcess& proc_;
    const double ww_;
  };

  struct PhaseSpaceGenerator2to3 final : public PhaseSpaceGenerator {
    using PhaseSpaceGenerator::PhaseSpaceGenerator;

    void initialise() override;
    bool generate() override;

    double y_c_, pt_c_;
  };

  struct PhaseSpaceGenerator2to4 final : public PhaseSpaceGenerator {
    using PhaseSpaceGenerator::PhaseSpaceGenerator;

    void initialise() override;
    bool generate() override;

    double y_c1_, y_c2_, pt_diff_, phi_pt_diff_;
  };
}  // namespace cepgen

#endif
