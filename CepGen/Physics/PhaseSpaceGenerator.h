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

#include <vector>

#include "CepGen/Physics/Momentum.h"
#include "CepGen/Physics/PDG.h"

namespace cepgen {
  namespace proc {
    class Process;
  }
  /// A phase space mapping utility
  class PhaseSpaceGenerator {
  public:
    explicit PhaseSpaceGenerator(proc::Process&, const std::vector<pdgid_t>&);

    /// Initialise the phase space generator
    void initialise();
    /// Generate the output particles' momenta
    const std::vector<Momentum>& generate();

  private:
    proc::Process& proc_;
    const std::vector<pdgid_t> output_;
    std::vector<Momentum> momenta_;
  };
}  // namespace cepgen

#endif
