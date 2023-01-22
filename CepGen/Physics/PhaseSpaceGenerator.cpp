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

#include "CepGen/Physics/PhaseSpaceGenerator.h"
#include "CepGen/Process/Process.h"

namespace cepgen {

  template <>
  void PhaseSpaceGenerator<2>::initialise() {
    proc_.defineVariable(vars_.at(0),  // y_c1_,
                         proc::Process::Mapping::linear,
                         proc_.kinematics().cuts().central.rapidity_single,
                         {-6., 6.},
                         "First outgoing particle rapidity");
    proc_.defineVariable(vars_.at(1),  // y_c2_,
                         proc::Process::Mapping::linear,
                         proc_.kinematics().cuts().central.rapidity_single,
                         {-6., 6.},
                         "Second outgoing particle rapidity");
    proc_.defineVariable(vars_.at(3),  // pt_diff_,
                         proc::Process::Mapping::linear,
                         proc_.kinematics().cuts().central.pt_diff,
                         {0., 500.},
                         "Final state particles transverse momentum difference");
    proc_.defineVariable(vars_.at(4),  // phi_pt_diff_,
                         proc::Process::Mapping::linear,
                         proc_.kinematics().cuts().central.phi_diff,
                         {0., 2. * M_PI},
                         "Final state particles azimuthal angle difference");
  }

  template <>
  const PhaseSpaceGenerator<2>::Momenta& PhaseSpaceGenerator<2>::generate() {
    return momenta_;
  }
}  // namespace cepgen
