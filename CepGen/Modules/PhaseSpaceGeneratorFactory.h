/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2024  Laurent Forthomme
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
/** \file */

#ifndef CepGen_Modules_PhaseSpaceGeneratorFactory_h
#define CepGen_Modules_PhaseSpaceGeneratorFactory_h

#include "CepGen/Modules/ModuleFactory.h"

/// Add a central phase space generator to the list of handled modules
#define REGISTER_PSGEN(name, obj)                                                       \
  namespace cepgen {                                                                    \
    struct BUILDERNM(obj) {                                                             \
      BUILDERNM(obj)() { PhaseSpaceGeneratorFactory::get().registerModule<obj>(name); } \
    };                                                                                  \
    static const BUILDERNM(obj) gPhaseSpaceGen##obj;                                    \
  }                                                                                     \
  static_assert(true, "")

namespace cepgen {
  class PhaseSpaceGenerator;
  /// A phase space mapping algorithms factory
  DEFINE_FACTORY(std::string, PhaseSpaceGeneratorFactory, PhaseSpaceGenerator, "Phase space generator factory");
}  // namespace cepgen

#endif
