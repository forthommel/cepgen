/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2020-2023  Laurent Forthomme
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

#ifndef CepGenAddOns_MadGraphWrapper_MadGraphProcessBuilder_h
#define CepGenAddOns_MadGraphWrapper_MadGraphProcessBuilder_h

#include "CepGen/Core/SteeredObject.h"
#include "CepGenAddOns/MadGraphWrapper/MadGraphProcess.h"

namespace cepgen {
  class MadGraphProcessBuilder : public SteeredObject<MadGraphProcessBuilder> {
  public:
    explicit MadGraphProcessBuilder(const ParametersList&);

    static ParametersDescription description();

  protected:
    const std::string params_card_;
    std::shared_ptr<MadGraphProcess> mg5_proc_;
  };
}  // namespace cepgen

#endif
