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

#include "CepGen/Core/Exception.h"
#include "CepGen/Generator.h"
#include "CepGen/Utils/AbortHandler.h"
#include "CepGenAddOns/MadGraphWrapper/MadGraphInterface.h"
#include "CepGenAddOns/MadGraphWrapper/MadGraphProcessBuilder.h"

namespace cepgen {
  MadGraphProcessBuilder::MadGraphProcessBuilder(const ParametersList& params)
      : SteeredObject(params), params_card_(steer<std::string>("parametersCard")) {
    utils::AbortHandler();
    try {
      const auto& lib_file = steer<std::string>("lib");
      if (!lib_file.empty())
        loadLibrary(lib_file);
      else {
        const MadGraphInterface interf(params);
        loadLibrary(interf.run());
      }
    } catch (const utils::RunAbortedException&) {
      CG_FATAL("MadGraphProcessBuilder") << "MadGraph_aMC process generation aborted.";
    }
    // once MadGraph process library is loaded into runtime environment, can define its wrapper object
    mg5_proc_.reset(new MadGraphProcess);
  }

  ParametersDescription MadGraphProcessBuilder::description() {
    auto desc = ParametersDescription();
    desc.setDescription("MadGraph_aMC process builder");
    desc.add<std::string>("lib", "").setDescription("Precompiled library for this process definition");
    desc.add<std::string>("parametersCard", "param_card.dat").setDescription("Runtime MadGraph parameters card");
    desc += MadGraphInterface::description();
    return desc;
  }
}  // namespace cepgen
