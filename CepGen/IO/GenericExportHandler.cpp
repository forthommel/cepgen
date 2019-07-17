#include "CepGen/IO/GenericExportHandler.h"

#include "CepGen/Processes/GenericProcess.h"
#include "CepGen/StructureFunctions/StructureFunctions.h"

#include "CepGen/Physics/Constants.h"
#include "CepGen/Physics/FormFactors.h"

#include "CepGen/Parameters.h"
#include "CepGen/Version.h"

#include <sstream>

namespace cepgen
{
  namespace io
  {
    GenericExportHandler::GenericExportHandler( const std::string& name ) :
      name_( name ), event_num_( 0. )
    {}

    GenericExportHandler::~GenericExportHandler()
    {}

    std::string
    GenericExportHandler::banner( const Parameters& params, const std::string& prep )
    {
      std::ostringstream os;
      os
        << prep << "  ***** Sample generated with CepGen v" << version() << " *****\n"
        << prep << "  * process: " << params.processName() << " (" << params.process()->mode() << ")\n"
        << prep << "  * beams:\n"
        << prep << "  *   - " << params.kinematics.incoming_beams.first << "\n"
        << prep << "  *   - " << params.kinematics.incoming_beams.second << "\n";
      if ( params.process()->mode() != KinematicsMode::ElasticElastic && !params.hadroniserName().empty() )
        os << prep << "  * hadroniser: " << params.hadroniserName() << "\n";
      os
        << prep << "  *--- incoming state\n";
      if ( params.kinematics.cuts.initial.q2.valid() )
        os
          << prep << "  * Q2 range (GeV2): "
          << params.kinematics.cuts.initial.q2 << "\n";
      if ( params.process()->mode() != KinematicsMode::ElasticElastic
        && params.kinematics.cuts.remnants.mass_single.valid() )
        os
          << prep << "  * remnants mass range (GeV/c2): "
          << params.kinematics.cuts.remnants.mass_single << "\n";
      os << prep << "  *--- central system\n";
      if ( params.kinematics.cuts.central.pt_single.valid() )
        os
          << prep << "  * single particle pt (GeV/c): "
          << params.kinematics.cuts.central.pt_single << "\n";
      if ( params.kinematics.cuts.central.energy_single.valid() )
        os
          << prep << "  * single particle energy (GeV): "
          << params.kinematics.cuts.central.energy_single << "\n";
      if ( params.kinematics.cuts.central.eta_single.valid() )
        os
          << prep << "  * single particle eta: "
          << params.kinematics.cuts.central.eta_single << "\n";
      if ( params.kinematics.cuts.central.pt_sum.valid() )
        os
          << prep << "  * total pt (GeV/c): "
          << params.kinematics.cuts.central.mass_sum << "\n";
      if ( params.kinematics.cuts.central.mass_sum.valid() )
        os
          << prep << "  * total invariant mass (GeV/c2): "
          << params.kinematics.cuts.central.mass_sum << "\n";
      os
        << prep << "  **************************************************";
      return os.str();
    }
  }
}

