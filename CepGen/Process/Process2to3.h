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

#ifndef CepGen_Process_Process2to3_h
#define CepGen_Process_Process2to3_h

#include "CepGen/Process/KTProcess.h"

namespace cepgen {
  namespace proc {
    /// A 2-to-3 (or 2-to-1 central) process
    class Process2to3 : public KTProcess {
    public:
      /// Initialise a 2-to-3 process
      /// \param[in] params Collection of user-defined steering parameters
      /// \param[in] partons Incoming hard scattering particles
      /// \param[in] cs_id Central particles PDG id
      explicit Process2to3(const ParametersList& params, std::array<pdgid_t, 2> partons, pdgid_t cs_id);

      static ParametersDescription description();

    protected:
      /// Set all cuts for the single outgoing particle phase space definition
      void setCuts(const cuts::Central& single);

      void preparePhaseSpace() override;
      void fillCentralParticlesKinematics() override;
      double computeKTFactorisedMatrixElement() override;

      /// Conform all kinematics variables to the user-defined phase space
      virtual void prepareProcessKinematics() = 0;
      /// Computation rule for the central matrix element
      virtual double computeCentralMatrixElement() const = 0;

      cuts::Central single_limits_;  ///< Limits to be applied on single central system's particle
    };
  }  // namespace proc
}  // namespace cepgen

#endif
