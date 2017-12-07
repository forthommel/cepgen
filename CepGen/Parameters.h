#ifndef CepGen_Parameters_h
#define CepGen_Parameters_h

#include "CepGen/Core/TamingFunction.h"
#include "CepGen/Core/Integrator.h"
#include "CepGen/Physics/Kinematics.h"
#include "CepGen/Processes/GenericProcess.h"
#include "CepGen/Hadronisers/GenericHadroniser.h"

#include <memory>

namespace CepGen
{
  class Event;
  /// List of parameters used to start and run the simulation job
  class Parameters
  {
    public:
      Parameters();
      /// Copy constructor (transfers ownership to the process/hadroniser!)
      Parameters( Parameters& );
      /// Const copy constructor (all but the process and the hadroniser)
      Parameters( const Parameters& );
      ~Parameters();
      /// Set the polar angle range for the produced leptons
      /// \param[in] thetamin The minimal value of \f$\theta\f$ for the outgoing leptons
      /// \param[in] thetamax The maximal value of \f$\theta\f$ for the outgoing leptons
      void setThetaRange( float thetamin, float thetamax );
      /// Dump the input parameters in the console
      void dump( std::ostream& os = Logger::get().outputStream, bool pretty = true ) const;

      //----- process to compute

      /// Process for which the cross-section will be computed and the events will be generated
      Process::GenericProcess* process() { return process_.get(); }
      /// Name of the process considered
      std::string processName() const;
      /// Set the process to study
      void setProcess( Process::GenericProcess* proc ) { process_.reset( proc ); }

      //----- events kinematics

      /// Events kinematics for phase space definition
      Kinematics kinematics;

      //----- VEGAS

      /// Collection of integrator parameters
      struct IntegratorParameters
      {
        IntegratorParameters() :
          type( Integrator::Vegas ), ncvg( 100000 ), itvg( 10 ),
          npoints( 100 ), first_run( true ), seed( 0 ) {}
        Integrator::Type type;
        /// Number of function calls to be computed for each point
        unsigned int ncvg; // ??
        /// Number of iterations for the integration
        unsigned int itvg;
        /// Number of points to "shoot" in each integration bin by the algorithm
        unsigned int npoints;
        /// Is it the first time the integrator is run?
        bool first_run;
        /// Random number generator seed
        unsigned long seed;
      };
      /// Integrator parameters
      IntegratorParameters integrator;

      //----- events generation

      /// Collection of events generation parameters
      struct Generation
      {
        Generation() :
          enabled( false ), maxgen( 0 ),
          symmetrise( false ), ngen( 0 ), gen_print_every( 1 ) {}
        /// Are we generating events ? (true) or are we only computing the cross-section ? (false)
        bool enabled;
        /// Maximal number of events to generate in this run
        unsigned int maxgen;
        /// Pointer to the last event produced in this run
        std::shared_ptr<Event> last_event;
        /// Do we want the events to be symmetrised with respect to the \f$z\f$-axis ?
        bool symmetrise;
        /// Number of events already generated in this run
        unsigned int ngen;
        /// Frequency at which the events are displayed to the end-user
        unsigned int gen_print_every;
      };
      /// Events generation parameters
      Generation generation;

      /// Specify if the generated events are to be stored
      void setStorage( bool store ) { store_ = store; }
      /// Are the events generated in this run to be stored in the output file ?
      bool storage() const { return store_; }

      //----- hadronisation algorithm

      /// Hadronisation algorithm to use for the proton(s) fragmentation
      Hadroniser::GenericHadroniser* hadroniser() { return hadroniser_.get(); }
      /// Set the hadronisation algorithm
      void setHadroniser( Hadroniser::GenericHadroniser* hadr ) { hadroniser_.reset( hadr ); }
      /// Maximal number of trials for the hadronisation of the proton(s) remnants
      unsigned int hadroniser_max_trials;

      //----- taming functions

      /// Functionals to be used to account for rescattering corrections (implemented within the process)
      TamingFunctionsCollection taming_functions;

    private:
      std::unique_ptr<Process::GenericProcess> process_;
      std::unique_ptr<Hadroniser::GenericHadroniser> hadroniser_;
      bool store_;
  };
}

#endif
