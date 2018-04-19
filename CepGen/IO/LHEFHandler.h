#ifndef CepGen_Export_LHEFHandler_h
#define CepGen_Export_LHEFHandler_h

#include "CepGen/IO/ExportHandler.h"
#include "CepGen/IO/HepMCHandler.h"

#ifndef LIBHEPMC
# ifndef PYTHIA8
#  pragma message( "HepMC/Pythia8 are not linked to this instance!" )
# endif
#else
# ifndef HEPMC_VERSION3
#  ifdef PYTHIA8
#   include "Pythia8/Pythia.h"
#   define PYTHIA_LHEF 1
#  else
#   pragma message( "HepMC v3 or Pythia8 are required for the LHEF export!" )
#  endif
# else
#  include "HepMC/LHEF.h"
#  define HEPMC_LHEF 1
# endif
#endif

namespace CepGen
{
  class Event;
  namespace OutputHandler
  {
    /**
     * \brief Handler for the LHE file output
     * \author Laurent Forthomme <laurent.forthomme@cern.ch>
     * \date Sep 2016
     */
    class LHEFHandler : public ExportHandler
    {
      public:
        /// Class constructor
        /// \param[in] filename Output file path
        explicit LHEFHandler( const char* filename );
        ~LHEFHandler() override;
        void initialise( const Parameters& ) override;
        /// Writer operator
        void operator<<( const Event& ) override;
        void setCrossSection( double, double ) override;

      private:
#if defined ( HEPMC_LHEF )
        /// Writer object (from HepMC)
        std::unique_ptr<LHEF::Writer> lhe_output_;
        LHEF::HEPRUP run_;
#elif defined ( PYTHIA_LHEF )
        std::unique_ptr<Pythia8::Pythia> pythia_;
        struct LHAevent : Pythia8::LHAup
        {
          explicit LHAevent();
          void initialise( const Parameters& );
          void setCrossSection( unsigned short, double, double );
          void feedEvent( unsigned short proc_id, const Event& ev );
          bool setInit() override { return true; }
          bool setEvent( int ) override { return true; }
          void addComments( const std::string& );
        };
        std::unique_ptr<LHAevent> lhaevt_;
#endif
    };
  }
}

#endif

