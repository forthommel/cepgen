#include "CepGen/Hadronisers/GenericHadroniser.h"

#include "CepGen/Event/Event.h"
#include "CepGen/Event/Particle.h"

#include "CepGen/Physics/PDG.h"

#include "CepGen/Core/EventModifierHandler.h"
#include "CepGen/Core/ParametersList.h" //FIXME
#include "CepGen/Core/Exception.h"
#include "CepGen/Core/utils.h"

#include <algorithm>

extern "C"
{
  void hweini_();
  void hwudat_();
  void hwuine_();
  void hwbgen_();
  void hwcfor_();
  void hwcdec_();
  void hwdhad_();
  void hwufne_();
  void upinit_() {}
  void upevnt_() {}
  void hwaend_() { CG_INFO( "Herwig6Hadroniser" ) << "End of run"; }
  void timel_(double) {}
  /// Particles content of the event
  extern struct
  {
    /// Event number
    int nevhep;
    /// Number of particles in the event
    int nhep;
    /// Particles' status code
    int isthep[4000];
    /// Particles' PDG id
    int idhep[4000];
    /// Particles' mothers
    int jmohep[4000][2];
    /// Particles' daughters
    int jdahep[4000][2];
    /// Particles' kinematics, in GeV (px, py, pz, E, M)
    double phep[4000][5];
    /// Primary vertex for the particles
    double vhep[4000][5];
  } hepevt_;
  extern struct
  {
    double avwgt, evwgt, gamwt, tlout;
    double wbigst, wgtmax, wgtsum, wsqsum;
    int idhw[4000], ierror, istat;
    int lwevt, maxer, maxpr, nowgt, nrn[2];
    int numer, numeru, nwgts, gensof;
  } hwevnt_;
}

namespace cepgen
{
  namespace hadr
  {
    /// Herwig 6 hadronisation algorithm
    class Herwig6Hadroniser : public GenericHadroniser
    {
      public:
        explicit Herwig6Hadroniser( const ParametersList& );

        void setParameters( const Parameters& ) override {}
        void readString( const char* param ) override;
        void init() override;
        bool run( Event& ev, double& weight, bool full ) override;

        void setCrossSection( double xsec, double xsec_err ) override {}

      private:
        static const std::unordered_map<Particle::Status,int> kStatusMatchMap;

        void fillParticle( size_t id, const Particle& part );
        unsigned long num_events_;
    };

    Herwig6Hadroniser::Herwig6Hadroniser( const ParametersList& plist ) :
      GenericHadroniser( plist, "herwig6" ), num_events_( 0ul )
    {}

    const std::unordered_map<Particle::Status,int>
    Herwig6Hadroniser::kStatusMatchMap = {
      { Particle::Status::PrimordialIncoming, 101 },
      { Particle::Status::FinalState, 1 },
      { Particle::Status::Unfragmented, 167 },
      { Particle::Status::Undecayed, 2 },
      { Particle::Status::Fragmented, 184 },
      { Particle::Status::Propagator, 3 },
      { Particle::Status::Incoming, 135 },
    };

    void
    Herwig6Hadroniser::readString( const char* param )
    {
    }

    void
    Herwig6Hadroniser::init()
    {
      CG_WARNING( "Herwig6Hadroniser" )
        << "Branching fraction not yet implemented in this hadroniser.\n\t"
        << "You will have to specify manually the multiplication factor according\n\t"
        << "to your list of open channels.";
      hweini_();
      hwudat_();
    }

    bool
    Herwig6Hadroniser::run( Event& ev, double& weight, bool full )
    {
      weight = 1.;
      hepevt_.nevhep = num_events_;
      const auto& evt_comp = ev.compressed();
      hepevt_.nhep = 0;
      for ( const auto& part : evt_comp.particles() )
        fillParticle( hepevt_.nhep++, part );
      std::cout << "--->before:" << hepevt_.nhep << std::endl;

      hwuine_();
      hwbgen_();
      hwcfor_();
      hwcdec_();
      hwdhad_();
      hwufne_();
      if ( hwevnt_.ierror != 0 )
        return false;

      std::cout << "--->after:" << hepevt_.nhep << std::endl;
      num_events_++;
      return true;
    }

    void
    Herwig6Hadroniser::fillParticle( size_t id, const Particle& part )
    {
      hepevt_.isthep[id] = (int)part.status();
      hepevt_.idhep[id] = part.integerPdgId();
      const auto& moth = part.mothers();
      hepevt_.jmohep[id][0] = ( moth.size() > 0 ) ? *moth.begin()+1 : 0; // parent 1
      hepevt_.jmohep[id][1] = ( moth.size() > 1 ) ? *moth.rbegin()+1 : 0; // parent 2
      const auto& daug = part.daughters();
      hepevt_.jdahep[id][0] = ( daug.size() > 0 ) ? *daug.begin()+1 : 0; // daughter 1
      hepevt_.jdahep[id][1] = ( daug.size() > 1 ) ? *daug.rbegin()+1 : 0; // daughter 2
      for ( int j = 0; j < 4; ++j )
        hepevt_.phep[id][j] = part.momentum()[j]; // momentum
      hepevt_.phep[id][4] = part.mass();
      for ( int j = 0; j < 5; ++j )
        hepevt_.vhep[id][j] = 0.; // vertex position
    }
  }
}

// register hadroniser
REGISTER_MODIFIER( "herwig6", Herwig6Hadroniser )

