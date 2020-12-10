#include "CepGen/Physics/Hadroniser.h"
#include "CepGen/Modules/EventModifierFactory.h"

#include "CepGen/Core/ParametersList.h" //FIXME
#include "CepGen/Core/Exception.h"

#include "CepGen/Event/Event.h"
#include "CepGen/Event/Particle.h"

#include "CepGen/Physics/PDG.h"

#include <algorithm>

extern "C"
{
  void hwigin_(); // initialise other common blocks
  void hweini_(); // initialise elementary process
  void hwudat_();
  void hwepro_();
  void hwuine_(); // initialise event
  void hwbgen_(); // generate parton cascades
  void hwcfor_(); // perform cluster formation
  void hwcdec_(); // perform cluster decay
  void hwdhad_(); // perform unstable particles decay
  void hwufne_(); // finalise event
  void hwuepr_(); // print event content
  void hwuinc_(); // compute constants and lookup tables
  void hwhrem_( int&, int& );
  void hwhsct_( int& report, int& firstc, int& jmueo, double& ptjim );
  void upinit_() {}
  void upevnt_() {}
  void hwaend_() { CG_INFO( "Herwig6Hadroniser" ) << "End of run"; }
  void timel_(double& tres) { tres = 1.e10; }
  extern struct
  {
    int ipart1, ipart2;
  } hwbeam_;
  extern struct
  {
    char part1[8], part2[8];
  } hwbmch_;
  extern struct
  {
    double ebeam1, ebeam2, pbeam1, pbeam2;
    int iproc, maxev;
  } hwproc_;
  const int np = 4000;
  /// Particles content of the event
  extern struct
  {
    /// Event number
    int nevhep;
    /// Number of particles in the event
    int nhep;
    /// Particles' status code
    int isthep[np];
    /// Particles' PDG id
    int idhep[np];
    /// Particles' mothers
    int jmohep[np][2];
    /// Particles' daughters
    int jdahep[np][2];
    /// Particles' kinematics, in GeV (px, py, pz, E, M)
    double phep[np][5];
    /// Primary vertex for the particles
    double vhep[np][5];
  } hepevt_;
  const int maxpup = 100;
  extern struct
  {
    int idbmup[2];
    double ebmup[2];
    int pdfgup[2], pdfsup[2], idwtup, nprup;
    double xsecup[maxpup], xerrup[maxpup], xmaxup[maxpup];
    int lprup[maxpup];
  } heprup_;

  extern struct
  {
    double asfixd, clq[6][7], coss, costh, ctmax, disf[2][13];
    double emlst, emmax, emmin, empow, emsca, epoln[3];
    double gcoef[7], gpoln, omega0, phomas, ppoln[3];
    double ptmax, ptmin, ptpow;
    double q2max, q2min, q2pow;
    double q2wwmn, q2wwmx;
    double qlim, sins, thmax, y4jt, tmnisr, tqwt;
    double xx[2], xlmin, xxmin;
    double ybmax, ybmin, yjmax, yjmin;
    double ywwmax, ywwmin;
    double whmin, zjmax, zmxisr;
    int iaphig, ibrn[2], ibsh, ico[10], idcmf, idn[10], iflmax, iflmin, ihpro, ipro;
    int mapq, maxfl, bgshat;
    int colisr, fstevt, fstwgt, genev, hvfcen, tpol, durham; // logical
  } hwhard_;
  extern struct
  {
    double avwgt, evwgt, gamwt, tlout;
    double wbigst, wgtmax, wgtsum, wsqsum;
    int idhw[np], ierror, istat;
    int lwevt, maxer, maxpr, nowgt, nrn[2];
    int numer, numeru, nwgts, gensof;
  } hwevnt_;
}

namespace cepgen
{
  namespace hadr
  {
    /// Herwig 6 hadronisation algorithm
    class Herwig6Hadroniser : public Hadroniser
    {
      public:
        explicit Herwig6Hadroniser( const ParametersList& );

        void setParameters( const Parameters& ) override {}
        void readString( const char* param ) override {}
        void init() override;
        bool run( Event& ev, double& weight, bool full ) override;

        void setCrossSection( double xsec, double xsec_err ) override {}

      private:
        static const std::unordered_map<Particle::Status,int> kStatusMatchMap;

        void fillParticle( size_t id, const Particle& part );
        unsigned long num_events_;
    };

    Herwig6Hadroniser::Herwig6Hadroniser( const ParametersList& plist ) :
      Hadroniser( plist ), num_events_( 0ul )
    {
      hwudat_();
      hwhard_.ibrn[0] = seed_;
      hwhard_.ibrn[1] = 2*seed_;
    }

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
    Herwig6Hadroniser::init()
    {
      hwuinc_();
      heprup_.idbmup[0] = heprup_.idbmup[1] = 2212;
      heprup_.ebmup[0] = heprup_.ebmup[1] = 6500.;
      heprup_.pdfgup[0] = heprup_.pdfgup[1] = -1;
      heprup_.pdfsup[0] = heprup_.pdfsup[1] = 0;
      heprup_.idwtup = 0;
      heprup_.nprup = 1;
      heprup_.xsecup[0] = 1.;
      heprup_.xerrup[0] = 0.;
      heprup_.xmaxup[0] = 1.;
      heprup_.lprup[0] = 0;
      /*const std::string p( "P" );
      std::copy( p.begin(), p.end(), hwbmch_.part1 );
      std::copy( p.begin(), p.end(), hwbmch_.part2 );
      //hwproc_.iproc = 1500;*/
      hwproc_.iproc = 0;
      CG_WARNING("");
      hwigin_();
      CG_WARNING("");
      hweini_();
      CG_WARNING( "Herwig6Hadroniser" )
        << "Branching fraction not yet implemented in this hadroniser.\n\t"
        << "You will have to specify manually the multiplication factor according\n\t"
        << "to your list of open channels.";
    }

    bool
    Herwig6Hadroniser::run( Event& ev, double& weight, bool full )
    {
      hwevnt_.avwgt = weight * 1.e-3;
      hepevt_.nevhep = num_events_;
      const auto& evt_comp = ev.compress();
      hepevt_.nhep = 0;
      for ( const auto& part : evt_comp.particles() )
        fillParticle( hepevt_.nhep++, part );
      std::cout << "--->before:" << hepevt_.nhep << std::endl;
      hwuepr_();
      hwhard_.genev = true;
      CG_WARNING("");
      hwuine_();
      CG_WARNING("");
      /*int report, firstc, jmueo;
      double ptjim;
      hwhsct_( report, firstc, jmueo, ptjim );
      CG_INFO("")<<report<<":"<<firstc<<":"<<jmueo<<":"<<ptjim;*/
      //hwbgen_();
      /*int ib1, ib2;
      hwhrem_( ib1, ib2 );
      CG_INFO("")<<ib1<<":"<<ib2;*/
      hwepro_();
      CG_WARNING("");
      hwcfor_();
      hwcdec_();
      /*hwdhad_();
      hwufne_();*/
      if ( hwevnt_.ierror != 0 )
        return false;
      hwuepr_();

      std::cout << "--->after:" << hepevt_.nhep << std::endl;
      num_events_++;
      weight = hwevnt_.avwgt * 1.e3;
      return true;
    }

    void
    Herwig6Hadroniser::fillParticle( size_t id, const Particle& part )
    {
      hepevt_.isthep[id] = (int)part.status();
      /*//if ( part.role() == Particle::Role::IncomingBeam1 )
      if ( part.role() == Particle::Role::OutgoingBeam1 )
      //if ( part.role() == Particle::Role::Parton1 )
        hepevt_.isthep[id] = 147;
      //if ( part.role() == Particle::Role::IncomingBeam2 )
      if ( part.role() == Particle::Role::OutgoingBeam2 )
      //if ( part.role() == Particle::Role::Parton2 )
        hepevt_.isthep[id] = 148;*/
      hepevt_.idhep[id] = part.integerPdgId();
      const auto& moth = part.mothers();
      hepevt_.jmohep[id][0] = ( moth.size() > 0 ) ? *moth.begin()+1 : 0; // parent 1
      hepevt_.jmohep[id][1] = ( moth.size() > 1 ) ? *moth.rbegin()+1 : 0; // parent 2
      const auto& daug = part.daughters();
      hepevt_.jdahep[id][0] = ( daug.size() > 0 ) ? *daug.begin()+1 : 0; // daughter 1
      hepevt_.jdahep[id][1] = ( daug.size() > 1 ) ? *daug.rbegin()+1 : 0; // daughter 2
      const auto& mom_vec = part.momentum().pVector();
      for ( int j = 0; j < 4; ++j )
        hepevt_.phep[id][j] = mom_vec.at( j ); // momentum
      hepevt_.phep[id][4] = part.mass();
      for ( int j = 0; j < 5; ++j )
        hepevt_.vhep[id][j] = 0.; // vertex position
    }
  }
}
// register hadroniser
REGISTER_MODIFIER( "herwig6", Herwig6Hadroniser )

