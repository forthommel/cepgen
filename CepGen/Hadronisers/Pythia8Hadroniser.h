#ifndef Pythia8Hadroniser_h
#define Pythia8Hadroniser_h

#include "GenericHadroniser.h"

#ifdef PYTHIA8
#include <Pythia8/Pythia.h>
#include <memory>
#endif

namespace CepGen
{
  class Parameters;
  class Particle;
  class LHAEvent : public Pythia8::LHAup
  {
    public:
      LHAEvent( const Parameters* );
      void feedEvent( const Event& ev, bool full );
      bool setInit() override;
      bool setEvent( int ) override;
      void setCrossSection( int id, double xsec, double xsec_err );
      void setProcess( int id, double xsec, double q2_scale, double alpha_qed, double alpha_qcd );

      unsigned short cgPart( unsigned short py_id ) const;
      unsigned short pyPart( unsigned short cg_id ) const;
      void addCorresp( unsigned short py_id, unsigned short cg_id );
      void dumpCorresp() const;

      static constexpr unsigned short invalid_id = 999;
    private:
      static const double mp_, mp2_;
      void fragmentState( const Particle& p_x, double xbj );
      std::vector<std::pair<unsigned short, unsigned short> > py_cg_corresp_;
      const Parameters* params_;
  };

  namespace Hadroniser
  {
    /**
     * Full interface to the Pythia8 hadronisation algorithm. It can be used in a single particle decay mode as well as a full event hadronisation using the string model, as in Jetset.
     * \brief Pythia8 hadronisation algorithm
     */
    class Pythia8Hadroniser : public GenericHadroniser
    {
      public:
        explicit Pythia8Hadroniser( const Parameters& );
        ~Pythia8Hadroniser();

        bool run( Event& ev, double& weight, bool full ) override;
        void setSeed( long long seed ) override;
        void setCrossSection( double xsec, double xsec_err ) override;

#ifdef PYTHIA8
        bool init();
        void readString( const char* param );
        void readString( const std::string& param ) { readString( param.c_str() ); }
#endif

      private:
        static constexpr unsigned short invalid_idx_ = 999;
        unsigned short max_attempts_;
        std::vector<unsigned short> min_ids_;
        std::map<short,short> py_cg_corresp_, cg_py_corresp_;
#ifdef PYTHIA8
        bool launchPythia( Event& ev );
        void updateEvent( Event& ev, double& weight );
        /// A Pythia8 core to be wrapped
        std::unique_ptr<Pythia8::Pythia> pythia_;
        std::shared_ptr<LHAEvent> lhaevt_;
#endif
    };
  }
}

#endif
