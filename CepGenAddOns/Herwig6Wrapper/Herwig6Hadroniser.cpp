#include <algorithm>

#include "CepGen/Core/Exception.h"
#include "CepGen/Core/ParametersList.h"  //FIXME
#include "CepGen/Event/Event.h"
#include "CepGen/Event/Particle.h"
#include "CepGen/Modules/EventModifierFactory.h"
#include "CepGen/Physics/Hadroniser.h"
#include "CepGen/Physics/PDG.h"
#include "CepGenAddOns/Herwig6Wrapper/HerwigInterface.h"

extern "C" {
extern void herwig_init_();
}
namespace cepgen {
  namespace hadr {
    /// Herwig 6 hadronisation algorithm
    class Herwig6Hadroniser : public Hadroniser {
    public:
      explicit Herwig6Hadroniser(const ParametersList&);

      void setRuntimeParameters(const Parameters&) override {}
      void readString(const char* param) override {}
      void init() override;
      bool run(Event& ev, double& weight, bool full) override;

      void setCrossSection(double xsec, double xsec_err) override {}

    private:
      static const std::unordered_map<Particle::Status, int> kStatusMatchMap;

      void fillParticle(size_t id, const Particle& part);
      unsigned long num_events_;
    };

    Herwig6Hadroniser::Herwig6Hadroniser(const ParametersList& plist) : Hadroniser(plist), num_events_(0ul) {
      herwig_init_();
      //hwudat_();
      hwhard_.ibrn[0] = seed_;
      hwhard_.ibrn[1] = 2 * seed_;
    }

    const std::unordered_map<Particle::Status, int> Herwig6Hadroniser::kStatusMatchMap = {
        {Particle::Status::PrimordialIncoming, 101},
        {Particle::Status::FinalState, 1},
        {Particle::Status::Unfragmented, 167},
        {Particle::Status::Undecayed, 2},
        {Particle::Status::Fragmented, 184},
        {Particle::Status::Propagator, 3},
        {Particle::Status::Incoming, 135},
    };

    void Herwig6Hadroniser::init() {
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
      hwuinc_();
      /*const std::string p( "P" );
      std::copy( p.begin(), p.end(), hwbmch_.part1 );
      std::copy( p.begin(), p.end(), hwbmch_.part2 );
      //hwproc_.iproc = 1500;*/
      hwproc_.iproc = 0;
      CG_WARNING("");
      hwigin_();
      CG_WARNING("");
      hweini_();
      CG_WARNING("Herwig6Hadroniser") << "Branching fraction not yet implemented in this hadroniser.\n\t"
                                      << "You will have to specify manually the multiplication factor according\n\t"
                                      << "to your list of open channels.";
    }

    bool Herwig6Hadroniser::run(Event& ev, double& weight, bool full) {
      hwevnt_.avwgt = weight * 1.e-3;
      hepevt_.nevhep = num_events_;
      const auto& evt_comp = ev.compress();
      hepevt_.nhep = 0;
      for (const auto& part : evt_comp.particles())
        fillParticle(hepevt_.nhep++, part);
      CG_LOG << "--->before:" << hepevt_.nhep;
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
      if (hwevnt_.ierror != 0)
        return false;
      hwuepr_();

      CG_LOG << "--->after:" << hepevt_.nhep;
      num_events_++;
      weight = hwevnt_.avwgt * 1.e3;
      return true;
    }

    void Herwig6Hadroniser::fillParticle(size_t id, const Particle& part) {
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
      hepevt_.jmohep[id][0] = (moth.size() > 0) ? *moth.begin() + 1 : 0;   // parent 1
      hepevt_.jmohep[id][1] = (moth.size() > 1) ? *moth.rbegin() + 1 : 0;  // parent 2
      const auto& daug = part.daughters();
      hepevt_.jdahep[id][0] = (daug.size() > 0) ? *daug.begin() + 1 : 0;   // daughter 1
      hepevt_.jdahep[id][1] = (daug.size() > 1) ? *daug.rbegin() + 1 : 0;  // daughter 2
      const auto& mom_vec = part.momentum().pVector();
      for (int j = 0; j < 4; ++j)
        hepevt_.phep[id][j] = mom_vec.at(j);  // momentum
      hepevt_.phep[id][4] = part.mass();
      for (int j = 0; j < 5; ++j)
        hepevt_.vhep[id][j] = 0.;  // vertex position
    }
  }  // namespace hadr
}  // namespace cepgen
// register hadroniser
REGISTER_MODIFIER("herwig6", Herwig6Hadroniser)
