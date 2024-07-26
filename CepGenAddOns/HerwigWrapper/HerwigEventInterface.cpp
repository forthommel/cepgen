/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2018-2024  Laurent Forthomme
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

#include <ThePEG/Interface/InterfacedBase.h>
#include <ThePEG/Interface/Parameter.h>
#include <ThePEG/LesHouches/LesHouches.h>
#include <ThePEG/LesHouches/LesHouchesReader.h>

#include "CepGen/Core/Exception.h"
#include "CepGen/Core/RunParameters.h"
#include "CepGen/Event/Event.h"
#include "CepGen/Physics/Constants.h"
#include "CepGen/Physics/PDG.h"
#include "CepGen/Utils/Value.h"

extern std::shared_ptr<cepgen::Event> kCepGenEvent;
extern std::shared_ptr<cepgen::RunParameters> kCepGenParameters;

namespace ThePEG {
  /// ThePEG/Herwig interface to CepGen run and event structure
  class CepGenInterface : public LesHouchesReader {
  public:
    explicit CepGenInterface() {}
    /// Register the module in a ThePEG repository scope
    static void Init();
    /// Set the cross section for the process
    /// \param[in] xsec Process cross section and uncertainty, in pb
    inline void setCrossSection(const cepgen::Value& xsec) { xsec_ = xsec; }

  protected:
    /// Make a simple clone of this object
    /// \return a pointer to the new object
    inline IBPtr clone() const override { return new_ptr(*this); }
    /// Make a clone of this object, possibly modifying the cloned object to make it sane
    /// \return a pointer to the new object
    inline IBPtr fullclone() const override { return clone(); }

  private:
    void open() override;
    void close() override {}

    void fillParticle(unsigned short id, const cepgen::Particle&);
    void dumpEvent() const;

    double getEvent() override;
    bool doReadEvent() override;
    inline double eventWeight() { return 1.; }
    double reweight() {
      preweight = 1.;
      return preweight;
    }

    inline long scan() override { return 1; }
    inline std::vector<std::string> optWeightsNamesFunc() override { return {"No weight name defined"}; }

    static const double mp_, mp2_;
    bool init_{false};
    cepgen::Value xsec_;
    long printEvery_;
    static ClassDescription<CepGenInterface> initCepGenInterface;
  };

  const double CepGenInterface::mp_ = cepgen::PDG::get().mass(cepgen::PDG::proton);
  const double CepGenInterface::mp2_ = CepGenInterface::mp_ * CepGenInterface::mp_;

  void CepGenInterface::open() {
    auto params = kCepGenParameters;
    if (!params) {
      CG_ERROR("CepGenInterface") << "Parameters block is not set.";
      throw Stop();
    }
    //--- beam particles information
    const auto& pos_beam = params->kinematics().incomingBeams().positive();
    const auto& neg_beam = params->kinematics().incomingBeams().negative();
    heprup.IDBMUP = std::pair<long, long>{(long)pos_beam.integerPdgId(), (long)neg_beam.integerPdgId()};
    heprup.EBMUP = std::pair<double, double>{pos_beam.momentum().pz(), neg_beam.momentum().pz()};
    //heprup.PDFGUP = std::pair<int,int>{ 0, 0 }; //FIXME
    //--- subprocesses considered
    heprup.NPRUP = 1;                       // number of subprocesses
    heprup.LPRUP = {0};                     // subprocess code
    heprup.XSECUP = {(double)xsec_};        // subprocess cross section
    heprup.XERRUP = {xsec_.uncertainty()};  // subprocess cross section error
    heprup.XMAXUP = {1000.};                // maximum event weight
    //--- how the generator envisages the events weights should be interpreted
    heprup.IDWTUP = 3;

    CG_DEBUG("CepGenInterface") << "CepGen parameters successfully parsed!";
  }

  bool CepGenInterface::doReadEvent() {
    reset();

    //--- interfacing magic is done here...
    auto evt = kCepGenEvent;
    if (!evt)
      throw Stop();

    hepeup.XPDWUP = {-1, -1};
    hepeup.IDPRUP = 0;
    hepeup.XWGTUP = 1.;
    hepeup.SCALUP = (*evt)[cepgen::Particle::Intermediate][0].get().momentum().mass();
    hepeup.AQEDUP = cepgen::constants::ALPHA_EM;
    hepeup.AQCDUP = cepgen::constants::ALPHA_QCD;

    const auto& cs = (*evt)[cepgen::Particle::CentralSystem];
    const auto& part1 = (*evt)[cepgen::Particle::Parton1];
    const auto& part2 = (*evt)[cepgen::Particle::Parton2];

    hepeup.resize(part1.size() + part2.size() + cs.size());  // incoming partons + central system
    unsigned short i = 0;
    std::vector<unsigned short> parton1_ids, parton2_ids;

    //--------------------------------------------------------------------------
    // add the partons
    //--------------------------------------------------------------------------

    for (auto p : part1)
      fillParticle(i, p.get().setStatus(4)), parton1_ids.emplace_back(i++);
    for (auto p : part2)
      fillParticle(i, p.get().setStatus(4)), parton2_ids.emplace_back(i++);

    //--------------------------------------------------------------------------
    // add the final state system
    //--------------------------------------------------------------------------

    for (const auto& p : cs)
      fillParticle(i, p), hepeup.MOTHUP[i++] = std::make_pair(parton1_ids.at(0), parton2_ids.at(0));

    //--------------------------------------------------------------------------
    // final sanity check
    //--------------------------------------------------------------------------

    if (i != hepeup.NUP)
      throw CG_FATAL("CepGenInterface") << "Failed to add all particles in HEPEUP block!";

    CG_DEBUG("CepGenInterface") << "Event successfully built!";
    dumpEvent();

    return true;
  }

  void CepGenInterface::fillParticle(unsigned short id, const cepgen::Particle& part) {
    hepeup.IDUP[id] = part.integerPdgId();
    hepeup.ISTUP[id] = (int)part.status();
    hepeup.ICOLUP[id] = std::pair<int, int>{0, 0};  //FIXME
    const auto& mom = part.momentum();
    hepeup.PUP[id] = std::array<double, 5>{mom.px(), mom.py(), mom.pz(), mom.energy(), mom.mass()};
    hepeup.VTIMUP[id] = -1.;
    hepeup.SPINUP[id] = 9;
  }

  double CepGenInterface::getEvent() {
    if (!doReadEvent())
      return 0.;

    if (!init_) {
      if (!(inPDF.first && inPDF.second))
        initPDFs();
      //std::cout << "--> " << getEBeamA() << "|" << getEBeamB() << std::endl;
      //createPartonBinInstances();
      if (!checkPartonBin())
        throw CG_FATAL("CepGenInterface") << "Found event which cannot be handled by the assigned PartonExtractor!";
      init_ = true;
    }

    fillEvent();
    getSubProcess();

    return eventWeight();
  }

  void CepGenInterface::dumpEvent() const {
    CG_INFO("CepGenInterface").log([&](auto& log) {
      log << "Event from subprocess " << hepeup.IDPRUP << "."
          << "\n\tScale: " << hepeup.SCALUP << " GeV, couplings: QED=" << hepeup.AQEDUP << ", QCD=" << hepeup.AQCDUP;
      for (unsigned short i = 0; i < hepeup.NUP; ++i)
        log << "\n"
            << i << ": "
            << "(pdgid=" << std::setw(6) << hepeup.IDUP[i] << ", "
            << "status=" << std::setw(3) << hepeup.ISTUP[i] << ", "
            << "m=" << std::setw(8) << hepeup.PUP[i][4] << ", "
            << "moth=" << hepeup.MOTHUP[i].first << "|" << hepeup.MOTHUP[i].second << ") (p|E)=(" << std::setw(12)
            << std::setprecision(5) << hepeup.PUP[i][0] << "," << std::setw(12) << std::setprecision(5)
            << hepeup.PUP[i][1] << "," << std::setw(12) << std::setprecision(5) << hepeup.PUP[i][2] << "|"
            << std::setw(12) << std::setprecision(5) << hepeup.PUP[i][3] << ") GeV";
    });
  }

  //--- ThePEG plugin registration

  /// Book a description of this library
  ClassDescription<CepGenInterface> CepGenInterface::initCepGenInterface;

  void CepGenInterface::Init() {
    static ClassDocumentation<CepGenInterface> documentation("Interfacing object between CepGen and ThePEG");
    static Parameter<CepGenInterface, long> interfacePrintEvery("PrintEvery",
                                                                "Periodicity of event printout",
                                                                &CepGenInterface::printEvery_,
                                                                -1,
                                                                -1,
                                                                10000,
                                                                false,
                                                                false,
                                                                Interface::lowerlim);  //FIXME not yet implemented
  }

  /// Specify the base class for the CepGen interfacing library
  template <>
  struct BaseClassTrait<CepGenInterface, 1> : public ClassTraitsType {
    /// Alias for the first base class
    typedef LesHouchesReader NthBase;
  };

  /// Information on the CepGen interfacing library
  template <>
  struct ClassTraits<CepGenInterface> : public ClassTraitsBase<CepGenInterface> {
    /// Platform-independent class name
    static std::string className() { return "ThePEG::CepGenInterface"; }
    /// Dynamic library location where the class CepGenInterface is implemented
    static std::string library() { return "libCepGenHerwig.so"; }
  };
}  // namespace ThePEG
