/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2013-2023  Laurent Forthomme
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

#include <unordered_map>

#include "CepGen/Core/Exception.h"
#include "CepGen/Event/Event.h"
#include "CepGen/Modules/EventModifierFactory.h"
#include "CepGen/Parameters.h"
#include "CepGen/Physics/Hadroniser.h"
#include "CepGen/Physics/Kinematics.h"
#include "CepGen/Utils/String.h"
#include "CepGen/Utils/Value.h"
#include "CepGenAddOns/Pythia8Wrapper/CepGenEventInterface.h"

namespace cepgen {
  namespace hadr {
    /// Interface to the Pythia8 hadronisation algorithm
    /// \note It can be used in a single particle decay mode as well as a full event hadronisation using the string model, as in Jetset.
    class Pythia8Hadroniser : public Hadroniser {
    public:
      explicit Pythia8Hadroniser(const ParametersList& plist)
          : Hadroniser(plist),
            pythia_(new Pythia8::Pythia),
            cg_evt_interface_(new Pythia8::CepGenEventInterface),
            correct_central_(steer<bool>("correctCentralSystem")),
            debug_(steer<bool>("debug")),
            debug_lhef_(steer<bool>("debugLHEF")),
            output_config_(steer<std::string>("outputConfig")) {
        if (debug_)
          pythia_->settings.mode("Print:verbosity", 3);
      }

      virtual ~Pythia8Hadroniser() {
        if (!output_config_.empty())
          pythia_->settings.writeFile(output_config_, false);
        if (debug_lhef_)
          cg_evt_interface_->closeLHEF(true);
      }

      static ParametersDescription description();

      void readString(const std::string& param) override {
        if (!pythia_->readString(param))
          throw CG_FATAL("Pythia8Hadroniser") << "The Pythia8 core failed to parse the following setting:\n\t" << param;
      }
      void initialise() override;
      bool run(Event& ev, double& weight, bool fast) override;

      void setCrossSection(const Value& cross_section) override {
        cg_evt_interface_->setCrossSection(0, cross_section, cross_section.uncertainty());
      }

    private:
      void* enginePtr() override { return (void*)pythia_.get(); }

      pdgids_t min_ids_;
      std::unordered_map<short, short> py_cg_corresp_;

      /// Pythia 8 core to be wrapped
      const std::unique_ptr<Pythia8::Pythia> pythia_;
      /// Event interface between CepGen and Pythia
      const std::shared_ptr<Pythia8::CepGenEventInterface> cg_evt_interface_;

      const bool correct_central_;
      const bool debug_, debug_lhef_;
      const std::string output_config_;
      bool res_decay_{true};
      bool enable_hadr_{false};
    };

    void Pythia8Hadroniser::initialise() {
      const auto& kin = runParameters().kinematics();
      const auto& beams = kin.incomingBeams();
      cg_evt_interface_->initialise(beams.positive(), beams.negative());
#if defined(PYTHIA_VERSION_INTEGER) && PYTHIA_VERSION_INTEGER < 8300
      pythia_->setLHAupPtr(cg_evt_interface_.get());
#else
      pythia_->setLHAupPtr(cg_evt_interface_);
#endif
      pythia_->settings.parm("Beams:idA", (long)beams.positive().pdgId());
      pythia_->settings.parm("Beams:idB", (long)beams.negative().pdgId());
      pythia_->settings.mode("Beams:frameType", 5);  // specify we will be using a LHA input
      pythia_->settings.parm("Beams:eCM", beams.sqrtS());
      pythia_->settings.flag("LesHouches:matchInOut", false);
      pythia_->settings.flag("BeamRemnants:primordialKT", false);
      min_ids_ = kin.minimumFinalState();
      if (debug_lhef_)
        cg_evt_interface_->openLHEF("debug.lhe");
      if (pythia_->settings.flag("ProcessLevel:all") != enable_hadr_)
        pythia_->settings.flag("ProcessLevel:all", enable_hadr_);

      if (seed_ == -1ll)
        pythia_->settings.flag("Random:setSeed", false);
      else {
        pythia_->settings.flag("Random:setSeed", true);
        pythia_->settings.mode("Random:seed", seed_);
      }

#if defined(PYTHIA_VERSION_INTEGER) && PYTHIA_VERSION_INTEGER >= 8226
      switch (beams.mode()) {
        case mode::Kinematics::ElasticElastic: {
          pythia_->settings.mode("BeamRemnants:unresolvedHadron", 3);
          pythia_->settings.flag("PartonLevel:MPI", false);
        } break;
        case mode::Kinematics::InelasticElastic: {
          pythia_->settings.mode("BeamRemnants:unresolvedHadron", 2);
          pythia_->settings.flag("PartonLevel:MPI", false);
        } break;
        case mode::Kinematics::ElasticInelastic: {
          pythia_->settings.mode("BeamRemnants:unresolvedHadron", 1);
          pythia_->settings.flag("PartonLevel:MPI", false);
        } break;
        case mode::Kinematics::InelasticInelastic:
        default: {
          pythia_->settings.mode("BeamRemnants:unresolvedHadron", 0);
        } break;
      }
#else
      CG_WARNING("Pythia8Hadroniser") << "Beam remnants framework for this version of Pythia "
                                      << "(" << utils::format("%.3f", pythia_->settings.parm("Pythia:versionNumber"))
                                      << ")\n\t"
                                      << "does not support mixing of unresolved hadron states.\n\t"
                                      << "The proton remnants output might hence be wrong.\n\t"
                                      << "Please update the Pythia version or disable this part.";
#endif
      if (!pythia_->init())
        throw CG_FATAL("Pythia8Hadroniser") << "Failed to initialise the Pythia8 core!\n\t"
                                            << "See the message above for more details.";

      res_decay_ = pythia_->settings.flag("ProcessLevel:resonanceDecays");
      if (correct_central_ && res_decay_)
        CG_WARNING("Pythia8Hadroniser") << "Central system's kinematics correction enabled while resonances are\n\t"
                                        << "expected to be decayed. Please check that this is fully intended.";
      if (debug_lhef_)
        cg_evt_interface_->initLHEF();
    }

    bool Pythia8Hadroniser::run(Event& ev, double& weight, bool fast) {
      // initialise the event weight before running any decay algorithm
      weight = 1.;

      // only launch the event modification if:
      // 1) the full event kinematics (i.e. with remnants) is to be specified,
      // 2) the remnants are to be fragmented, or
      // 3) the resonances are to be decayed.
      if (fast && !res_decay_)
        return true;
      if (!fast && !remn_fragm_ && !res_decay_)
        return true;

      //--- switch full <-> partial event
      if (fast == enable_hadr_) {
        enable_hadr_ = !fast;
        initialise();
      }

      // convert the CepGen event into a custom LHA format
      unsigned int type = Pythia8::CepGenEventInterface::Type::central;
      if (fast)
        type |= Pythia8::CepGenEventInterface::Type::partonsKT;
      else
        type |= Pythia8::CepGenEventInterface::Type::partonsCollinear;
      //type |= Pythia8::CepGenEventInterface::Type::beamRemnants;
      cg_evt_interface_->feedEvent(ev, type);
      if (debug_lhef_ && !fast)
        cg_evt_interface_->eventLHEF();
      cg_evt_interface_->listEvent();

      // launch the hadronisation / resonances decays, and update the event accordingly
      auto& num_hadr_trials = ev.metadata["pythia8:num_hadronisation_trials"];
      num_hadr_trials = 0;
      while (++num_hadr_trials <= max_trials_) {
        // run the hadronisation/fragmentation algorithm
        if (!pythia_->next()) {
          CG_DEBUG("Pythia8Hadroniser") << "Failed hadronisation for event.";
          pythia_->stat();
          continue;
        }
        // from this point, hadronisation is successful
        CG_DEBUG("Pythia8Hadroniser") << "Pythia 8 hadronisation performed successfully.\n\t"
                                      << "Number of trials: " << num_hadr_trials << "/" << max_trials_ << ".\n\t"
                                      << "Particles multiplicity: " << ev.particles().size() << " (CepGen) → "
                                      << pythia_->event.size() << " (Pythia 8).";
        // update the event content with Pythia's output
        cg_evt_interface_->updateEvent(pythia_->event, ev, weight, correct_central_);
        return true;
      }
      CG_WARNING("Pythia8Hadroniser") << "Pythia 8 hadronisation failed after "
                                      << utils::s("attempt", num_hadr_trials, true) << " for this event.";
      return false;
    }

    ParametersDescription Pythia8Hadroniser::description() {
      auto desc = Hadroniser::description();
      desc.setDescription("Interface to the Pythia 8 string hadronisation/fragmentation algorithm");
      desc.add<bool>("correctCentralSystem", false)
          .setDescription("Correct the kinematics of the central system whenever required");
      desc.add<bool>("debugLHEF", false).setDescription("Switch on the dump of each event into a debugging LHEF file");
      desc.add<std::string>("outputConfig", "last_pythia_config.cmd")
          .setDescription("Output filename for a backup of the last Pythia configuration snapshot");
      return desc;
    }
  }  // namespace hadr
}  // namespace cepgen
// register hadroniser
using cepgen::hadr::Pythia8Hadroniser;
REGISTER_MODIFIER("pythia8", Pythia8Hadroniser);
