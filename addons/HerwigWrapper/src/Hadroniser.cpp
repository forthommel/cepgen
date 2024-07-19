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

#include <Herwig/API/HerwigAPI.h>
#include <Herwig/API/HerwigUI.h>
#include <ThePEG/EventRecord/Event.h>
#include <ThePEG/Handlers/SamplerBase.h>
#include <ThePEG/Persistency/PersistentIStream.h>
#include <ThePEG/Repository/EventGenerator.h>
#include <ThePEG/Repository/Repository.h>
#include <ThePEG/Utilities/DynamicLoader.h>
#include <ThePEG/Vectors/HepMCTraits.h>

#include <iostream>
#include <memory>
#include <sstream>

#include "CepGen/Core/Exception.h"
#include "CepGen/Core/ParametersList.h"
#include "CepGen/Core/RunParameters.h"
#include "CepGen/Event/Event.h"
#include "CepGen/Modules/EventModifierFactory.h"
#include "CepGen/Physics/Hadroniser.h"
#include "CepGen/Utils/Filesystem.h"

namespace cepgen::herwig {
  std::shared_ptr<Event> kCepGenEvent;  ///< Last event produced by the generator
  RunParameters* kCepGenParameters;     ///< Generator running parameters

  /// Interface to the Herwig hadronisation algorithm
  /// \note It can be used in a single particle decay mode as well as a full event hadronisation using the cluster model.
  class Hadroniser : public cepgen::hadr::Hadroniser, private Herwig::HerwigUI {
  public:
    explicit Hadroniser(const ParametersList& params)
        : cepgen::hadr::Hadroniser(params),
          repo_location_(steerPath("herwigPath")),
          run_(steer<std::string>("run")),
          generator_(steer<std::string>("generator")),
          repository_(fullPath(steer<std::string>("repository"))),
          in_file_(fullPath("defaults/HerwigDefaults.in")),
          prep_read_dir_(std::vector<std::string>{repo_location_ / "lib"}) {
      ThePEG::Repository::exitOnError() = steer<bool>("exitOnError");
      ThePEG::Repository::load(repository_);
      //Herwig::API::init( *this );
      ThePEG::SamplerBase::setIntegratePerJob(1);
      ThePEG::SamplerBase::setIntegrationJobs(steer<int>("numParallelJobs"));
      CG_INFO("herwig:Hadroniser") << "Initialising the Herwig core.\n"
                                   << ThePEG::Repository::banner() << "Base path:\n  " << repo_location_ << "\n"
                                   << "Repository: " << steer<std::string>("repository");
    }
    virtual ~Hadroniser() {
      if (thepeg_)
        thepeg_->finalize();
      ThePEG::Repository::cleanup();
    }

    static ParametersDescription description() {
      auto desc = cepgen::hadr::Hadroniser::description();
      desc.setDescription("Interface to the Herwig C++ utilitaries");
      desc.add<int>("numParallelJobs", 1).setDescription("number of jobs to run in parallel");
      desc.add<bool>("exitOnError", true).setDescription("exit the CepGen run if an error is encountered?");
      desc.add<std::string>("herwigPath", "").setDescription("path to the Herwig installation");
      desc.add<std::string>("run", "").setDescription("name of the Herwig run");
      desc.add<std::string>("generator", "").setDescription("name of the generator to call");
      desc.add<std::string>("repository", "HerwigDefaults.rpo").setDescription("location to the repository");
      return desc;
    }

    /// \name CepGen UI part
    //\{
    void readString(const std::string& param) override {
      if (const std::string out = ThePEG::Repository::exec(param, std::cerr); !out.empty())
        throw CG_FATAL("herwig:Hadroniser") << "Herwig/ThePEG error:\n" << out;
    }
    void initialise() override;
    bool run(Event&, double&, bool) override;
    void setCrossSection(const Value&) override {}
    //\}

    /// \name Herwig UI part
    //\{
    inline Herwig::RunMode::Mode runMode() const override { return run_mode_; }
    inline std::string repository() const override { return repository_; }
    inline std::string inputfile() const override { return in_file_; }
    inline std::string setupfile() const override { return setup_file_; }
    inline bool resume() const override { return false; }
    inline bool tics() const override { return true; }
    inline std::string tag() const override { return ""; }                         //FIXME
    inline std::string integrationList() const override { return "integration"; }  //FIXME
    inline const std::vector<std::string>& appendReadDirectories() const override { return prep_read_dir_; }
    inline const std::vector<std::string>& prependReadDirectories() const override { return app_read_dir_; }

    inline long N() const { return 1l; }
    inline int seed() const { return seed_; }
    inline int jobs() const { return 1; }
    inline unsigned int jobSize() const { return 1; }
    inline unsigned int maxJobs() const { return 1; }
    inline void quitWithHelp() const override {
      CG_ERROR("herwig:Hadroniser") << "An error occured...";
      quit();
    }
    void quit() const override {
      ThePEG::Repository::cleanup();
      CG_INFO("herwig:Hadroniser") << "Cleanup of the hadroniser";
    }
    inline std::ostream& outStream() const override { return *utils::Logger::get().output(); }
    inline std::ostream& errStream() const override { return std::cerr; }
    inline std::istream& inStream() const override { return ss_; }
    //\}

  private:
    inline std::string fullPath(const std::string& path) const { return repo_location_ / "share" / "Herwig" / path; }

    ThePEG::EGPtr thepeg_{nullptr};
    ThePEG::EventPtr evt_;
    Herwig::RunMode::Mode run_mode_{Herwig::RunMode::READ};
    mutable std::stringstream ss_;
    const fs::path repo_location_;
    const std::string run_;
    const std::string generator_, repository_, in_file_, setup_file_;
    const std::vector<std::string> prep_read_dir_, app_read_dir_;
  };

  void Hadroniser::initialise() {
    std::cerr.setstate(std::ios_base::badbit);  //FIXME avoid to fill the error stream
    ThePEG::Repository::update();
    if (CG_LOG_MATCH("herwig:Hadroniser", debug)) {
      std::ostringstream oss, oss2;
      ThePEG::Repository::stats(oss);
      for (const auto& path : ThePEG::DynamicLoader::allPaths())
        if (path != "." && path != "/")
          oss2 << "\n *) " << path;
      CG_DEBUG("herwig:Hadroniser") << "ThePEG configuration:\n"
                                    << "==================================\n"
                                    << oss.str() << "==================================\n"
                                    << "Paths loaded in the dynamic loader:" << oss2.str() << "\n"
                                    << "==================================";
    }
    //--- building the environment
    if (generator_.empty())
      throw CG_FATAL("herwig:Hadroniser") << "Empty event generator!";
    try {
      ThePEG::BaseRepository::CheckObjectDirectory(generator_);
      ThePEG::EGPtr tmp = ThePEG::BaseRepository::GetObject<ThePEG::EGPtr>(generator_);
      if (!tmp)
        throw CG_FATAL("herwig:Hadroniser") << "Event generator could not be initialised!";
      kCepGenParameters = const_cast<RunParameters*>(&runParameters());
      ThePEG::SamplerBase::setRunLevel(ThePEG::SamplerBase::RunMode);
      thepeg_ = ThePEG::Repository::makeRun(tmp, run_);
      thepeg_->setSeed((long)seed_);
    } catch (const Exception&) {
      throw;
    } catch (const ThePEG::Exception& e) {
      throw CG_FATAL("herwig:Hadroniser") << "ThePEG exception caught:\n\t" << e.what();
    } catch (const std::exception& e) {
      throw CG_ERROR("herwig:Hadroniser") << "Core exception caught:\n\t" << e.what();
    } catch (...) {
      throw CG_FATAL("herwig:Hadroniser") << "Unknown exception caught!";
    }
    CG_INFO("herwig:Hadroniser") << "Event generator successfully initialised.";

    /*switch (runParameters().kinematics.mode) {
      case Kinematics::Mode::ElasticElastic:
        break;
      case Kinematics::Mode::ElasticInelastic:
        break;
      case Kinematics::Mode::InelasticElastic:
        break;
      case Kinematics::Mode::InelasticInelastic:
      default:
        break;
    }*/
  }

  bool Hadroniser::run(Event& ev, double& weight, bool) {
    weight = 1.;
    kCepGenEvent.reset(new Event(ev));

    //kCepGenEvent->dump();
    try {
      std::cout << "hihi" << std::endl;
      ThePEG::EventPtr evt = thepeg_->shoot();
      ThePEG::tSubProPtr proc = evt->primarySubProcess();
      if (!proc) {
        CG_WARNING("herwig:Hadroniser") << "Failed to retrieve the primary subprocess";
        return false;
      }
      proc->printMe(std::cerr);
      for (const auto& ip : proc->collision()->getRemnants())
        std::cout << ip->id() << std::endl;
      exit(0);
      return true;
    } catch (ThePEG::Exception& e) {
      throw CG_FATAL("herwig:Hadroniser") << "ThePEG exception caught:\n\t" << e.what();
    } catch (std::exception& e) {
      throw CG_FATAL("herwig:Hadroniser") << "Core exception caught:\n\t" << e.what();
    } catch (const char* what) {
      throw CG_FATAL("herwig:Hadroniser") << "Other exception caught:\n\t" << what;
    }
    return true;
  }
}  // namespace cepgen::herwig
using HerwigHadroniser = cepgen::herwig::Hadroniser;
REGISTER_MODIFIER("herwig", HerwigHadroniser);
