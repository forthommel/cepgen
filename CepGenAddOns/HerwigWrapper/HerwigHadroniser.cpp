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
#include "CepGen/EventFilter/EventModifier.h"
#include "CepGen/Modules/EventModifierFactory.h"
#include "CepGen/Physics/Hadroniser.h"
#include "CepGen/Utils/Filesystem.h"

std::shared_ptr<cepgen::Event> kCepGenEvent;               ///< Last event produced by the generator
std::shared_ptr<cepgen::RunParameters> kCepGenParameters;  ///< Generator running parameters

/// Interface to the Herwig hadronisation algorithm
/// \note It can be used in a single particle decay mode as well as a full event hadronisation using the cluster model.
class HerwigHadroniser : public cepgen::hadr::Hadroniser, private Herwig::HerwigUI {
public:
  explicit HerwigHadroniser(const cepgen::ParametersList& params)
      : cepgen::hadr::Hadroniser(params),
        repo_location_(steerPath("herwigPath")),
        generator_(steer<std::string>("generator")),
        repository_(fullPath(steer<std::string>("repository"))),
        in_file_(fullPath("defaults/HerwigDefaults.in")),
        prep_read_dir_(std::vector<std::string>{repo_location_ / "lib"}) {
    ThePEG::Repository::exitOnError() = steer<bool>("exitOnError");
    ThePEG::Repository::load(repository_);
    //Herwig::API::init( *this );
    ThePEG::SamplerBase::setIntegratePerJob(1);
    ThePEG::SamplerBase::setIntegrationJobs(steer<int>("numParallelJobs"));
    CG_INFO("HerwigHadroniser") << "Initialising the Herwig core.\n"
                                << ThePEG::Repository::banner() << "Base path:\n  " << repo_location_ << "\n"
                                << "Repository: " << steer<std::string>("repository");
  }
  virtual ~HerwigHadroniser() {
    if (thepeg_)
      thepeg_->finalize();
    ThePEG::Repository::cleanup();
  }

  static cepgen::ParametersDescription description() {
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
  inline void readString(const std::string& param) override {
    if (const std::string out = ThePEG::Repository::exec(param, std::cerr); !out.empty())
      throw CG_FATAL("HerwigHadroniser") << "Herwig/ThePEG error:\n" << out;
  }
  void initialise() override;
  bool run(cepgen::Event&, double&, bool) override;
  inline void setCrossSection(const cepgen::Value&) override {}
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
    CG_ERROR("HerwigHadroniser") << "An error occured...";
    quit();
  }
  inline void quit() const override {
    ThePEG::Repository::cleanup();
    CG_INFO("HerwigHadroniser") << "Cleanup of the hadroniser";
  }
  inline std::ostream& outStream() const override { return *cepgen::utils::Logger::get().output(); }
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
  const std::string generator_, repository_, in_file_, setup_file_;
  const std::vector<std::string> prep_read_dir_, app_read_dir_;
};

void HerwigHadroniser::initialise() {
  std::cerr.setstate(std::ios_base::badbit);  //FIXME avoid to fill the error stream
  ThePEG::Repository::update();
  if (CG_LOG_MATCH("HerwigHadroniser", debug)) {
    std::ostringstream oss, oss2;
    ThePEG::Repository::stats(oss);
    for (const auto& path : ThePEG::DynamicLoader::allPaths())
      if (path != "." && path != "/")
        oss2 << "\n *) " << path;
    CG_DEBUG("HerwigHadroniser") << "ThePEG configuration:\n"
                                 << "==================================\n"
                                 << oss.str() << "==================================\n"
                                 << "Paths loaded in the dynamic loader:" << oss2.str() << "\n"
                                 << "==================================";
  }
  //--- building the environment
  if (generator_.empty())
    throw CG_FATAL("HerwigHadroniser") << "Empty event generator!";
  try {
    ThePEG::BaseRepository::CheckObjectDirectory(generator_);
    ThePEG::EGPtr tmp = ThePEG::BaseRepository::GetObject<ThePEG::EGPtr>(generator_);
    if (!tmp)
      throw CG_FATAL("HerwigHadroniser") << "Event generator could not be initialised!";
    kCepGenParameters.reset(new cepgen::RunParameters(runParameters()));
    ThePEG::SamplerBase::setRunLevel(ThePEG::SamplerBase::RunMode);
    thepeg_ = ThePEG::Repository::makeRun(tmp, steer<std::string>("run"));
    thepeg_->setSeed((long)seed_);
  } catch (const ThePEG::Exception& e) {
    throw CG_FATAL("HerwigHadroniser") << "ThePEG exception caught:\n\t" << e.what();
  } catch (const cepgen::Exception&) {
    throw;
  } catch (const std::exception& e) {
    throw CG_ERROR("HerwigHadroniser") << "Core exception caught:\n\t" << e.what();
  } catch (...) {
    throw CG_FATAL("HerwigHadroniser") << "Unknown exception caught!";
  }
  CG_INFO("HerwigHadroniser") << "Event generator successfully initialised.";

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

bool HerwigHadroniser::run(cepgen::Event& ev, double& weight, bool) {
  weight = 1.;
  kCepGenEvent.reset(new cepgen::Event(ev));
  try {
    ThePEG::EventPtr evt = thepeg_->shoot();
    ThePEG::tSubProPtr proc = evt->primarySubProcess();
    if (!proc) {
      CG_WARNING("HerwigHadroniser") << "Failed to retrieve the primary subprocess";
      return false;
    }
    proc->printMe(std::cerr);
    for (const auto& ip : proc->collision()->getRemnants())
      std::cout << ip->id() << std::endl;
    exit(0);
    return true;
  } catch (ThePEG::Exception& e) {
    throw CG_FATAL("HerwigHadroniser") << "ThePEG exception caught:\n\t" << e.what();
  } catch (std::exception& e) {
    throw CG_FATAL("HerwigHadroniser") << "Core exception caught:\n\t" << e.what();
  } catch (const char* what) {
    throw CG_FATAL("HerwigHadroniser") << "Other exception caught:\n\t" << what;
  }
  return true;
}

REGISTER_MODIFIER("herwig", HerwigHadroniser);
