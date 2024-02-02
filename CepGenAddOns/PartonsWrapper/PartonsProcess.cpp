/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2024  Laurent Forthomme
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

#include <ElementaryUtils/logger/CustomException.h>
#include <partons/Partons.h>

#include <cstring>

#include "CepGen/Core/Exception.h"
#include "CepGen/Event/Event.h"
#include "CepGen/Modules/ProcessFactory.h"
#include "CepGen/Physics/PDG.h"
#include "CepGen/Process/Process.h"
#include "CepGen/Utils/Environment.h"
#include "CepGen/Utils/Filesystem.h"
#include "CepGen/Utils/String.h"

using namespace cepgen;

static PARTONS::Partons* partons_ = nullptr;

class PartonsProcess final : public cepgen::proc::Process {
public:
  explicit PartonsProcess(const ParametersList& params) : proc::Process(params) {
    if (!partons_)
      partons_ = PARTONS::Partons::getInstance();
    std::vector<std::string> args{fs::canonical("/proc/self/exe").parent_path().parent_path() / "partons" /
                                  "partons.properties"};
    const auto arguments = steer<ParametersList>("arguments");
    for (const auto& key : arguments.keys())
      args.emplace_back(utils::format("--%s=%s", key.data(), arguments.getString(key).data()));
    std::vector<char*> argv;
    std::transform(args.begin(), args.end(), std::back_inserter(argv), [](const std::string& str) -> char* {
      char* out = new char[str.size() + 1];
      std::strcpy(out, str.data());
      return out;
    });
    try {
      partons_->init(args.size(), argv.data());
    } catch (const ElemUtils::CustomException& exc) {
      throw CG_FATAL("PartonsProcess") << "Fatal Partons exception: " << exc.what();
    }
  }

  virtual ~PartonsProcess() { partons_->close(); }

  proc::ProcessPtr clone() const override { return proc::ProcessPtr(new PartonsProcess(*this)); }

  void addEventContent() override {
    proc::Process::setEventContent({{Particle::IncomingBeam1, {PDG::proton}},
                                    {Particle::IncomingBeam2, {PDG::proton}},
                                    {Particle::Parton1, {PDG::photon}},
                                    {Particle::Parton2, {PDG::photon}},
                                    {Particle::OutgoingBeam1, {PDG::proton}},
                                    {Particle::OutgoingBeam2, {PDG::proton}},
                                    {Particle::CentralSystem, {PDG::muon, PDG::muon}}});
  }
  double computeWeight() override { return 0.; }
  void prepareKinematics() override {
    //--- variables mapping
    defineVariable(m_u_t1_, Mapping::linear, {0., 1.}, "u_t1");
  }
  void fillKinematics() override {}

  static ParametersDescription description() {
    auto desc = proc::Process::description();
    desc.setDescription("Partons process");
    return desc;
  }

private:
  // mapped variables
  double m_u_t1_{0.};
};
// register process
REGISTER_PROCESS("partons", PartonsProcess);
