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

#include <Pythia8/PhaseSpace.h>

#include <vector>

#include "CepGen/Core/Exception.h"
#include "CepGen/Modules/PhaseSpaceGeneratorFactory.h"
#include "CepGen/Physics/PDG.h"
#include "CepGen/Process/FactorisedProcess.h"
#include "CepGen/Process/PartonsCollinearPhaseSpaceGenerator.h"
#include "CepGen/Process/PartonsKTPhaseSpaceGenerator.h"
#include "CepGen/Utils/Math.h"

using namespace cepgen;

template <typename T>
class RamboPhaseSpaceGenerator : public cepgen::PhaseSpaceGenerator {
public:
  explicit RamboPhaseSpaceGenerator(const ParametersList& params)
      : PhaseSpaceGenerator(params),
        part_psgen_(new T(params)),
        random_engine_(new Pythia8RandomWrapper(coords_)),
        random_(steer<unsigned long long>("seed")),
        rambo_(&random_),
        central_(steer<std::vector<int> >("ids")) {
    random_.rndmEnginePtr(random_engine_);
    //for (const auto& id : partons_)
    //  masses_.emplace_back(PDG::get().mass(id));
    for (const auto& id : central_)
      masses_.emplace_back(PDG::get().mass(id));
  }

  static ParametersDescription description() {
    auto desc = PhaseSpaceGenerator::description();
    desc.setDescription("Rambo central phase space generator");
    desc.add<unsigned long long>("seed", 42);
    return desc;
  }

  pdgids_t partons() const override {
    CG_ASSERT(part_psgen_);
    return pdgids_t{part_psgen_->positiveFlux().partonPdgId(), part_psgen_->negativeFlux().partonPdgId()};
  }

  void setCentral(const std::vector<int>& central) override { central_ = central; }

  std::vector<int> central() const override { return central_; }

  void initialise(proc::FactorisedProcess* process) override {
    CG_ASSERT(part_psgen_);
    part_psgen_->initialise(process);
    proc_ = process;
    const auto ndim = masses_.size() * 4;
    coords_.resize(ndim);
    for (size_t i = 0; i < ndim; ++i)
      proc_->defineVariable(coords_[i], proc::Process::Mapping::linear, {0., 1.}, "x" + std::to_string(i));
  }

  bool generate() override {
    CG_ASSERT(part_psgen_);
    if (!part_psgen_->generatePartonKinematics())
      return false;
    dynamic_cast<Pythia8RandomWrapper*>(random_engine_.get())->reset();
    central_weight_ = rambo_.genPoint((proc_->q1() + proc_->q2()).mass(), masses_, vecs_);
    if (!utils::positive(central_weight_))
      return false;
    if (vecs_.size() < masses_.size())
      throw CG_FATAL("RamboPhaseSpaceGenerator:generate") << "Not enough momenta were generated by Rambo.";
    const auto convert_momentum = [](const Pythia8::Vec4& vec) {
      return Momentum::fromPxPyPzE(vec.px(), vec.py(), vec.pz(), vec.e());
    };
    for (size_t i = 0; i < central_.size(); ++i)
      proc_->pc(i) = convert_momentum(vecs_.at(i));
    return true;
  }

  double weight() const override {
    const auto fluxes_weight = part_psgen_->fluxes();
    if (!utils::positive(fluxes_weight))
      return 0.;
    return fluxes_weight * central_weight_;
  }

  //FIXME only works for 2-to-4
  double that() const override {
    return 0.5 * ((proc_->q1() - proc_->pc(0)).mass2() + (proc_->q2() - proc_->pc(1)).mass2());
  }

  // FIXME only works for 2-to-4
  double uhat() const override {
    return 0.5 * ((proc_->q1() - proc_->pc(1)).mass2() + (proc_->q2() - proc_->pc(0)).mass2());
  }

private:
  class Pythia8RandomWrapper final : public Pythia8::RndmEngine {
  public:
    explicit Pythia8RandomWrapper(std::vector<double>& coords) : coords_(coords) {}
    void reset() { index_ = 0; }
    double flat() override {
      if (index_ >= coords_.size()) {
        CG_WARNING("Pythia8RandomWrapper:flat")
            << "Coordinate index " << index_ << " exceeds the total number of dimensions, " << coords_.size()
            << ", allocated to the Rambo phase space definition.";
        //FIXME need to ensure the Rambo loop is not exceeding the number of dimensions of
        // the phase space definition ; this requires a bit of treatment to the number of
        // possible trials that can be performed and the multiplicity of external particles
        // to the process.
        return 0.;
      }
      return coords_.at(index_++);
    }

  private:
    std::vector<double>& coords_;
    size_t index_{0};
  };

  std::vector<double> coords_;
  proc::FactorisedProcess* proc_;  //NOT owning

  const std::unique_ptr<PartonsPhaseSpaceGenerator> part_psgen_;
  const Pythia8::RndmEnginePtr random_engine_;
  Pythia8::Rndm random_;
  Pythia8::Rambo rambo_;
  std::vector<int> central_;

  std::vector<double> masses_;
  std::vector<Pythia8::Vec4> vecs_;
  double central_weight_{0.};
};
typedef RamboPhaseSpaceGenerator<cepgen::PartonsKTPhaseSpaceGenerator> KTRambo;
typedef RamboPhaseSpaceGenerator<cepgen::PartonsCollinearPhaseSpaceGenerator> CollRambo;
REGISTER_PSGEN("ktrambo", KTRambo);
REGISTER_PSGEN("collrambo", CollRambo);
