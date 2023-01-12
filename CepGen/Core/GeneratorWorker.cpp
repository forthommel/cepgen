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

#include "CepGen/Core/Exception.h"
#include "CepGen/Core/GeneratorWorker.h"
#include "CepGen/Event/Event.h"
#include "CepGen/EventFilter/EventExporter.h"
#include "CepGen/EventFilter/EventModifier.h"
#include "CepGen/Generator.h"
#include "CepGen/Integration/GridParameters.h"
#include "CepGen/Integration/Integrator.h"
#include "CepGen/Integration/ProcessIntegrand.h"
#include "CepGen/Parameters.h"
#include "CepGen/Process/Process.h"
#include "CepGen/Utils/AbortHandler.h"
#include "CepGen/Utils/ProgressBar.h"
#include "CepGen/Utils/String.h"
#include "CepGen/Utils/TimeKeeper.h"

namespace cepgen {
  GeneratorWorker::GeneratorWorker(const Parameters* params, const Integrator* integr)
      : params_(params), integrator_(integr), integrand_(new ProcessIntegrand(params_)) {
    if (!integrator_)
      throw CG_FATAL("GeneratorWorker") << "Invalid integrator object at 0x" << integrator_ << " specified!";

    coords_.resize(integrand_->size());

    CG_DEBUG("GeneratorWorker:integrator")
        << "New generator worker initialised for integration/event generation.\n\t"
        << "Parameters at " << (void*)params_ << ".\n\t"
        << "Dim-" << integrand_->size() << " " << integrator_->name() << " integrator.";
  }

  GeneratorWorker::~GeneratorWorker() {
    CG_DEBUG("GeneratorWorker") << "Generator worker destructed. Releasing the parameters at " << (void*)params_ << ".";
  }

  //-----------------------------------------------------------------------------------------------
  // events generation part
  //-----------------------------------------------------------------------------------------------

  std::thread GeneratorWorker::thread(size_t num_events, Event::callback callback) {
    utils::gRunMode = 1;
    return std::thread([=] { generate(num_events, callback); });
  }

  void GeneratorWorker::generate(size_t num_events, Event::callback callback) {
    if (!params_)
      throw CG_FATAL("GeneratorWorker:generate") << "[thread#" << std::this_thread::get_id() << "] "
                                                 << "No steering parameters specified!";

    if (num_events < 1)
      num_events = params_->generation().maxGen();

    while (params_->numGeneratedEvents() < num_events && (int)utils::gRunMode >= 0)
      next(callback);
  }

  bool GeneratorWorker::next(Event::callback callback) {
    CG_TICKER(ThreadSafe<Parameters>(params_)->timeKeeper());

    if (!integrator_)
      throw CG_FATAL("GeneratorWorker:generate") << "[thread#" << std::this_thread::get_id() << "] "
                                                 << "No integrator object handled!";
    if (!grid_)
      throw CG_FATAL("GeneratorWorker:generate") << "[thread#" << std::this_thread::get_id() << "] "
                                                 << "No grid object handled!";

    // apply correction cycles if required from previous event
    if (ps_bin_ != UNASSIGNED_BIN) {
      bool store = false;
      while (!correctionCycle(store)) {
      }
      if (store)
        return storeEvent(callback);
    }

    //--- normal generation cycle

    double weight = 0.;
    while ((int)utils::gRunMode >= 0) {
      double y = -1.;
      // select a function value and reject if fmax is too small
      do {
        // ...
        ps_bin_ = integrator_->uniform(0., grid_->size());
        y = integrator_->uniform(0., grid_->globalMax());
        { ThreadSafe<GridParameters>(grid_)->increment(ps_bin_); }
      } while (y > grid_->maxValue(ps_bin_));
      // shoot a point x in this bin
      grid_->shoot(integrator_, ps_bin_, coords_);
      // get weight for selected x value
      weight = integrator_->eval(*integrand_, coords_);
      if (weight > y)
        break;
    }

    if (weight > grid_->maxValue(ps_bin_)) {
      // if weight is higher than local or global maximum,
      // init correction cycle for the next event
      ThreadSafe<GridParameters>(grid_)->initCorrectionCycle(ps_bin_, weight);
    } else  // no grid correction needed for this bin
      ps_bin_ = UNASSIGNED_BIN;

    // return with an accepted event
    return storeEvent(callback);
  }

  bool GeneratorWorker::correctionCycle(bool& store) {
    CG_TICKER(ThreadSafe<Parameters>(params_)->timeKeeper());

    CG_DEBUG_LOOP("GeneratorWorker:correction") << "[thread#" << std::this_thread::get_id() << "] "
                                                << "Correction cycles are started.\n\t"
                                                << "bin = " << ps_bin_ << "\n\t"
                                                << "correction value = " << grid_->correctionValue(ps_bin_) << ".";

    if (grid_->correctionValue(ps_bin_) >= 1.) {
      ThreadSafe<GridParameters>(grid_)->setCorrectionValue(ps_bin_, grid_->correctionValue(ps_bin_) - 1.);
    }

    if (integrator_->uniform() < grid_->correctionValue(ps_bin_)) {
      { ThreadSafe<GridParameters>(grid_)->setCorrectionValue(ps_bin_, -1.); }
      // select x values in phase space bin
      grid_->shoot(integrator_, ps_bin_, coords_);
      const double weight = integrator_->eval(*integrand_, coords_);
      // parameter for correction of correction
      { ThreadSafe<GridParameters>(grid_)->rescale(ps_bin_, weight); }
      // accept event
      if (weight >= grid_->maxHist(integrator_, ps_bin_)) {
        store = true;
        return true;
      }
      return false;
    }
    // correction if too big weight is found while correction
    // (all your bases are belong to us...)
    bool correct = false;
    { correct = ThreadSafe<GridParameters>(grid_)->correct(ps_bin_); }
    return correct;
  }

  bool GeneratorWorker::storeEvent(Event::callback callback) {
    CG_TICKER(ThreadSafe<Parameters>(params_)->timeKeeper());

    if (!integrand_->process().hasEvent())
      return true;

    const auto ngen = params_->numGeneratedEvents();
    if ((ngen + 1) % params_->generation().printEvery() == 0)
      CG_INFO("GeneratorWorker:store") << "[thread#" << std::this_thread::get_id() << "] "
                                       << utils::s("event", ngen + 1, true) << " generated.";
    {
      std::lock_guard<std::mutex>(Generator::mutex);
      if (callback)
        callback(integrand_->process().event(), ngen);
      for (auto& mod : params_->eventExportersSequence())
        *mod << integrand_->process().event();
    }
    { ThreadSafe<Parameters>(params_)->addGenerationTime(integrand_->process().event().time_total); }
    num_generated_events_++;
    return true;
  }

  //-----------------------------------------------------------------------------------------------
  // initial preparation run before the generation of unweighted events
  //-----------------------------------------------------------------------------------------------

  void GeneratorWorker::computeGenerationParameters() {
    if (!grid_)
      throw CG_FATAL("GeneratorWorker:setGen") << "No generator grid object specified.";
    if (!params_)
      throw CG_FATAL("GeneratorWorker:setGen") << "No steering parameters specified.";

    integrand_->setStorage(false);

    CG_INFO("GeneratorWorker:setGen") << "Preparing the grid ("
                                      << utils::s("point", params_->generation().numPoints(), true) << "/bin) "
                                      << "for the generation of unweighted events.";

    const double inv_num_points = 1. / params_->generation().numPoints();
    std::vector<double> point_coord(integrand_->size(), 0.);
    if (point_coord.size() < grid_->n(0).size())
      throw CG_FATAL("GridParameters:shoot") << "Coordinates vector multiplicity is insufficient!";

    // ...
    double sum = 0., sum2 = 0., sum2p = 0.;

    utils::ProgressBar prog_bar(grid_->size(), 5);

    //--- main loop
    for (unsigned int i = 0; i < grid_->size(); ++i) {
      double fsum = 0., fsum2 = 0.;
      for (size_t j = 0; j < params_->generation().numPoints(); ++j) {
        grid_->shoot(integrator_, i, point_coord);
        const double weight = integrator_->eval(*integrand_, point_coord);
        { ThreadSafe<GridParameters>(grid_)->setValue(i, weight); }
        fsum += weight;
        fsum2 += weight * weight;
      }
      const double av = fsum * inv_num_points, av2 = fsum2 * inv_num_points;
      const double sig2 = av2 - av * av;
      sum += av;
      sum2 += av2;
      sum2p += sig2;

      // per-bin debugging loop
      CG_DEBUG_LOOP("GeneratorWorker:setGen").log([&](auto& dbg) {
        const double sig = sqrt(sig2);
        const double eff = (grid_->maxValue(i) != 0.) ? av / grid_->maxValue(i) : 0.;
        dbg << "n-vector for bin " << i << ": " << grid_->n(i) << "\n\t"
            << "av   = " << av << "\n\t"
            << "sig  = " << sig << "\n\t"
            << "fmax = " << grid_->maxValue(i) << "\n\t"
            << "eff  = " << eff;
      });
      prog_bar.update(i + 1);
    }  // end of main loop

    const double inv_max = 1. / grid_->size();
    sum *= inv_max;
    sum2 *= inv_max;
    sum2p *= inv_max;

    const double sig = sqrt(sum2 - sum * sum), sigp = sqrt(sum2p);

    double eff1 = 0.;
    for (unsigned int i = 0; i < grid_->size(); ++i)
      eff1 += sum / grid_->size() * grid_->maxValue(i);
    const double eff2 = sum / grid_->globalMax();

    CG_DEBUG("GeneratorWorker:setGen") << "Average function value         = " << sum << "\n\t"
                                       << "Average squared function value = " << sum2 << "\n\t"
                                       << "Overall standard deviation     = " << sig << "\n\t"
                                       << "Average standard deviation     = " << sigp << "\n\t"
                                       << "Maximum function value         = " << grid_->globalMax() << "\n\t"
                                       << "Average inefficiency           = " << eff1 << "\n\t"
                                       << "Overall inefficiency           = " << eff2;
    { ThreadSafe<GridParameters>(grid_)->setPrepared(true); }
    //--- from now on events will be stored
    integrand_->setStorage(true);

    CG_INFO("GeneratorWorker:setGen") << "Grid prepared! Now launching the production.";
  }
}  // namespace cepgen
