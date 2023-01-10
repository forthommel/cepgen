/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2013-2022  Laurent Forthomme
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

#ifndef CepGen_Core_GeneratorWorker_h
#define CepGen_Core_GeneratorWorker_h

#include <memory>
#include <thread>
#include <vector>

#include "CepGen/Event/Event.h"
#include "CepGen/Generator.h"

namespace cepgen {
  class Integrator;
  class Parameters;
  class ProcessIntegrand;
  class GridParameters;

  /// Monte-Carlo generator instance
  class GeneratorWorker {
  public:
    /// Book the memory slots and structures for the generator
    /// \param[in] integ Integrator instance handled by the mother generator
    explicit GeneratorWorker(const Parameters*, const Integrator*);
    virtual ~GeneratorWorker();

    void setGrid(const GridParameters* grid) { grid_ = grid; }
    ThreadSafe<GridParameters> grid();
    const GridParameters* grid() const { return grid_; }

    /// Launch the event generation
    /// \param[in] num_events Number of events to generate
    /// \param[in] callback The callback function applied on every event generated
    void generate(size_t num_events = 0, Event::callback callback = nullptr);
    std::thread thread(size_t num_events = 0, Event::callback callback = nullptr) {
      return std::thread([=] { generate(num_events, callback); });
    }
    /// Generate a single event
    /// \param[in] callback The callback function applied on every event generated
    bool next(Event::callback callback = nullptr);
    /// Function evaluator
    ProcessIntegrand& integrand() { return *integrand_; }
    /// Function evaluator
    const ProcessIntegrand& integrand() const { return *integrand_; }
    /// Number of events generated by this worker
    size_t numGeneratedEvents() const { return num_generated_events_; }

  private:
    /// Placeholder for invalid bin indexing
    static constexpr int UNASSIGNED_BIN = -999;

    /// Store the event in the output file
    /// \param[in] callback The callback function for every event generated
    /// \return A boolean stating whether or not the event was successfully saved
    bool storeEvent(Event::callback);
    /// Apply a correction cycle to the grid
    bool correctionCycle(bool&);
    /// Prepare the object for event generation
    void computeGenerationParameters();

    /// Steering parameters for the event generation
    /// \note NOT owning
    const Parameters* params_{nullptr};
    /// Pointer to the mother-handled integrator instance
    /// \note NOT owning
    const Integrator* integrator_{nullptr};
    /// Pointer to the generation grid
    /// \note NOT owning
    const GridParameters* grid_{nullptr};
    /// Local event weight evaluator
    std::unique_ptr<ProcessIntegrand> integrand_;
    /// Selected bin at which the function will be evaluated
    int ps_bin_{UNASSIGNED_BIN};  ///< Last bin to be corrected
    std::vector<double> coords_;  ///< Phase space coordinates being evaluated
    size_t num_generated_events_{0};
  };
}  // namespace cepgen

#endif
