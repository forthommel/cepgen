/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2016-2025  Laurent Forthomme
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

#ifndef CepGenPythia8_CepGenEvent_h
#define CepGenPythia8_CepGenEvent_h

#include <Pythia8/Pythia.h>

#include <unordered_map>

namespace cepgen {
  class RunParameters;
  class Event;
  class Particle;
}  // namespace cepgen

namespace Pythia8 {
  /// Interfacing between CepGen and Pythia8 event definitions
  class CepGenEvent final : public LHAup {
  public:
    explicit CepGenEvent();
    void initialise(const cepgen::RunParameters&);  ///< Initialise this conversion object with CepGen parameters

    enum struct Type {
      centralAndPartons,          ///< only include initiators and central system
      centralAndBeamRemnants,     ///< include non-dissociated beam remnants and central system
      centralAndFullBeamRemnants  ///< include dissociated beam remnants and central system
    };  ///< List of particles to be included to the event content
    /// Feed a new CepGen event to this conversion object
    /// \param[in] ev CepGen event to be fed
    /// \param[in] type Type of storage
    /// \param[in] event_weight Weight associated to this event in the LHA record
    /// \return false if the event could not be mapped onto a collinear LHA record
    bool feedEvent(const cepgen::Event& ev, const Type& type, double event_weight = 1.);
    /// Lorentz transformation bringing the (collinear) Les Houches central system back to the CepGen (kT-dependent) one
    /// \note Identity unless the last event was fed with the Type::centralAndBeamRemnants mode
    inline const RotBstMatrix& lhaToCepGenBoost() const { return boost_lha_to_cg_; }
    inline void setWeightedEvents() { setStrategy(4); }  ///< Events have a mean weight equal to the cross-section
    /// Set the cross-section for a given process
    /// \param[in] id Process identifier
    /// \param[in] cross_section Process cross-section, in pb
    /// \param[in] cross_section_err Uncertainty on process cross-section, in pb
    void setCrossSection(int id, double cross_section, double cross_section_err);
    /// Specify new process attributes
    /// \param[in] id Process identifier
    /// \param[in] cross_section Process cross-section, in pb
    /// \param[in] q2_scale Hard event scale \f$Q^2\f$, in GeV\f$^2\f$
    /// \param[in] alpha_qed \f$\alpha_{\rm em}\f$ for this process
    /// \param[in] alpha_qcd \f$\alpha_{\rm s}\f$ for this process
    void setProcess(int id, double cross_section, double q2_scale, double alpha_qed, double alpha_qcd);

    void addComments(const std::string& comments);  ///< Feed comments to the LHEF block

    /// Retrieve the CepGen particle index given its Pythia8 event id
    /// \param[in] pythia_id Pythia8 particle id
    /// \return CepGen particle id
    unsigned short cepgenId(unsigned short pythia_id) const;
    /// Retrieve the Pythia8 particle index given its CepGen event id
    /// \param[in] cepgen_id CepGen particle id
    /// \return Pythia8 particle id
    unsigned short pythiaId(unsigned short cepgen_id) const;
    /// Add a CepGen particle to the event content
    void addCepGenParticle(const cepgen::Particle& part,
                           int status = INVALID_ID,
                           const std::pair<int, int>& mothers = {0, 0},
                           const std::pair<int, int>& colours = {0, 0});
    /// Register a new Pythia8 / CepGen particle mapping
    /// \param[in] pythia_id Pythia8 particle id
    /// \param[in] cepgen_id CepGen particle id
    void addCorresp(unsigned short pythia_id, unsigned short cepgen_id);
    void dumpCorresp() const;  ///< Print all Pythia8 / CepGen Particles correspondences

    static constexpr unsigned short INVALID_ID = 999;        ///< Invalid id association
    static constexpr unsigned short MIN_COLOUR_INDEX = 501;  ///< Minimal colour indexing number

    inline bool setInit() override { return true; }
#if defined(PYTHIA_VERSION_INTEGER) && PYTHIA_VERSION_INTEGER >= 8200
    bool setEvent(int) override { return true; }
#else
    bool setEvent(int, double) override { return true; }
#endif

  private:
    /// Build a collinear initial state from kT-dependent partons emitted by the beams
    ///
    /// \note Pythia 8 requires massless, collinear hard-process initiators, while space-like partons are generated in CepGen.
    ///   Two strategies:
    ///   - dissociative beam: kT absorbed by explicit recoil -> incoming massless valence quark \f$k_i = x_i\cdot P_i\f$,
    ///     outgoing massless quark \f$k_i - q_i\f$ -> untouched parton \f$q_i\f$ in hard system.
    ///     \f$(k_i - q_i)^2 = 0\f\Rightarrow x_i = q_i^2 / (P_i^+ q_i^-)\f\rightarrow q_i^+/P_i^+\f$ (\f$k_T \to 0\f$).
    ///   - elastic beam: collinear parton, \f$x_i\cdot P_i\f$, and left to Pythia as unresolved hadron.
    ///     \f$x_i\f$ such that hard system keeps central-system invariant mass, and central system Lorentz-transformed to
    ///     collinear frame (transformation stored for the way back).
    ///   Energy-momentum is conserved exactly in the LHA record in all beam configurations.
    bool feedCollinearInitialState(const cepgen::Event& ev, double scale, unsigned short& colour_index);
    std::pair<int, int> findMothers(const cepgen::Event& cepgen_event, const cepgen::Particle& cepgen_particle) const;
    const double mp_, mp2_;
    bool inel1_{false}, inel2_{false};
    std::unordered_map<unsigned short, unsigned short> py_cg_corresp_;
    const cepgen::RunParameters* params_{nullptr};  //NOT owning
    RotBstMatrix boost_lha_to_cg_;                  ///< collinear Les Houches frame -> CepGen frame (central system)
    RotBstMatrix boost_cg_to_lha_;                  ///< inverse transformation
    bool collinearised_{false};
  };
}  // namespace Pythia8
#endif
