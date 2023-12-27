/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2020-2023  Laurent Forthomme
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

#ifndef CepGenAddOns_EventInterfaces_CepGenEventInterface_h
#define CepGenAddOns_EventInterfaces_CepGenEventInterface_h

#include <Pythia8/Pythia.h>

#include <unordered_map>

namespace cepgen {
  class Beam;
  class Event;
  class Particle;
}  // namespace cepgen

namespace Pythia8 {
  /// Interfacing between CepGen and Pythia8 event definitions
  class CepGenEventInterface : public LHAup {
  public:
    /// List of particles to be included to the event content
    enum Type : unsigned int {
      central = (1 << 0),           ///< include central system
      partonsKT = (1 << 1),         ///< include kt-dependent initiators
      partonsCollinear = (1 << 2),  ///< include purely collinear initiators
      beamRemnants = (1 << 3)       ///< include undissociated beam remnants
    };
    explicit CepGenEventInterface();
    /// Initialise this conversion object with CepGen parameters
    void initialise(const cepgen::Beam& pos_beam, const cepgen::Beam& neg_beam);
    /// Feed a new CepGen event to this conversion object
    /// \param[in] ev CepGen event to be fed
    /// \param[in] type Type of storage
    void feedEvent(const cepgen::Event& ev, unsigned int type);
    /// Set the cross section for a given process
    /// \param[in] id Process identifier
    /// \param[in] cross_section Process cross section, in pb
    /// \param[in] cross_section_err Uncertainty on process cross section, in pb
    void setCrossSection(int id, double cross_section, double cross_section_err);
    /// Specify new process attributes
    /// \param[in] id Process identifier
    /// \param[in] cross_section Process cross section, in pb
    /// \param[in] q2_scale Hard event scale \f$Q^2\f$, in GeV\f$^2\f$
    /// \param[in] alpha_qed \f$\alpha_{\rm em}\f$ for this process
    /// \param[in] alpha_qcd \f$\alpha_{\rm s}\f$ for this process
    void setProcess(int id, double cross_section, double q2_scale, double alpha_qed, double alpha_qcd);
    /// Update the CepGen event content with the Pythia record
    void updateEvent(const Pythia8::Event&, cepgen::Event&, double& weight, bool) const;

    /// Feed comments to the LHEF block
    void addComments(const std::string& comments);

    /// Add a CepGen particle to the event content
    void addCepGenParticle(const cepgen::Particle& part,
                           int status = INVALID_ID,
                           const std::pair<int, int>& mothers = {0, 0},
                           const std::pair<int, int>& colours = {0, 0});
    /// Register a new Pythia8 / CepGen particle mapping
    /// \param[in] py_id Pythia8 particle id
    /// \param[in] cg_id CepGen particle id
    void addCorresp(unsigned short py_id, unsigned short cg_id);
    /// Print all Pythia8 / CepGen Particles correspondences
    void dumpCorresp(std::ostream&) const;

    static constexpr unsigned short INVALID_ID = 999;        ///< Invalid id association
    static constexpr unsigned short MIN_COLOUR_INDEX = 501;  ///< Minimal colour indexing number

    inline bool setInit() override { return true; }
#if defined(PYTHIA_VERSION_INTEGER) && PYTHIA_VERSION_INTEGER >= 8200
    bool setEvent(int) override { return true; }
#else
    bool setEvent(int, double) override { return true; }
#endif

  private:
    static constexpr unsigned short PYTHIA_STATUS_IN_BEAM = 12;
    static constexpr unsigned short PYTHIA_STATUS_IN_PARTON_KT = 61;

    cepgen::Particle& addParticleToEvent(cepgen::Event&, cepgen::Particle::Role, const Particle&, const Vec4&) const;
    /// Check if particle identifier is already in the list of handled PDGs ; if not, add it
    void checkPDGid(const Particle& py_part) const;
    cepgen::Particle::Role findRole(const cepgen::Event&, const Event&, const Particle&) const;
    std::pair<int, int> findMothers(const cepgen::Event&, const cepgen::Particle&) const;
    /// Retrieve the CepGen particle index given its Pythia8 event id
    /// \param[in] py_id Pythia8 particle id
    /// \return CepGen particle id
    unsigned short cepgenId(unsigned short py_id) const;
    /// Retrieve the Pythia8 particle index given its CepGen event id
    /// \param[in] cg_id CepGen particle id
    /// \return Pythia8 particle id
    unsigned short pythiaId(unsigned short cg_id) const;

    const double mp_, mp2_;
    bool inel1_{false}, inel2_{false};
    std::unordered_map<unsigned short, unsigned short> py_cg_corresp_;
    mutable short offset_{0};  //FIXME
    mutable bool first_evt_{true};
  };
}  // namespace Pythia8
#endif
