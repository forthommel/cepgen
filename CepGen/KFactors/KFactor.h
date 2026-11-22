/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2013-2026  Laurent Forthomme
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

#ifndef CepGen_KFactors_KFactor_h
#define CepGen_KFactors_KFactor_h

#include <atomic>
#include <memory>

#include "CepGen/Modules/NamedModule.h"
#include "CepGen/Utils/Limits.h"

namespace cepgen {
  class Coupling;
  class Event;
}  // namespace cepgen

namespace cepgen {
  /// A generic higher-order (QCD) correction \f$K=\sigma_{\rm NLO}/\sigma_{\rm LO}\f$ modelling
  /// \note Evaluated once per phase space point at the integrand level, on the event content filled by the process
  ///  (like the gap survival factor), hence usable with any process carrying an event. Both the cross section and the
  ///  (unweighted) events distributions account for the rescaling. It is a purely inclusive reweighting: it does not
  ///  generate the additional radiation of the real-emission contribution.
  /// \par Scale choice
  ///  The renormalisation scale is \f$\mu=f\,\mu_0\f$, where \f$f\f$ is \a muFactor (for scale variations) and \f$\mu_0\f$
  ///  is selected by \a scale: "fixed" (value of \a mu), "mass" (central system invariant mass), or "pt" (mean
  ///  transverse momentum of the central particles). The scale is never allowed below \a muMin.
  /// \par Validity range
  ///  The correction is only evaluated for a central system invariant mass inside \a mass, and wherever the
  ///  implementation deems itself valid. Elsewhere, \a outOfRange selects between "unity" (K=1), and "reject"
  ///  (the phase space point is discarded). A warning is issued once.
  /// \date Nov 2026
  class KFactor : public NamedModule<KFactor> {
  public:
    explicit KFactor(const ParametersList&);
    ~KFactor() override;

    static ParametersDescription description();

    /// Compute the correction factor for a given event
    /// \return K-factor (1 if out of validity range and policy is "unity", 0 if policy is "reject")
    double operator()(const Event&) const;
    bool singleThreaded() const override;

  protected:
    /// Raise an error if the event content does not match the one the correction is defined for (run once)
    virtual void checkApplicability(const Event&) const {}
    /// Is the correction valid for this event (in addition to the central mass range)?
    virtual bool inValidityRange(const Event&) const { return true; }
    /// Compute the correction at a given renormalisation scale (in GeV)
    virtual double compute(const Event&, double mu) const = 0;

    double alphaS(double mu) const;  ///< Strong coupling at a given scale (GeV)

  public:
    enum class Scale { fixed, mass, pt };  ///< Renormalisation scale base choice

  private:
    double scale(const Event&) const;

    const Scale scale_;
    const double mu_, mu_factor_, mu_min_;
    const Limits mass_range_;
    const bool reject_;
    std::unique_ptr<Coupling> alpha_s_;
    mutable std::atomic<bool> checked_{false}, warned_{false};
  };
}  // namespace cepgen

#endif
