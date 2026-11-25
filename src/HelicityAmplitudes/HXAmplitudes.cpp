/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2026  Laurent Forthomme
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

#include <cmath>
#include <map>
#include <tuple>
#include <vector>

#include "CepGen/Core/Exception.h"
#include "CepGen/HelicityAmplitudes/HelicityAmplitudes.h"

using namespace cepgen;
using namespace cepgen::amplitudes;

/// Model-independent partial-wave helicity amplitudes for \f$\gamma^*\gamma^*\to HX\f$, with \f$H\f$ a spin-0
/// (parity-even) Higgs boson and \f$X\f$ a neutral, C-even particle of integer spin \f$s\f$ and intrinsic parity
/// \f$\eta_X\f$. Only transverse photons are considered.
///
/// Conventions (as in the other models: photon 1 along +z, photon 2 along -z, scattering plane azimuth set to 0):
///  - \f$\Lambda=\lambda_1-\lambda_2\f$ is the spin projection on the beam axis; \f$\lambda_X\f$ the helicity of X,
///    so that the projection on the H direction is \f$-\lambda_X\f$
///  - \f$A_{\lambda_1\lambda_2;\lambda_X}(\theta)=\sum_J(2J+1)a^J_{\lambda_1\lambda_2;\lambda_X}\,d^J_{\Lambda,-\lambda_X}(\theta)\f$
///  - parity: \f$a^J_{-\lambda_1,-\lambda_2;-\lambda_X}=\eta_X(-1)^s a^J_{\lambda_1\lambda_2;\lambda_X}\f$ (J-independent,
///    up to a convention-dependent phase), imposed automatically on the user-provided amplitudes
///  - the contact (non-resonant) amplitudes are returned as \f$M=g\,\hat s\,F(\hat s)\,A\f$, with \f$g=c/\Lambda^2\f$
///    (`coupling`, GeV^-2) and an optional dipole form factor \f$F=(1+\hat s/\Lambda_{FF}^2)^{-2}\f$
///  - s-channel resonances \f$R\f$ (spin \f$J\f$, parity \f$\eta_R\f$, mass \f$M_R\f$, width \f$\Gamma_R\f$) are added coherently:
///    \f$a^J_{\lambda_1\lambda_2;\lambda_X}=c_{\lambda_1\lambda_2}\,b_{\lambda_X}\,M_R^2/(\hat s-M_R^2+iM_R\Gamma_R)\f$ where
///    \f$c\f$ are the \f$\gamma\gamma\to R\f$ and \f$b\f$ the \f$R\to HX\f$ helicity amplitudes (dimensionless).
///    Parity gives \f$c_{--}=\eta_R(-1)^Jc_{++}\f$, \f$c_{-+}=\eta_R(-1)^Jc_{+-}\f$ and
///    \f$b_{-\lambda_X}=\eta_R\eta_X(-1)^{J-s}b_{\lambda_X}\f$ (so \f$b_0\f$ requires \f$\eta_R\eta_X(-1)^{J-s}=+1\f$)
class HXAmplitudes final : public HelicityAmplitudes {
public:
  explicit HXAmplitudes(const ParametersList& params)
      : HelicityAmplitudes(params),
        spin_(steer<int>("spin")),
        eta_(steer<int>("parity")),
        jmax_(steer<int>("jmax")),
        landau_yang_(steer<bool>("landauYang")),
        coupling_(steer<double>("coupling")),
        ff_scale2_(std::pow(steer<double>("formFactorScale"), 2)) {
    if (spin_ < 0)
      throw CG_FATAL("HXAmplitudes") << "X must have a non-negative integer spin (single production of a half-integer "
                                     << "spin particle from γγ → H X would violate fermion number conservation).";
    if (eta_ != 1 && eta_ != -1)
      throw CG_FATAL("HXAmplitudes") << "The intrinsic parity of X must be ±1, got " << eta_ << ".";
    if (jmax_ < 0)
      throw CG_FATAL("HXAmplitudes") << "The maximal partial wave must be positive, got " << jmax_ << ".";
    for (const auto& amplitude : steer<std::vector<ParametersList> >("amplitudes")) {
      const auto l1 = amplitude.get<int>("l1"), l2 = amplitude.get<int>("l2");
      if (std::abs(l1) != 1 || std::abs(l2) != 1)
        throw CG_FATAL("HXAmplitudes") << "Only transverse photons are supported (l1, l2 = ±1).";
      const auto lx = amplitude.get<int>("lX");
      if (std::abs(lx) > spin_)
        throw CG_FATAL("HXAmplitudes") << "|λ_X| = " << std::abs(lx) << " exceeds the spin of X (" << spin_ << ").";
      const auto j = amplitude.get<int>("J");
      if (j < 0 || j > jmax_)
        throw CG_FATAL("HXAmplitudes") << "Partial wave J=" << j << " out of range [0," << jmax_ << "].";
      if (j < std::abs(l1 - l2) || j < std::abs(lx))
        throw CG_FATAL("HXAmplitudes") << "J=" << j << " < max(|Λ|,|λ_X|) is not allowed.";
      if (landau_yang_ && l1 == l2 && j % 2 == 1)
        throw CG_FATAL("HXAmplitudes") << "Landau-Yang: odd J is forbidden for λ₁=λ₂.";
      amplitudes_[{j, +l1, +l2, +lx}] = amplitude.get<std::complex<double> >("amplitude", 0.);
      amplitudes_[{j, -l1, -l2, -lx}] =
          amplitudes_.at({j, +l1, +l2, +lx}) * (eta_ * (spin_ % 2 ? -1. : 1.));  // parity partner
    }
    for (const auto& r : steer<std::vector<ParametersList> >("resonances"))
      addResonance(r);
  }

  static ParametersDescription description() {
    auto desc = HelicityAmplitudes::description();
    desc.setDescription("Model-independent γγ → H X partial-wave helicity amplitudes");
    desc.add("spin", 0).setDescription("integer spin of X");
    desc.add("parity", +1).setDescription("intrinsic parity of X (±1); natural parity if equal to (-1)^spin");
    desc.add("jmax", 0).setDescription("highest partial wave included");
    desc.add("landauYang", false).setDescription("forbid odd J for λ₁=λ₂ (strictly valid for on-shell photons only)");
    desc.add("coupling", 1.e-3)
        .setDescription("contact terms overall scale c/Λ² multiplying ŝ in the amplitude (GeV⁻²)");
    desc.add("formFactorScale", 0.).setDescription("dipole cut-off on √ŝ (GeV) for the contact terms; <=0 disables it");
    desc.add("resonances", std::vector<ParametersList>{})
        .setDescription(
            "s-channel resonances {J, parity, mass, width, c0, c2, b=[{lX, amplitude}, ...]}: "
            "c0/c2 are the γγ → R amplitudes with Λ=0/2, b the R → H X ones (lX >= 0, partners implicit)");
    desc.add("amplitudes", std::vector<ParametersList>{})
        .setDescription("partial-wave amplitudes {J, l1, l2, lX, amplitude (complex)}; parity partners are implicit");
    return desc;
  }

  bool hasAAHX() const override { return true; }

  std::complex<double> amplitude(const Kinematics& kin, const HelicityState& hel) const override {
    const int l1 = hel[0], l2 = hel[1], lh = hel[2], lx = hel[3];
    if (std::abs(l1) != 1 || std::abs(l2) != 1 || lh != 0 || std::abs(lx) > spin_)  // transverse photons, spin-0 H
      return 0.;
    const auto theta = std::acos(kin.cos_theta);
    const int lambda = l1 - l2;
    std::complex<double> out = 0.;
    {  // contact (non-resonant) partial waves
      std::complex<double> contact = 0.;
      for (int j = 0; j <= jmax_; ++j)
        if (const auto it = amplitudes_.find({j, l1, l2, lx}); it != amplitudes_.end())
          contact += (2. * j + 1.) * it->second * wignerD(j, lambda, -lx, theta);
      if (contact != 0.) {
        const auto form_factor = ff_scale2_ > 0. ? 1. / std::pow(1. + kin.shat / ff_scale2_, 2) : 1.;
        out += coupling_ * kin.shat * form_factor * contact;
      }
    }
    for (const auto& res : resonances_) {  // s-channel resonances (coherent sum)
      const auto c = res.production(l1, l2), b = res.decay(lx);
      if (c == 0. || b == 0.)
        continue;
      const auto bw = res.mass * res.mass / std::complex<double>(kin.shat - res.mass * res.mass, res.mass * res.width);
      out += (2. * res.j + 1.) * c * b * bw * wignerD(res.j, lambda, -lx, theta);
    }
    return out;
  }

private:
  /// s-channel resonance R → H X, coupled to γγ
  struct Resonance : SteeredObject<Resonance> {
    explicit Resonance(const ParametersList& params)
        : SteeredObject(params),
          j(steer<int>("J")),
          eta(steer<int>("parity")),
          mass(steer<double>("mass")),
          width(steer<double>("width")),
          c0(steer<std::complex<double> >("c0")),
          c2(steer<std::complex<double> >("c2")) {
      if (j < 0)
        throw CG_FATAL("HXAmplitudes") << "Resonance spin must be a non-negative integer.";
      if (eta != 1 && eta != -1)
        throw CG_FATAL("HXAmplitudes") << "Resonance parity must be ±1.";
      if (mass <= 0. || width <= 0.)
        throw CG_FATAL("HXAmplitudes") << "Resonance mass&width must be positive, got " << mass << ", " << width << ".";
      if (c2 != 0. && j < 2)
        throw CG_FATAL("HXAmplitudes") << "Λ=2 production (c2) requires a resonance with J >= 2.";
      if (c0 == 0. && c2 == 0.)
        throw CG_FATAL("HXAmplitudes") << "Resonance has no γγ production coupling (c0, c2).";
    }

    int j{0}, eta{1};
    double mass{0.}, width{0.};
    std::complex<double> c0{0., 0.}, c2{0., 0.};  ///< γγ → R with Λ=0 (λ₁=λ₂=+1) and Λ=2 (λ₁=+1, λ₂=-1)
    std::map<int, std::complex<double> > b;       ///< R → H X, for λ_X >= 0
    /// γγ → R helicity amplitude (parity partners from c_{-λ₁,-λ₂} = η_R (-1)^J c_{λ₁λ₂})
    std::complex<double> production(int l1, int l2) const {
      const auto base = (l1 == l2) ? c0 : c2;
      return l1 > 0 ? base : (eta * (j % 2 ? -1. : 1.)) * base;
    }
    /// R → H X helicity amplitude (parity partner from b_{-λ} = η_R η_X (-1)^{J-s} b_λ)
    std::complex<double> decay(int lx) const {
      const auto it = b.find(std::abs(lx));
      if (it == b.end())
        return 0.;
      return lx >= 0 ? it->second : it->second * partner_sign;
    }
    double partner_sign{1.};
  };

  void addResonance(const ParametersList& r) {
    Resonance res(r);
    if (landau_yang_ && res.j % 2 == 1 && res.c0 != 0.)
      throw CG_FATAL("HXAmplitudes") << "Landau-Yang: odd-J resonance cannot have Λ=0 production (c0).";
    res.partner_sign = res.eta * eta_ * (((res.j - spin_) % 2 == 0) ? 1. : -1.);
    for (const auto& bb : r.get<std::vector<ParametersList> >("b", {})) {
      const auto lx = bb.get<int>("lX");
      const auto val = bb.get<std::complex<double> >("amplitude", 0.);
      if (lx < 0 || lx > std::min(spin_, res.j))
        throw CG_FATAL("HXAmplitudes") << "Resonance decay helicity lX=" << lx << " must be in [0, min(s,J)].";
      if (lx == 0 && val != 0. && res.partner_sign < 0.)
        throw CG_FATAL("HXAmplitudes") << "R → H X with λ_X=0 vanishes by parity unless η_R η_X (-1)^(J-s) = +1.";
      res.b[lx] = val;
    }
    if (res.b.empty())
      throw CG_FATAL("HXAmplitudes") << "Resonance has no decay coupling b.";
    resonances_.push_back(res);
  }

  /// Wigner small-d function \f$d^j_{m'm}(\beta)\f$ for integer arguments
  static double wignerD(int j, int mp, int m, double beta) {
    if (std::abs(mp) > j || std::abs(m) > j)
      return 0.;
    const auto fact = [](int n) { return std::tgamma(n + 1.); };
    const auto c = std::cos(0.5 * beta), s = std::sin(0.5 * beta);
    const auto pref = std::sqrt(fact(j + mp) * fact(j - mp) * fact(j + m) * fact(j - m));
    double sum = 0.;
    for (int k = 0; k <= 2 * j + 1; ++k) {
      const int a = j + m - k, b = mp - m + k, d = j - mp - k;
      if (a < 0 || b < 0 || d < 0)
        continue;
      sum += (((mp - m + k) % 2 == 0) ? 1. : -1.) * std::pow(c, 2 * j + m - mp - 2 * k) * std::pow(s, mp - m + 2 * k) /
             (fact(a) * fact(k) * fact(b) * fact(d));
    }
    return pref * sum;
  }

  const int spin_;
  const int eta_;
  const int jmax_;
  const bool landau_yang_;
  const double coupling_;
  const double ff_scale2_;

  std::map<std::tuple<int, int, int, int>, std::complex<double> > amplitudes_;
  std::vector<Resonance> resonances_;
};
REGISTER_HELICITY_AMPLITUDES("hxModel", HXAmplitudes);
