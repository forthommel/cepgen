/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2013-2021  Laurent Forthomme
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

#include <Math/Integrator.h>
#include <Math/IntegratorMultiDim.h>

#include "CepGen/Core/Exception.h"
#include "CepGen/Integration/Integrand.h"
#include "CepGen/Integration/Integrator.h"
#include "CepGen/Modules/IntegratorFactory.h"

namespace cepgen {
  /// FOAM general-purpose integration algorithm
  class IntegratorROOT final : public Integrator {
  public:
    explicit IntegratorROOT(const ParametersList&);

    static ParametersDescription description();

    void integrate(double&, double&) override;

  private:
    std::function<double(const double*)> func_;
    std::function<double(double)> func_1d_;
    std::vector<double> min_, max_;
    /// integration type (adaptive, MC methods, etc..)
    const std::string type_;
    /// desired absolute Error
    const double absTol_;
    /// desired relative Error
    const double relTol_;
    /// maximum number of sub-intervals
    const unsigned int size_;

    std::unique_ptr<ROOT::Math::IntegratorMultiDim> integr_;
    std::unique_ptr<ROOT::Math::IntegratorOneDim> integr_1d_;
  };

  IntegratorROOT::IntegratorROOT(const ParametersList& params)
      : Integrator(params),
        func_([=](const double* x) -> double {
          return integrand_->eval(std::vector<double>(x, x + integrand_->size()));
        }),
        func_1d_([=](double x) -> double { return integrand_->eval(std::vector<double>{x}); }),
        type_(steer<std::string>("type")),
        absTol_(steer<double>("absTol")),
        relTol_(steer<double>("relTol")),
        size_(steer<int>("size")) {
    {
      auto type = ROOT::Math::IntegratorMultiDim::Type::kDEFAULT;
      if (type_ == "adaptive")
        type = ROOT::Math::IntegratorMultiDim::Type::kADAPTIVE;
      else if (type_ == "plain")
        type = ROOT::Math::IntegratorMultiDim::Type::kPLAIN;
      else if (type_ == "miser")
        type = ROOT::Math::IntegratorMultiDim::Type::kMISER;
      else if (type_ == "vegas")
        type = ROOT::Math::IntegratorMultiDim::Type::kVEGAS;
      integr_.reset(new ROOT::Math::IntegratorMultiDim(type, absTol_, relTol_, size_));
    }
    {
      auto type = ROOT::Math::IntegratorOneDim::Type::kDEFAULT;
      if (type_ == "gauss")
        type = ROOT::Math::IntegratorOneDim::Type::kGAUSS;
      else if (type_ == "legendre")
        type = ROOT::Math::IntegratorOneDim::Type::kLEGENDRE;
      else if (type_ == "adaptive")
        type = ROOT::Math::IntegratorOneDim::Type::kADAPTIVE;
      else if (type_ == "adaptiveSingular")
        type = ROOT::Math::IntegratorOneDim::Type::kADAPTIVESINGULAR;
      else if (type_ == "nonAdaptive")
        type = ROOT::Math::IntegratorOneDim::Type::kNONADAPTIVE;
      integr_1d_.reset(new ROOT::Math::IntegratorOneDim(type, absTol_, relTol_, size_));
    }
    //--- a bit of printout for debugging
    CG_DEBUG("Integrator:build") << "ROOT generic integrator built\n\t"
                                 << "N-dimensional type: " << integr_->Name() << ",\n\t"
                                 << "1-dimensional type: " << integr_1d_->Name() << ",\n\t"
                                 << "Absolute tolerance: " << absTol_ << ",\n\t"
                                 << "Relative tolerance: " << relTol_ << ",\n\t"
                                 << "Number of sub-intervals: " << size_ << ".";
  }

  void IntegratorROOT::integrate(double& result, double& abserr) {
    if (!initialised_) {
      if (integrand_->size() == 1)
        integr_1d_->SetFunction(func_1d_);
      else
        integr_->SetFunction(func_, integrand_->size());
      min_.clear();
      max_.clear();
      for (const auto& lim : limits_) {
        min_.emplace_back(lim.min());
        max_.emplace_back(lim.max());
      }
      initialised_ = true;
    }

    if (integrand_->size() == 1) {
      result_ = result = integr_1d_->Integral(limits_.at(0).min(), limits_.at(0).max());
      err_result_ = abserr = integr_1d_->Error();
    } else {
      result_ = result = integr_->Integral(min_.data(), max_.data());
      err_result_ = abserr = integr_->Error();
    }
  }

  ParametersDescription IntegratorROOT::description() {
    auto desc = Integrator::description();
    desc.setDescription("ROOT general purpose MC integrator");
    desc.add<std::string>("type", "default");
    desc.add<double>("absTol", -1.);
    desc.add<double>("relTol", -1.);
    desc.add<int>("size", 0);
    return desc;
  }
}  // namespace cepgen

REGISTER_INTEGRATOR("ROOT", IntegratorROOT)
