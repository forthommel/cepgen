/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2023  Laurent Forthomme
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

#ifndef CepGen_Utils_Value_h
#define CepGen_Utils_Value_h

#include <iosfwd>

namespace cepgen {
  /// A scalar value with its uncertainty
  template <typename T = double>
  class Value {
  public:
    explicit Value(const T& val = 0., const T& unc = 0.);

    friend std::ostream& operator<<(std::ostream&, const Value&);

    /// Central value extraction
    operator T() const { return val_; }
    /// Absolute uncertainty around the central value
    T uncertainty() const { return unc_; }
    /// Relative uncertainty around the central value
    T relativeUncertainty() const;

    /// Comparison operator
    bool operator<(const Value& oth) const { return val_ < oth.val_; }

    //--- error propagation operators

    Value operator+(const Value&) const;
    template <typename U>
    inline Value operator+(const U& cst) const {
      return Value{val_ + cst, unc_};
    }
    Value operator-(const Value&) const;
    template <typename U>
    inline Value operator-(const U& cst) const {
      return Value{val_ - cst, unc_};
    }
    Value operator*(const Value&) const;
    template <typename U>
    inline Value operator*(const U& cst) const {
      return Value{val_ * cst, unc_ * cst};
    }
    Value operator/(const Value&) const;
    template <typename U>
    inline Value operator/(const U& cst) const {
      return Value{val_ / cst, unc_ / cst};
    }

  private:
    T val_;  ///< Central value
    T unc_;  ///< Uncertainty on value
  };
}  // namespace cepgen

#endif
