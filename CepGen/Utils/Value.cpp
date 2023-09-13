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

#include <cmath>
#include <iostream>

#include "CepGen/Utils/Value.h"

namespace cepgen {
  template <typename T>
  Value::Value(const T& val, const T& unc) : val_(val), unc_(unc) {}

  T Value::relativeUncertainty() const { return val_ == 0. ? 0. : unc_ / val_; }

  Value Value::operator+(const Value& oth) const { return Value{val_ + oth.val_, std::hypot(unc_, oth.unc_)}; }

  Value Value::operator-(const Value& oth) const { return Value{val_ - oth.val_, std::hypot(unc_, oth.unc_)}; }

  Value Value::operator*(const Value& oth) const {
    const auto prod = val_ * oth.val_;
    return Value{prod, prod * std::hypot(relativeUncertainty(), oth.relativeUncertainty())};
  }

  Value Value::operator/(const Value& oth) const {
    const auto ratio = val_ / oth.val_;
    return Value{ratio, ratio * std::hypot(relativeUncertainty(), oth.relativeUncertainty())};
  }

  std::ostream& operator<<(std::ostream& os, const Value& value) { return os << value.val_ << " +/- " << value.unc_; }
}  // namespace cepgen
