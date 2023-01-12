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

#ifndef CepGen_Utils_ThreadSafeObject_h
#define CepGen_Utils_ThreadSafeObject_h

#include <functional>
#include <mutex>

namespace cepgen {
  namespace utils {
    /// Helper base object containing its own mutex and a set of const-qualified methods
    /// to apply with a proper locking state
    /// \date Jan 2023
    template <typename T>
    class ThreadSafeObject {
    public:
      /// Apply a non-const qualified method to the base object
      const ThreadSafeObject<T>& apply(std::function<void(T&)> callback) const {
        mutex_.lock();
        callback((T&)*this);
        mutex_.unlock();
        return *this;
      }
      /// Apply a non-const qualified probing method to the base object and return its result
      template <typename U>
      U probe(std::function<U(T&)> callback) const {
        mutex_.lock();
        auto ret = callback((T&)*this);
        mutex_.unlock();
        return ret;
      }

    private:
      mutable std::mutex mutex_;
    };
  }  // namespace utils
}  // namespace cepgen

#endif
