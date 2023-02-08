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

#ifndef CepGenAddOns_MartyWrapper_MartyModelFactory_h
#define CepGenAddOns_MartyWrapper_MartyModelFactory_h

#include "CepGen/Modules/ModuleFactory.h"

/// Add a Marty model definition to the factory
#define REGISTER_MARTY_MODEL(name, obj)                                       \
  namespace cepgen {                                                          \
    struct BUILDERNM(obj) {                                                   \
      BUILDERNM(obj)() { MartyModelFactory::get().registerModel<obj>(name); } \
    };                                                                        \
    static const BUILDERNM(obj) gMartyModel##obj;                             \
  }

namespace mty {
  class Model;
}

namespace cepgen {
  class MartyModelFactory {
  public:
    static MartyModelFactory& get() {
      static MartyModelFactory factory;
      return factory;
    }
    std::vector<std::string> models() const {
      std::vector<std::string> models;
      for (const auto& mod : map_)
        models.emplace_back(mod.first);
      return models;
    }
    /// Register a named model in the database
    /// \tparam T Class to register (inherited from mty::Model base class)
    template <typename T>
    void registerModel(const std::string& name) {
      static_assert(std::is_base_of<mty::Model, T>::value,
                    "\n\n  *** Failed to register an object with improper inheritance into the factory. ***\n");
      if (map_.count(name) > 0)
        throw std::invalid_argument("\n\n  *** detected a duplicate Marty model registration for name \"" + name +
                                    "\"! ***\n");
      map_.insert(std::make_pair(name, &buildModel<T>));
    }
    /// Build one instance of a named model
    /// \param[in] name Model name to retrieve
    std::unique_ptr<mty::Model> build(const std::string& name) const {
      if (map_.count(name) == 0)
        throw std::invalid_argument("\n\n  *** failed to retrieve a Marty model with name '" + name + "'! ***\n");
      return map_.at(name)();
    }

  private:
    MartyModelFactory() = default;
    /// Marty model constructor
    template <typename T>
    static std::unique_ptr<mty::Model> buildModel() {
      return std::unique_ptr<mty::Model>(new T);
    }
    typedef std::unique_ptr<mty::Model> (*Builder)();
    std::unordered_map<std::string, Builder> map_;
  };
}  // namespace cepgen

#endif
