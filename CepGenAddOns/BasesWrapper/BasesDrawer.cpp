/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2024  Laurent Forthomme
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

#include <cstring>

#include "CepGen/Modules/DrawerFactory.h"
#include "CepGen/Utils/Drawer.h"
#include "CepGen/Utils/Graph.h"
#include "CepGen/Utils/Histogram.h"
#include "CepGen/Utils/Message.h"
#include "CepGenAddOns/BasesWrapper/BasesCommonBlocks.h"

namespace cepgen {
  namespace utils {
    class BasesDrawer : public Drawer {
    public:
      explicit BasesDrawer(const ParametersList& params) : Drawer(params) {}

      static ParametersDescription description() {
        auto desc = Drawer::description();
        desc.setDescription("Bases histogram drawer");
        return desc;
      }

      const BasesDrawer& draw(const Graph1D&, const Mode&) const override { return *this; }
      const BasesDrawer& draw(const Graph2D&, const Mode&) const override { return *this; }
      const BasesDrawer& draw(const Hist1D& hist, const Mode&) const override {
        int id = ++num_hist_, nbins = hist.nbins();
        double xmin = hist.range().min(), xmax = hist.range().max();
        int hist_name_len = hist.name().size() - 1;
        std::vector<char> hist_name(hist_name_len + 1);
        std::strcpy(hist_name.data(), hist.name().c_str());
        xhinit_(id, xmin, xmax, nbins, hist_name.data(), hist_name_len);
        for (int i = 0; i < nbins; ++i) {
          auto dx = hist.binRange(i).x(0.5), fx = (double)hist.value(i);
          xhfill_(id, dx, fx);
        }
        int lu = 6, ifg = 1;
        xhplot_(lu, ifg, id);
        return *this;
      }
      const BasesDrawer& draw(const Hist2D&, const Mode&) const override { return *this; }

      const BasesDrawer& draw(const DrawableColl&,
                              const std::string& name = "",
                              const std::string& title = "",
                              const Mode& mode = Mode::none) const override {
        return *this;
      }

    private:
      mutable int num_hist_{0};
    };
  }  // namespace utils
}  // namespace cepgen
using cepgen::utils::BasesDrawer;
REGISTER_DRAWER("bases", BasesDrawer);
