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
        int id = num_hist_++, nbins = hist.nbins();
        double xmin = hist.range().min(), xmax = hist.range().max();
        auto hist_name = hist.name();
        //xhinit_(id, xmin, xmax, nbins, &(hist_name[0]));
        auto* title = const_cast<char*>(" ");
        xhinit_(id, xmin, xmax, nbins, title);
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
