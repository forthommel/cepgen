/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2022-2023  Laurent Forthomme
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

#include <plotter.h>

#include "CepGen/Core/Exception.h"
#include "CepGen/Modules/DrawerFactory.h"
#include "CepGen/Utils/Drawable.h"
#include "CepGen/Utils/Drawer.h"

namespace cepgen {
  namespace utils {
    class DrawerLibPlot : public Drawer {
    public:
      explicit DrawerLibPlot(const ParametersList&);

      const DrawerLibPlot& draw(const Graph1D&, const Mode&) const override;
      const DrawerLibPlot& draw(const Graph2D&, const Mode&) const override;
      const DrawerLibPlot& draw(const Hist1D&, const Mode&) const override;
      const DrawerLibPlot& draw(const Hist2D&, const Mode&) const override;

      const DrawerLibPlot& draw(const DrawableColl&,
                                const std::string& name = "",
                                const std::string& title = "",
                                const Mode& mode = Mode::none) const override;
    };

    DrawerLibPlot::DrawerLibPlot(const ParametersList& params) : Drawer(params) {}

    const DrawerLibPlot& DrawerLibPlot::draw(const Graph1D& graph, const Mode&) const {
      CG_WARNING("DrawerLibPlot:draw") << "Not yet implemented.";
      return *this;
    }

    const DrawerLibPlot& DrawerLibPlot::draw(const Graph2D& graph, const Mode&) const {
      CG_WARNING("DrawerLibPlot:draw") << "Not yet implemented.";
      return *this;
    }

    const DrawerLibPlot& DrawerLibPlot::draw(const Hist1D& hist, const Mode&) const {
      CG_WARNING("DrawerLibPlot:draw") << "Not yet implemented.";
      return *this;
    }

    const DrawerLibPlot& DrawerLibPlot::draw(const Hist2D& hist, const Mode&) const {
      CG_WARNING("DrawerLibPlot:draw") << "Not yet implemented.";
      return *this;
    }

    const DrawerLibPlot& DrawerLibPlot::draw(const DrawableColl& objs,
                                             const std::string& name,
                                             const std::string& title,
                                             const Mode& mode) const {
      /*TMultiGraph mg;
      THStack hs;
      ROOTCanvas canv(name.c_str(), "");
      size_t i = 0;
      for (const auto* obj : objs) {
        if (obj->isHist1D()) {
          auto* hist = new TH1D(convert(*dynamic_cast<const Hist1D*>(obj)));
          hist->SetLineColor(ROOTCanvas::colours.at(i++));
          hs.Add(hist);
          canv.AddLegendEntry(hist, obj->title(), "l");
        } else if (obj->isGraph1D()) {
          auto* gr = new TGraphErrors(convert(*dynamic_cast<const Graph1D*>(obj)));
          gr->SetLineColor(ROOTCanvas::colours.at(i++));
          mg.Add(gr);
          canv.AddLegendEntry(gr, obj->title(), "l");
        } else {
          CG_WARNING("DrawerROOT:draw") << "Cannot add drawable '" << obj->name() << "' to the stack.";
          continue;
        }
      }
      const bool has_hists = hs.GetHists() && !hs.GetHists()->IsEmpty();
      const bool has_graphs = mg.GetListOfGraphs() && !mg.GetListOfGraphs()->IsEmpty();
      if (has_hists)
        hs.Draw(mode & Mode::nostack ? "nostack" : "");
      if (has_graphs)
        mg.Draw((std::string("l") + (!has_hists ? "a" : "")).c_str());
      if (has_hists)
        canv.Prettify(hs.GetHistogram());
      else if (has_graphs)
        canv.Prettify(mg.GetHistogram());
      canv.Save("pdf");*/
      CG_WARNING("DrawerLibPlot:draw") << "Not yet implemented.";
      return *this;
    }
  }  // namespace utils
}  // namespace cepgen

REGISTER_DRAWER("libplot", DrawerLibPlot)
