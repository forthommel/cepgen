#include <matplotlibcpp.h>

#include "CepGen/Core/Exception.h"
#include "CepGen/Modules/DrawerFactory.h"
#include "CepGen/Utils/Drawer.h"
#include "CepGen/Utils/Graph.h"
#include "CepGen/Utils/Histogram.h"

namespace plt = matplotlibcpp;

namespace cepgen {
  namespace utils {
    class DrawerMatplotlib : public Drawer {
    public:
      explicit DrawerMatplotlib(const ParametersList&);

      static ParametersDescription description() {
        auto desc = Drawer::description();
        desc.setDescription("Matplotlib plotter");
        desc.add<bool>("tight", true).setDescription("use a compact layout with minimal margins");
        return desc;
      }

      const DrawerMatplotlib& draw(const Graph1D&, const Mode&) const override;
      const DrawerMatplotlib& draw(const Graph2D&, const Mode&) const override;
      const DrawerMatplotlib& draw(const Hist1D&, const Mode&) const override;
      const DrawerMatplotlib& draw(const Hist2D&, const Mode&) const override;

      const DrawerMatplotlib& draw(const DrawableColl&,
                                   const std::string& name = "",
                                   const std::string& title = "",
                                   const Mode& mode = Mode::none) const override;

    private:
      static void plot(const Graph1D&, const Mode&);
      static void plot(const Graph2D&, const Mode&);
      static void plot(const Hist1D&, const Mode&);
      void postDraw(const Drawable&, const Mode&) const;
      const bool tight_;
    };

    DrawerMatplotlib::DrawerMatplotlib(const ParametersList& params) : Drawer(params), tight_(steer<bool>("tight")) {}

    const DrawerMatplotlib& DrawerMatplotlib::draw(const Graph1D& graph, const Mode& mode) const {
      plt::figure();
      plot(graph, mode);
      postDraw(graph, mode);
      plt::title(graph.title());
      plt::save(graph.name() + ".pdf");
      return *this;
    }

    const DrawerMatplotlib& DrawerMatplotlib::draw(const Graph2D& graph, const Mode& mode) const {
      plt::figure();
      plot(graph, mode);
      postDraw(graph, mode);
      plt::title(graph.title());
      plt::save(graph.name() + ".pdf");
      return *this;
    }

    const DrawerMatplotlib& DrawerMatplotlib::draw(const Hist1D& hist, const Mode& mode) const {
      plt::figure();
      plot(hist, mode);
      postDraw(hist, mode);
      plt::title(hist.title());
      plt::save(hist.name() + ".pdf");
      return *this;
    }

    const DrawerMatplotlib& DrawerMatplotlib::draw(const Hist2D&, const Mode&) const {
      CG_WARNING("DrawerMatplotlib:draw") << "Not yet implemented.";
      return *this;
    }

    const DrawerMatplotlib& DrawerMatplotlib::draw(const DrawableColl& objs,
                                                   const std::string& name,
                                                   const std::string& title,
                                                   const Mode& mode) const {
      try {
        plt::figure();
        const Drawable* first_obj = nullptr;
        for (const auto* obj : objs) {
          if (obj->isHist1D()) {
            auto* hist = dynamic_cast<const Hist1D*>(obj);
            plot(*hist, mode);
            if (!first_obj)
              first_obj = hist;
          }
          if (obj->isGraph1D()) {
            auto* gr = dynamic_cast<const Graph1D*>(obj);
            plot(*gr, mode);
            if (!first_obj)
              first_obj = gr;
          }
        }
        if (!title.empty())
          plt::title(title);
        if (first_obj)
          postDraw(*first_obj, mode);
        if (objs.size() > 1)
          plt::legend();
        plt::save(name + ".pdf");
      } catch (const std::runtime_error& err) {
        CG_WARNING("DrawerMatplotlib:draw") << "Failed to draw a plots collection. Matplotlib error: " << err.what();
      }
      return *this;
    }

    void DrawerMatplotlib::plot(const Graph1D& gr, const Mode& mode) {
      std::vector<double> x, y, xerr, yerr;
      for (const auto& pt : gr.points()) {
        x.emplace_back(pt.first.value);
        y.emplace_back(pt.second.value);
        xerr.emplace_back(pt.first.value_unc);
        yerr.emplace_back(pt.second.value_unc);
      }
      if ((mode & Mode::logx) && (mode & Mode::logy))
        plt::named_loglog(gr.title(), x, y);
      else if (mode & Mode::logx)
        plt::named_semilogx(gr.title(), x, y);
      else if (mode & Mode::logy)
        plt::named_semilogy(gr.title(), x, y);
      else if (yerr != std::vector<double>(yerr.size(), 0.))
        plt::errorbar(x, y, yerr, {{"label", gr.title()}, {"linestyle", ""}});
      else
        plt::plot(x, y, {{"label", gr.title()}});
    }

    void DrawerMatplotlib::plot(const Graph2D& gr, const Mode&) {
      std::vector<std::vector<double> > x, y, z;
      for (const auto& xv : gr.points()) {
        const auto xval = xv.first.value;
        std::vector<double> xrow, yrow, zrow;
        for (const auto& yv : xv.second) {
          xrow.emplace_back(xval);
          yrow.emplace_back(yv.first.value);
          zrow.emplace_back(yv.second.value);
        }
        x.emplace_back(xrow);
        y.emplace_back(yrow);
        z.emplace_back(zrow);
      }
      plt::plot_surface(x, y, z, {{"label", gr.title()}});
      //plt::contour(x, y, z, {{"label", gr.title()}});
      plt::set_zlabel(gr.zAxis().label());
    }

    void DrawerMatplotlib::plot(const Hist1D& hist, const Mode& mode) {
      std::vector<double> x, y, yerr;
      for (size_t ibin = 0; ibin < hist.nbins(); ++ibin) {
        x.emplace_back(hist.binRange(ibin).x(0.5));
        y.emplace_back(hist.value(ibin));
        yerr.emplace_back(hist.valueUnc(ibin));
      }
      //plt::bar(x, y, "", "", 1., {{"label", hist.title()}});
      //plt::bar(x, y);
      std::map<std::string, std::string> plot_style = {{"label", hist.title()}, {"drawstyle", "steps"}};
      if ((mode & Mode::logx) && (mode & Mode::logy))
        plt::named_loglog(hist.title(), x, y, "o");
      else if (mode & Mode::logx)
        plt::named_semilogx(hist.title(), x, y, "o");
      else if (mode & Mode::logy)
        plt::named_semilogy(hist.title(), x, y, "o");
      else if (yerr != std::vector<double>(yerr.size(), 0.))
        plt::errorbar(x, y, yerr, plot_style);
      else
        plt::plot(x, y, plot_style);
    }

    void DrawerMatplotlib::postDraw(const Drawable& dr, const Mode& mode) const {
      if (mode & Mode::grid)
        plt::grid(true);
      const auto& yrange = dr.yAxis().range();
      if (yrange.valid()) {
        auto rng = plt::ylim();
        if (yrange.hasMin())
          rng[0] = yrange.min();
        if (yrange.hasMax())
          rng[1] = yrange.max();
        try {
          plt::ylim(rng.at(0), rng.at(1));
        } catch (const std::runtime_error& err) {
          CG_WARNING("DrawerMatplotlib:postDraw")
              << "Failed to set Y range to " << rng << ". Matplotlib error: " << err.what();
        }
      }
      plt::xlabel(dr.xAxis().label());
      plt::ylabel(dr.yAxis().label());
      if (tight_)
        plt::tight_layout();
    }
  }  // namespace utils
}  // namespace cepgen

REGISTER_DRAWER("matplotlib", DrawerMatplotlib)
