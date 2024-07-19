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

#include "CepGen/Cards/Handler.h"
#include "CepGen/Core/RunParameters.h"
#include "CepGen/Modules/CardsHandlerFactory.h"
#include "CepGen/Utils/ArgumentsParser.h"
#include "CepGen/Utils/Test.h"
#include "CepGenAddOns/PythonWrapper/Utils.h"

using namespace std;

int main(int argc, char* argv[]) {
  string path;
  cepgen::ArgumentsParser(argc, argv)
      .addOptionalArgument("path,p", "path to the Python card to steer", &path, "Cards/lpair_cfg.py")
      .parse();

  const auto card = cepgen::CardsHandlerFactory::get().build(".py");
  const auto [mod_path, mod_name] = cepgen::python::pythonPath(path);
  card->parseFile(path);
  CG_LOG << path << ":" << *card->runParameters();

  CG_TEST_SUMMARY;
}
