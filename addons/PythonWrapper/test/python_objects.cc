/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2022-2025  Laurent Forthomme
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

// clang-format off
#include "CepGenPython/Environment.h"
// clang-format on
#include "CepGen/Utils/ArgumentsParser.h"
#include "CepGen/Utils/Test.h"
#include "CepGenPython/ObjectPtr.h"

using namespace std;
using namespace std::string_literals;

#define TEST_TYPE(type, object)                            \
  {                                                        \
    auto py_obj = cepgen::python::ObjectPtr::make(object); \
    auto ret = py_obj.value<type>();                       \
    CG_TEST_EQUAL(ret, object, std::string(#object));      \
  }

int main(int argc, char* argv[]) {
  cepgen::ArgumentsParser(argc, argv).parse();

  cepgen::python::Environment env(cepgen::ParametersList{});
  TEST_TYPE(bool, true)
  TEST_TYPE(bool, false)
  TEST_TYPE(std::complex<double>, std::complex<double>(1.5, -2.5))
  TEST_TYPE(string, "Héhéhé, test @ ünıc0d€ 🐗"s)
  TEST_TYPE(cepgen::Limits, cepgen::Limits(-2., 3.1))
  TEST_TYPE(cepgen::ParametersList,
            cepgen::ParametersList()
                .set("foo", 42)
                .set("bar", M_PI)
                .set("baz", "héhé"s)
                .set("flag", true)
                .set("coupling", std::complex<double>(0.5, 1.25))
                .set("plist", cepgen::ParametersList().set("foo", 10).set("bar", 42.42).set("baz", "hîhî"s)))

  CG_TEST_SUMMARY;
}
