/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2026 Laurent Forthomme
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include <complex>
#include <thread>

#include "CepGen/Utils/ArgumentsParser.h"
#include "CepGen/Utils/Test.h"
#include "CepGenMadGraph/OffShell.h"

using namespace cepgen::mg5amc;

int main(int argc, char* argv[]) {
  cepgen::ArgumentsParser(argc, argv).parse();

  const double q1[4] = {10.2, 1.5, -0.7, 10.0}, q2[4] = {9.8, -0.4, 2.0, -9.5};
  const std::array<double, 4> eps1{0., 0.9, -0.42, 0.}, eps2{0., -0.2, 0.98, 0.};
  std::complex<double> vc[6];

  CG_TEST_EQUAL(offshell::fill(q1, -1, vc), false, "inactive: stock wavefunction");

  offshell::set(0, q1, eps1);
  offshell::set(1, q2, eps2);
  CG_TEST_EQUAL(offshell::fill(q1, -1, vc), true, "leg 0 override");
  CG_TEST_EQUIV(vc[0].real(), -q1[0], "leg 0 momentum (energy)");
  CG_TEST_EQUIV(vc[0].imag(), -q1[3], "leg 0 momentum (pz)");
  CG_TEST_EQUIV(vc[1].real(), -q1[1], "leg 0 momentum (px)");
  CG_TEST_EQUIV(vc[1].imag(), -q1[2], "leg 0 momentum (py)");
  for (size_t mu = 0; mu < 4; ++mu)
    CG_TEST_EQUIV(vc[2 + mu].real(), eps1.at(mu), "leg 0 polarisation component " + std::to_string(mu));
  CG_TEST_EQUAL(offshell::fill(q2, -1, vc), true, "leg 1 override");
  CG_TEST_EQUIV(vc[3].real(), eps2.at(1), "leg 1 polarisation (x)");
  CG_TEST_EQUAL(offshell::hits(0), (size_t)1, "leg 0 hits");
  CG_TEST_EQUAL(offshell::hits(1), (size_t)1, "leg 1 hits");

  const double other[4] = {5., 1., 1., 4.};
  CG_TEST_EQUAL(offshell::fill(other, -1, vc), false, "non-matching momentum falls through");
  CG_TEST_EQUAL(offshell::fill(q1, +1, vc), false, "outgoing leg falls through");

  bool other_thread_active = true;
  std::thread([&] { other_thread_active = offshell::fill(q1, -1, vc); }).join();
  CG_TEST_EQUAL(other_thread_active, false, "state is thread-local");

  offshell::reset();
  CG_TEST_EQUAL(offshell::fill(q1, -1, vc), false, "reset: stock wavefunction");

  CG_TEST_SUMMARY;
}
