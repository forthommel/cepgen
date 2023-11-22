#include "CepGen/Core/Exception.h"
#include "CepGen/Generator.h"
#include "CepGen/Modules/ProcessFactory.h"
#include "CepGen/Physics/PDG.h"
#include "CepGen/Process/Process.h"
#include "CepGen/Utils/ArgumentsParser.h"
#include "CepGenAddOns/MadGraphWrapper/MadGraphProcess.h"

int main(int argc, char* argv[]) {
  cepgen::ArgumentsParser(argc, argv).parse();
  cepgen::initialise();

  //1000013
  cepgen::PDG::get().define(cepgen::ParticleProperties(1000013, "mul", "slepton-2", 0., 100., 0., 3, false));

  auto mg5 = cepgen::ProcessFactory::get().build(
      "mg5_aMC",
      cepgen::ParametersList().set<std::string>("model", "MSSM_SLHA2").set<std::string>("process", "a a > mul+ mul-"));

  return 0;
}
