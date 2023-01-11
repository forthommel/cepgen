#include "CepGen/Core/GeneratorWorker.h"
#include "CepGen/Generator.h"
#include "CepGen/Integration/ProcessIntegrand.h"
#include "CepGen/Modules/ProcessFactory.h"
#include "CepGen/Parameters.h"
#include "CepGen/Process/Process.h"
#include "CepGen/Utils/ArgumentsParser.h"
#include "CepGen/Utils/Test.h"

using namespace std;

int main(int argc, char* argv[]) {
  int num_threads;

  cepgen::ArgumentsParser(argc, argv).addOptionalArgument("num-threads,n", "number of threads", &num_threads, 4).parse();

  cepgen::Generator gen;

  const cepgen::Kinematics kin(cepgen::ParametersList()
                                   .set<int>("mode", 1)
                                   .set<vector<int> >("pdgIds", {2212, 2212})
                                   .set<vector<double> >("pz", {6.8e3, 6.8e3})
                                   .set<double>("ptmin", 3.)
                                   .set<double>("mxmin", 10.));

  gen.parametersRef().setProcess(cepgen::proc::ProcessFactory::get().build("lpair"));
  gen.parametersRef().process().setKinematics(kin);

  CG_TEST(gen.parametersRef().process().kinematics() == kin, "Process kinematics");

  gen.parametersRef().generation().setNumThreads(num_threads);
  CG_TEST_EQUAL((int)gen.parametersRef().generation().numThreads(), num_threads, "Number of threads");

  {  // compute the cross section to initialise the workers
    double xsec, xsec_unc;
    gen.computeXsection(xsec, xsec_unc);
  }

  for (int i = 0; i < num_threads; ++i) {
    CG_TEST(&gen.worker(i).integrand().process() != &gen.parametersRef().process(),
            "Thread #" + to_string(i) + " process");
    CG_TEST(gen.worker(i).integrand().process().kinematics() == kin, "Thread #" + to_string(i) + " process kinematics");
  }

  CG_TEST_SUMMARY;
}
