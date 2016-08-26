#include <iostream>

#include "include/MCGen.h"
#include "export/EventWriter.h"

using namespace std;

/**
 * Main caller for this Monte Carlo generator. Loads the configuration files'
 * variables if set as an argument to this program, else loads a default
 * "LHC-like" configuration, then launches the cross-section computation and
 * the events generation.
 * \author Laurent Forthomme <laurent.forthomme@cern.ch>
 */
int main(int argc, char* argv[]) {
  MCGen mg;
  
  if (argc==1) InError("No config file provided.");

  Debugging(Form("Reading config file stored in %s", argv[1]));
  if (!mg.parameters->ReadConfigFile(argv[1])) {
    Information(Form("Error reading the configuration!\n\t"
                     "Please check your input file (%s)", argv[1]));
    return -1;
  }

  // We might want to cross-check visually the validity of our run
  mg.parameters->Dump();

  // Let there be cross-section...
  double xsec, err;
  mg.ComputeXsection(&xsec, &err);

  HepMC::GenCrossSection xs;
  xs.set_cross_section(xsec, err);

  if (!mg.parameters->generation) return 0;

  EventWriter::HepMC::output output("example.dat", std::ios::out);

  // The events generation starts here !
  for (int i=0; i<mg.parameters->maxgen; i++) {
    if (i%10000==0)
      cout << "Generating event #" << i+1 << endl;
    const Event ev = *mg.GenerateOneEvent();
    HepMC::GenEvent* hev = EventWriter::HepMC::Event(ev);
    hev->set_cross_section(xs);                                                                                                                                                                                                                                      
    hev->set_event_number(i); 

    output << (hev);
  }

  //mg.parameters->StoreConfigFile("lastrun.card");

  return 0;
}
