#ifndef CepGenAddOns_Herwig6Wrapper_HerwigInterface_h
#define CepGenAddOns_Herwig6Wrapper_HerwigInterface_h

namespace {
  extern "C" {
  void hwigin_();  // initialise other common blocks
  void hweini_();  // initialise elementary process
  void hwudat_();
  void hwefin_();
  void hwepro_();
  void hwuine_();  // initialise event
  void hwbgen_();  // generate parton cascades
  void hwcfor_();  // perform cluster formation
  void hwcdec_();  // perform cluster decay
  void hwdhad_();  // perform unstable particles decay
  void hwufne_();  // finalise event
  void hwuepr_();  // print event content
  void hwuidt_(int& iopt, int& ipdg, int& iwig, char nwig[8]);
  void hwuinc_();  // compute constants and lookup tables
  void hwhrem_(int&, int&);
  void hwhsct_(int& report, int& firstc, int& jmueo, double& ptjim);
  void upinit_() {}
  void upevnt_() {}
  void hwaend_() { CG_INFO("Herwig6Hadroniser") << "End of run"; }
  void timel_(double& tres) { tres = 1.e10; }
  extern struct { int ipart1, ipart2; } hwbeam_;
  extern struct { char part1[8], part2[8]; } hwbmch_;
  /// spin correlations in event
  extern struct {
    int ndecsy, nsearch, lrdec, lwdec, syspin, threeb, fourb;
    char taudec[6];
  } hwdspn_;
  extern struct {
    char autpdf[20][2];
    char bdecay[4];
  } hwprch_;
  extern struct {
    double ebeam1, ebeam2, pbeam1, pbeam2;
    int iproc, maxev;
  } hwproc_;
  const int np = 4000;
  /// Particles content of the event
  extern struct {
    /// Event number
    int nevhep;
    /// Number of particles in the event
    int nhep;
    /// Particles' status code
    int isthep[np];
    /// Particles' PDG id
    int idhep[np];
    /// Particles' mothers
    int jmohep[np][2];
    /// Particles' daughters
    int jdahep[np][2];
    /// Particles' kinematics, in GeV (px, py, pz, E, M)
    double phep[np][5];
    /// Primary vertex for the particles
    double vhep[np][5];
  } hepevt_;
  const int maxpup = 100;
  extern struct {
    int idbmup[2];
    double ebmup[2];
    int pdfgup[2], pdfsup[2], idwtup, nprup;
    double xsecup[maxpup], xerrup[maxpup], xmaxup[maxpup];
    int lprup[maxpup];
  } heprup_;

  extern struct {
    double asfixd, clq[6][7], coss, costh, ctmax, disf[2][13];
    double emlst, emmax, emmin, empow, emsca, epoln[3];
    double gcoef[7], gpoln, omega0, phomas, ppoln[3];
    double ptmax, ptmin, ptpow;
    double q2max, q2min, q2pow;
    double q2wwmn, q2wwmx;
    double qlim, sins, thmax, y4jt, tmnisr, tqwt;
    double xx[2], xlmin, xxmin;
    double ybmax, ybmin, yjmax, yjmin;
    double ywwmax, ywwmin;
    double whmin, zjmax, zmxisr;
    int iaphig, ibrn[2], ibsh, ico[10], idcmf, idn[10], iflmax, iflmin, ihpro, ipro;
    int mapq, maxfl, bgshat;
    int colisr, fstevt, fstwgt, genev, hvfcen, tpol, durham;  // logical
  } hwhard_;
  extern struct {
    double avwgt, evwgt, gamwt, tlout;
    double wbigst, wgtmax, wgtsum, wsqsum;
    int idhw[np], ierror, istat;
    int lwevt, maxer, maxpr, nowgt, nrn[2];
    int numer, numeru, nwgts, gensof;
  } hwevnt_;

  /// Basic parameters (and quantities derived from them)
  extern struct {
    double afch[2][16], alphem, b1l1m, betaf, btclm, cafac, cffac, clmax, clpow;
    double clsmr[2], cspeed, ensof, etamix, f0mix, f1mix, f2mix, gamh, gamw, gamz, gamzp;
    double gev2nb, h1mix, pdiqk, pgsmx, pgspl[4], phimx, pifac, prsof, pslpt[2], ptrms, pxrms;
    double qcdl3, qcdl5, qcdlam, qdiqk, qfch[16], qg, qspac, qv;
    double scabi, swein, tmtop, vfch[2][16], vckm[3][3], vgcut, vqcut, vpcut;
    double zbinm, effmin, omhmix, et2mix, ph3mix, gcutme;
    int ioprem, iprint, ispac;
    int lrsud, lwsud, modpdf[2], nbtry, ncolo, nctry, ndtry, netry, nvlav, ngspl, nstru, nstry, nzbin;
    int iop4jt[2], nprfmt, azsoft, azspin;
    int cldir[2], hardme, nospac, prndec, prvtx, softme, zprime, prndef, prntex, prnweb;
  } hwpram_;
  }
}  // namespace

#endif
