#ifndef eventCounter_H
#define eventCounter_H

#include <TROOT.h>
#include <iostream>

struct eventCounter_t {
  ULong_t numEvents;
  ULong_t numEventsPassedEvtTrigger;
  ULong_t numDielectronsUnweighted;
  //ULong_t numPassedGoodPV;
  // weighted counts
  double numDielectrons;
  double numDielectronsGoodEta;
  double numDielectronsGoodEt;
  double numDielectronsHLTmatched;
  double numDielectronsIDpassed;
  double numDielectronsGoodMass;
  double scale;

  eventCounter_t() :
    numEvents(0), numEventsPassedEvtTrigger(0),
    numDielectronsUnweighted(0),
    //numPassedGoodPV(0),
    numDielectrons(0), 
    numDielectronsGoodEta(0), numDielectronsGoodEt(0),
    numDielectronsHLTmatched(0),
    numDielectronsIDpassed(0),
    numDielectronsGoodMass(0),
    scale(1.)
  {}

  void clear() {
    numEvents=0; numEventsPassedEvtTrigger=0;
    numDielectronsUnweighted=0;
    numDielectrons=0;
    numDielectronsGoodEta=0; numDielectronsGoodEt=0;
    numDielectronsHLTmatched=0;
    numDielectronsIDpassed=0;
    numDielectronsGoodMass=0;
    scale=1.;
  }

  void numDielectrons_inc() { numDielectrons+=scale; }
  void numDielectronsGoodEta_inc() { numDielectronsGoodEta+=scale; }
  void numDielectronsGoodEt_inc() { numDielectronsGoodEt+=scale; }
  void numDielectronsHLTmatched_inc() { numDielectronsHLTmatched+=scale; }
  void numDielectronsIDpassed_inc() { numDielectronsIDpassed+=scale; }
  void numDielectronsGoodMass_inc() { numDielectronsGoodMass+=scale; }
  void setScale(double sc) { scale=sc; }

  void assign(const eventCounter_t &e) {
    numEvents=e.numEvents;
    numEventsPassedEvtTrigger=e.numEventsPassedEvtTrigger;
    numDielectronsUnweighted=e.numDielectronsUnweighted;
    numDielectrons=e.numDielectrons;
    numDielectronsGoodEta=e.numDielectronsGoodEta;
    numDielectronsGoodEt=e.numDielectronsGoodEt;
    numDielectronsHLTmatched=e.numDielectronsHLTmatched;
    numDielectronsIDpassed=e.numDielectronsIDpassed;
    numDielectronsGoodMass=e.numDielectronsGoodMass;
    scale=e.scale;
  }

  void add(const eventCounter_t &e) {
    numEvents+=e.numEvents;
    numEventsPassedEvtTrigger+=e.numEventsPassedEvtTrigger;
    numDielectronsUnweighted+=e.numDielectronsUnweighted;
    numDielectrons+=e.numDielectrons;
    numDielectronsGoodEta+=e.numDielectronsGoodEta;
    numDielectronsGoodEt+=e.numDielectronsGoodEt;
    numDielectronsHLTmatched+=e.numDielectronsHLTmatched;
    numDielectronsIDpassed+=e.numDielectronsIDpassed;
    numDielectronsGoodMass+=e.numDielectronsGoodMass;
    if (scale!=e.scale) scale=-1;
  }

  friend std::ostream& operator<<(std::ostream& out, const eventCounter_t &e) {
    const char *line="-----------------------------------------------\n";
    out << line;
    out << "eventCounter info:\n";
    out << Form("   numEvents                = %lu\n",e.numEvents);
    out << Form("   numEventsPassedEvtTrigger= %lu\n",e.numEventsPassedEvtTrigger);
    out << Form("   numDielectronsUnweighted = %lu\n",e.numDielectronsUnweighted);
    out << Form("   numDielectrons           = %9.2lf\n",e.numDielectrons);
    out << Form("   numDielectronsGoodEta    = %9.2lf\n",e.numDielectronsGoodEta);
    out << Form("   numDielectronsGoodEt     = %9.2lf\n",e.numDielectronsGoodEt);
    out << Form("   numDielectronsHLTmatched = %9.2lf\n",e.numDielectronsHLTmatched);
    out << Form("   numDielectronsIDpassed   = %9.2lf\n",e.numDielectronsIDpassed);
    out << Form("   numDielectronsGoodMass   = %9.2lf\n",e.numDielectronsGoodMass);
    out << Form("   scale = %9.4e\n",e.scale);
    out << line;
    return out;
  }
};


#endif
