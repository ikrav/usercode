#ifndef EtaEtaMass_H
#define EtaEtaMass_H

#include <TObject.h>
#include <ostream>

// ---------------------------------------------------------------
// Helper class
// ---------------------------------------------------------------

class EtaEtaMassData_t : public TObject {
public:
  EtaEtaMassData_t() : _Eta1(0.), _Eta2(0.), _Mass(0.), _Weight(1.),
		       _nGoodPV(0) {}
  ~EtaEtaMassData_t(){}

  double eta1() const { return _Eta1; }
  double eta2() const { return _Eta2; }
  double mass() const { return _Mass; }
  void mass(double m) { _Mass=m; }
  int nGoodPV() const { return _nGoodPV; }
  void nGoodPV(int nPV) { _nGoodPV=nPV; }
  double weight() const { return _Weight; }
  void weight(double w) { _Weight=w; }
  void Assign(const EtaEtaMassData_t &d) { _Eta1=d._Eta1; _Eta2=d._Eta2; _Mass=d._Mass; _Weight=d._Weight; _nGoodPV=d._nGoodPV; }
  void Assign(double eta1val, double eta2val, double massVal, double weightVal, int nPV) { _Eta1=eta1val; _Eta2=eta2val; _Mass=massVal; _Weight=weightVal; _nGoodPV=nPV; }

  friend std::ostream& operator<<(std::ostream& out, const EtaEtaMassData_t &data) {
    out << "EEM(" << data._Eta1 << ", " << data._Eta2 << ", " << data._Mass << ", " << data._Weight << "; nPV=" << data._nGoodPV << ")";
    return out;
  }

  // data fields
  double _Eta1,_Eta2,_Mass,_Weight;
  int _nGoodPV;

  ClassDef(EtaEtaMassData_t,3)
};


#endif
