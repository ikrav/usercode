#ifndef EtaEtaMass_H
#define EtaEtaMass_H

#include <TObject.h>
#include <ostream>

// ---------------------------------------------------------------
// Helper class
// ---------------------------------------------------------------

class EtaEtaMassData_t : public TObject {
public:
  EtaEtaMassData_t() : _Eta1(0.), _Eta2(0.), _Mass(0.), _Weight(1.) {}
  ~EtaEtaMassData_t(){}

  double eta1() const { return _Eta1; }
  double eta2() const { return _Eta2; }
  double mass() const { return _Mass; }
  void mass(double m) { _Mass=m; }
  double weight() const { return _Weight; }
  void weight(double w) { _Weight=w; }
  void Assign(const EtaEtaMassData_t &d) { _Eta1=d._Eta1; _Eta2=d._Eta2; _Mass=d._Mass; _Weight=d._Weight; }
  void Assign(double eta1val, double eta2val, double massVal, double weightVal=1.) { _Eta1=eta1val; _Eta2=eta2val; _Mass=massVal; _Weight=weightVal; }

  friend std::ostream& operator<<(std::ostream& out, const EtaEtaMassData_t &data) {
    out << "EEM(" << data._Eta1 << ", " << data._Eta2 << ", " << data._Mass << ", " << data._Weight << ")";
    return out;
  }

  // data fields
  double _Eta1,_Eta2,_Mass,_Weight;

  ClassDef(EtaEtaMassData_t,2)
};


#endif
