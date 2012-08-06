#ifndef HelpingTools_H
#define HelpingTools_H

#include <TString.h>
#include <string>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <vector>

//void HelpingTools();

//inline int PosOk(size_t pos) { return ((pos!=std::string::npos) ? 1:0); }

namespace helpers {

  void FileName(const std::string &str, std::string &name);

  inline void FileName(const TString &str, TString &name) {
    std::string s1=str.Data();
    std::string n1;
    FileName(s1,n1);
    name=n1;
  }

  std::string SimplifyDielectronFlagInfo(const std::string &info, int debug=1);

  inline std::string SimplifyDielectronFlagInfo(const TString &info, int debug=1) {
    std::string s=info.Data();
    return SimplifyDielectronFlagInfo(s,debug);
  }

};


// --------------------------------------------------------------

void PrintLines(std::ostream &out, const std::vector<std::string> &lines);
void PrintLines(std::ostream &out, const std::vector<std::vector<std::string>*> &lines);

inline void PrintLines(const std::vector<std::string> &lines) { return PrintLines(std::cout,lines); }
inline void PrintLines(const std::vector<std::vector<std::string>*> &lines) { PrintLines(std::cout,lines); }

// --------------------------------------------------------------


void PrintLinesHTML(std::vector<std::string> &html, const std::vector<std::string> &lines, const char *description=NULL, int start_table=1, int end_table=1, int start_doc=1, int end_doc=1); 

void PrintLinesHTML(std::ostream &out, const std::vector<std::string> &lines, const char *description=NULL, int start_table=1, int end_table=1, int start_doc=1, int end_doc=1); 
void PrintLinesHTML(std::ostream &out, const std::vector<std::vector<std::string>*> &lines, const std::vector<std::string> *descriptions);

inline void PrintLinesHTML(const std::vector<std::string> &lines, const std::string *description=NULL) { PrintLinesHTML(std::cout,lines,(description)?description->c_str():NULL); }
inline void PrintLinesHTML(const std::vector<std::vector<std::string>*> &lines, const std::vector<std::string> *names=NULL) { PrintLinesHTML(std::cout,lines,names); }

// --------------------------------------------------------------

template<class HistClass_t>
inline void SPrint(const HistClass_t &h, std::vector<std::string> &lines, int print_range, int print_error=0, int accuracy=2, int print_name=0, const char *print_label=NULL)  {
  std::cout << "SPrint: print_name=" << print_name << "\n";
  const int bin_count=h.GetNbinsX();
  //std::cout << "Contents of " << h.GetName() << ":\n";
  //std::cout << " (" << helpers::SimplifyDielectronFlagInfo(std::string(h.GetName())) << ")\n";

  int imin=(print_name) ? 1:0;
  if (print_name && print_label) imin++;

  if (int(lines.size())<bin_count) { std::string s="|  "; lines.insert(lines.end(),bin_count+1+imin,s); }
  if (print_name && print_label) {
    if (print_range) lines[0] += "       |  ";
    lines[0] += print_label;
    lines[0] +="  |";
  }
  if (print_name) {
    if (print_range) lines[imin-1] += "       |  ";
    lines[imin-1] += h.GetName();
    lines[imin-1] +="  |";
  }
  if (print_range) lines[imin] += "  range  |";
  lines[imin] += " ";
  lines[imin] += helpers::SimplifyDielectronFlagInfo(std::string(h.GetName()));
  lines[imin] += " |";
  char buf[50];
  char format[40], formatErr[40];

  int less100=1;
  int less1=1;
  for (Int_t ibin=1; (ibin<bin_count+1); ++ibin) {
    double v=h.GetBinContent(ibin);
    if (v>100) less100=0;
    if (less1 && (v>1.)) less1=0;
  }
  sprintf(format,(less100) ? "  %%7.%dlf |" : "  %%8.%dlf |",accuracy);
  sprintf(formatErr,"  %%7.%dlf+/-%%7.%dlf |",accuracy,accuracy);
  if (less1) {
    sprintf(format,"  %%6.4lf |");
    sprintf(formatErr,"  %%6.4lf+/-%%6.4lf |");
  }

  for (Int_t ibin=1; ibin<bin_count+1; ibin++) {
    const double hval= h.GetBinContent(ibin);
    if (print_range) {
      sprintf(buf,"  %5.1lf -- %5.1lf  |", h.GetBinLowEdge(ibin), (h.GetBinLowEdge(ibin)+h.GetBinWidth(ibin)));
      lines[ibin+imin] += buf;
    }
    if (print_error) {
      sprintf(buf, formatErr, hval, h.GetBinError(ibin));
      lines[ibin+imin] += buf;
    }
    else {
      sprintf(buf, format, hval);
      lines[ibin+imin] += buf;
    }
  }
}

// --------------------------------------------------------------

template<class HistClass_t>
inline void SPrint(const std::vector<HistClass_t*> &hV, std::vector<std::string> &lines, int print_range, int print_error=0, int accuracy=2, int print_name=0)  {
  for (int i=0; i<hV.size(); ++i) {
    int prnRange=(print_range && (i==0))?1:0;
    SPrint(*hV[i],lines,prnRange,print_error,accuracy,print_name);
  }
}

// --------------------------------------------------------------

template<class HistClass_t>
inline void SPrint(const std::vector<HistClass_t*> &hV, std::vector<std::string> &lines, int print_range, const std::vector<int> *print_error=NULL, const std::vector<int>* accuracy=NULL, int print_name=0)  {
  for (int i=0; i<hV.size(); ++i) {
    int prnRange=(print_range && (i==0))? 1:0;
    int prnErr=(print_error && (*print_error)[i]) ? 1:0;
    int prnAcc=(!accuracy) ? 2:(*accuracy)[i];
    SPrint(*hV[i],lines,prnRange,prnErr,prnAcc,print_name);
  }
}

// --------------------------------------------------------------

template<class HistoClass_t>
inline int Normalize(const std::vector<HistoClass_t*> &hV, std::vector<HistoClass_t*> &hNormV, std::vector<double> &norms, int do_normalize=1) {
  const unsigned int n=hV.size();
  hNormV.clear(); norms.clear();
  hNormV.reserve(n); norms.reserve(n);
  for (unsigned int i=0; i<n; ++i) {
    const HistoClass_t *h=hV[i];
    TString title=h->GetTitle();
    TString name=h->GetName();
    title.Append(" normed");
    name.Append("normed");
    norms.push_back(h->Integral());
    hNormV.push_back((HistoClass_t*)h->Clone(name.Data()));
    hNormV.back()->SetTitle(title.Data());
    if (do_normalize) hNormV.back()->Scale(1/norms.back());
  }
  return (hV.size())?1:0;
}

// --------------------------------------------------------------

template<class HistoClass_t>
inline int NormalizeTo(const std::vector<HistoClass_t*> &hV, const std::vector<HistoClass_t*> &hNormalizationTarget, std::vector<HistoClass_t*> &hNormV, std::vector<double> &norms, int do_normalize=1) {
  const unsigned int n=hV.size();
  hNormV.clear(); norms.clear();
  if (hNormalizationTarget.size()!=n) {
    std::cout << "NormalizeTo size difference: hV.size=" << hV.size() << ", hNormalizationTarget.size=" << hNormalizationTarget.size() << "\n";
    return 0;
  }
  hNormV.reserve(n); norms.reserve(n);
  for (unsigned int i=0; i<n; ++i) {
    const HistoClass_t *h=hV[i];
    TString title=h->GetTitle();
    TString name=h->GetName();
    title.Append(" normed");
    name.Append("normed");
    double the_norm=h->Integral();
    double scale_to_norm=hNormalizationTarget[i]->Integral();
    if (fabs(scale_to_norm)<1e-4) {
      std::cout << "Warning requested to normalize to " << scale_to_norm << " of a histogram <" << hNormalizationTarget[i]->GetName() << ">. Setting to 1\n";
      scale_to_norm=1;
    }
    double scale=(the_norm==0) ? 0 : scale_to_norm/the_norm;
    norms.push_back(the_norm);
    hNormV.push_back((HistoClass_t*)h->Clone(name.Data()));
    hNormV.back()->SetTitle(title.Data());
    if (do_normalize) hNormV.back()->Scale(scale);
  }
  return (hV.size())?1:0;
}

// --------------------------------------------------------------

template<class HistoClass_t>
inline int ChangeToAbsValsV(std::vector<HistoClass_t*> &hV) {
  for (unsigned int t=0; t<hV.size(); ++t) {
    HistoClass_t* h=hV[t];
    for (int i=1; i<h->GetNbinsX(); ++i) {
      double x=h->GetBinContent(i);
      if (x<0) {
	h->SetBinContent(i,fabs(x));
      }
    }
  }
  return (hV.size())?1:0;
}

// --------------------------------------------------------------

template<class HistoClass_t>
inline void CalculateMaxValues(const HistoClass_t &h, double & maxAbsY, double & sumY2) {
  const int bin_count=h.GetNbinsX();
  sumY2=0;
  for (int i=1; i<bin_count; ++i) {
    double y=h.GetBinContent(i);
    if (fabs(y)>fabs(maxAbsY)) maxAbsY=y;
    sumY2+= y*y;
  }
  return;
}

// --------------------------------------------------------------

template<class HistoClass_t>
inline int GetAsymmetryV(const std::vector<HistoClass_t*> &h1V, const std::vector<HistoClass_t*> &h2V, std::vector<HistoClass_t*> &asymV, int abs_vals=0) {
  const unsigned int n=h1V.size();
  if (n!=h2V.size()) {
    std::cout << "GetAsymmetry(V,V,V): vec1.size=" << n << ", vec2.size=" << h2V.size() << "\n";
    return 0;
  }
  asymV.clear(); asymV.reserve(n);
  for (unsigned int i=0; i<n; ++i) {
    double relw=h1V[i]->Integral() / h2V[i]->Integral();
    HistoClass_t *hnew=(HistoClass_t*)h1V[i]->GetAsymmetry(h2V[i],relw);
    hnew->GetYaxis()->SetTitle("asymmetry");
    asymV.push_back(hnew);
  }
  if (abs_vals) return ChangeToAbsValsV(asymV);
  return 1;
}

// --------------------------------------------------------------

inline void Pause(const char *msg=NULL) {
  char c;
  if (msg) std::cout << msg << "\n";
  std::cout << "waiting for a char...\n";
  std::cin >> c;
  if (c=='q') exit(0);
  std::cout << "Pause got c=" <<c << "\n";
}

// --------------------------------------------------------------

inline
std::ostream& operator<<(std::ostream &out, const TString &s) {
  out << s.Data();
  return out;
}

//------------------------------------------------------------------------------------------------------------------------

//template<class T>
//inline void PrintVec(const char *msg, const std::vector<T> &values);

template<class T>
inline void PrintVec(const char *msg, const std::vector<T> &arr) {
  if (msg) std::cout << msg;
  std::cout << " vec[" << arr.size() << "]: ";
  for (unsigned int i=0; i<arr.size(); ++i) {
    std::cout << " " << arr[i];
  }
  std::cout << "\n";
  return;
}

//------------------------------------------------------------------------------------------------------------------------

template<class T>
inline void PrintRooVec(const char *msg, const std::vector<T*> &arr) {
  if (msg) std::cout << msg;
  std::cout << " vec[" << arr.size() << "]: ";
  for (unsigned int i=0; i<arr.size(); ++i) {
    std::cout << " " << arr[i]->getVal();
  }
  std::cout << "\n";
  return;
}

//------------------------------------------------------------------------------------------------------------------------

template<class T, class T2>
  inline void PrintVec2(const char *msg, const std::vector<T> &arr, const std::vector<T2> &arr2) {
  if (msg) std::cout << msg;
  std::cout << " vec[" << arr.size() << "]: ";
  if (arr.size()!=arr2.size()) {
    std::cout << "size difference (" << arr.size() << " vs " << arr2.size() << ")\n";
    return;
  }
  for (unsigned int i=0; i<arr.size(); ++i) {
    std::cout << " (" << arr[i] << ',' << arr2[i] << ')';
  }
  std::cout << "\n";
  return;
}

//------------------------------------------------------------------------------------------------------------------------

template<class T>
inline void PrintVVecCounts(const char *msg, const std::vector<std::vector<T>*> &arr) {
  if (msg) std::cout << msg;
  std::cout << " vec[" << arr.size() << "] of vecs sizes: ";
  for (unsigned int i=0; i<arr.size(); ++i) {
    std::cout << " " << arr[i]->size();
  }
  std::cout << "\n";
  return;
}

//------------------------------------------------------------------------------------------------------------------------

template<class T>
inline
void ClearVec(std::vector<T*> &vec) {
  for (unsigned int i=0; i<vec.size(); ++i) delete vec[i];
  vec.clear();
}

//------------------------------------------------------------------------------------------------------------------------

template<class T>
inline
int AddToVec(std::vector<std::vector<T>*> &res, const std::vector<std::vector<T>*> &add) {
  if (res.size()!=add.size()) {
    std::cout << "AddToVec sizes: " << res.size() << " vs " << add.size() << "\n";
    return 0;
  }
  for (unsigned int i=0; i<res.size(); ++i) {
    std::vector<T>* vals= res[i];
    const std::vector<T> *av= add[i];
    if (av->size()) {
      vals->reserve( vals->size() + av->size() );
      vals->insert(vals->end(), av->begin(), av->end());
    }
  }
  return 1;
}

//------------------------------------------------------------------------------------------------------------------------

int RemoveWrongStringFromFile(const char *fname, const char *find_str, const char *replace=NULL);

//------------------------------------------------------------------------------------------------------------------------

template<class T>
inline void SWAP(T &a, T &b) { T c=a; a=b; b=c; }

template<class T>
inline void SWAP(T &a, T &b,  T &tmp) { tmp=a; a=b; b=tmp; }

//------------------------------------------------------------------------------------------------------------------------

template<class Type1_t, class Type2_t>
inline
int QuickSort(std::vector<Type1_t> &arr1, std::vector<Type2_t> &arr2, int calling=0, int print_err=1)
{
  // adapted from Numerical Recipes in C

  unsigned int count=arr1.size();
  if (arr1.size() != arr2.size()) {
    std::cout << "QuickSort: arr1.size=" << arr1.size() << ", arr2.size=" << arr2.size() << "\n";
    return 0;
  }

  const int insertion_sort_limit= 7;
  const int aux_index_size=100;
  if (count==0) return 1;

  long i,ir=count-1, j,k,l=0;
  std::vector<unsigned int> istack(aux_index_size);
  int jstack=0;
  Type1_t a1,tmp1;
  Type2_t a2,tmp2;
  
  for (;;) {
    if (ir-l<insertion_sort_limit) { // insertion sort when subarray is small enough
      for (j=l+1; j<=ir; ++j) {
	a1=arr1[j]; a2=arr2[j];
	for (i=j-1; i>=0; --i) {
	  if (a1>=arr1[i]) break;
	  arr1[i+1]=arr1[i];
	  arr2[i+1]=arr2[i];
	}
	arr1[i+1]=a1; arr2[i+1]=a2;
      }
      if (jstack<1) return 1;
      ir=istack[jstack-1]; // pop stack and begin a new round of partitioning
      l=istack[jstack-2];
      jstack -= 2;
    }
    else {
      // choose median of left, center and right elements as partitioning element a. Also rearrange so that a[l+1]<=a[l]<=a[ir]      
      k=(l+ir)>>1;
      SWAP(arr1[k],arr1[l+1],tmp1);
      SWAP(arr2[k],arr2[l+1],tmp2);
      if (arr1[l+1]>arr1[ir]) {
	SWAP(arr1[l+1],arr1[ir],tmp1);
	SWAP(arr2[l+1],arr2[ir],tmp2);
      }
      if (arr1[l]>arr1[ir]) { 
	SWAP(arr1[l],arr1[ir],tmp1);
	SWAP(arr2[l],arr2[ir],tmp2);
      }
      if (arr1[l+1]>arr1[l]) {
	SWAP(arr1[l+1],arr1[l],tmp1);
	SWAP(arr2[l+1],arr2[l],tmp2);
      }
      // initialize pointers of partitioning
      i=l+1;
      j=ir;
      // partitioning element
      a1=arr1[l]; a2=arr2[l];
      // innermost loop
      for (;;) {
	do ++i; while (a1>arr1[i]); // scan up to find element >a1
	do --j; while (arr1[j]>a1); // scan down to find element <a1
	if (j<i) break; // pointers crossed. Partitioning complete
	SWAP(arr1[i],arr1[j],tmp1);
	SWAP(arr2[i],arr2[j],tmp2);
      }
      arr1[l]=arr1[j]; arr1[j]=a1;
      arr2[l]=arr2[j]; arr2[j]=a2;
      jstack+=2;
      // push pointers to larger subarray on stack, process smaller subarray immediately.
      if (jstack>aux_index_size) {
	if (calling<5) return QuickSort(arr1,arr2,calling+1,print_err);
	if (print_err) printf("QuickSort: aux_index_size is too small in QSort(cnt=%d,type*,type*)",count);
	//return (print_err) ? WAM_infunction("QuickSort") : 0;
	return 0;
      }
      if (ir-l+1>=j-l) {
	istack[jstack-1]=ir;
	istack[jstack-2]=i;
	ir=j-1;
      }
      else {
	istack[jstack-1]=j-1;
	istack[jstack-2]=l;
	l=i;
      }
    }
  }
  return 1;
}

//------------------------------------------------------------------------------------------------------------------------

#endif
