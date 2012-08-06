#ifndef HelpingTools_cc
#define HelpingTools_cc


#include <iostream>
#include <math.h>

#ifndef HelpingToolsCodeOnly
#  include "HelpingTools.hh"
#else
namespace helpers{ 
  void FileName(std::string&, std::string&);
  std::string SimplifyDielectronFlagInfo(const std::string &info);
};
#endif

#ifndef __myLib__
inline int PosOk(size_t pos) { return ((pos!=std::string::npos) ? 1:0); }
#endif

//#include <strings.h>

// ------------------------------------------------------

void HelpingTools() {}

// ------------------------------------------------------
/*
void helpers::FileName(const TString &str, TString &name) {
  if ( &str != &name ) { name=str; }
  Ssiz_t s=name.Last('/');
  name.Remove(0,s+1);
}
*/

// ------------------------------------------------------

namespace helpers {
void FileName(const std::string &str, std::string &name) {
  if ( &str != &name ) { name=str; }
  size_t s=name.find_last_of('/');
  name.erase(0,s+1);
}
};

// ------------------------------------------------------
/*
namespace helpers {
std::string SimplifyDielectronFlagInfoOld(const std::string &info_orig, int debug) {
  const int debug_loc=0;
  std::string info=info_orig;
  if (debug_loc) std::cout << "simplifying info=<" << info << ">\n";
  int pass=0;
  for (int ip=0; ip<3; ++ip) {
    const char *t=(ip==2) ? "tot" : ((ip) ? "pass":"fail");
    const int tsize=strlen(t);
    size_t p=info.find(t);
    if (PosOk(p)) {
      info.erase(p,tsize);
      p=info.find(t);
      if (PosOk(p)) info.erase(p,tsize);
      pass=(2*ip-1);
    }
  }
  if ((pass!=0) && debug_loc) std::cout << "modified info=<" << info << ">, pass=" << pass << "\n";
  std::string s;
  int ratio=0;
  size_t div=info.find('_');
  { 
    size_t tt=info.find("HRatio");
    if (PosOk(tt)) div=info.find('_',div+1);
  }
  if (PosOk(div)) {
    ratio=1;
    int idxOne1=0;
    int idxOne2=0;
    size_t p=info.find('x');
    if (!PosOk(p)) { p=info.find('s'); if (PosOk(p)) idxOne1=1; }
    if (!PosOk(p)) { if (debug) { std::cout << " failed to find 'x' or 's' in a string with '_' in " << info << "\n"; } return s; }
    std::string eb1=info.substr(p,4);
    size_t p1=info.find('x',p+1);
    if (!PosOk(p1)) { p1=info.find('s'); if (PosOk(p1)) idxOne2=1; }
    if (idxOne1!=idxOne2) {
      if (debug) { std::cout << " failed to recognize properly 'x'&'s' distinction in " << info << "\n"; }
      return s; 
    }
    if (!PosOk(p1)) { if (debug) { std::cout << " failed to find 'x' or 's' in a string with '_' and 1('x' or 's') in " << info << "\n"; } return s; }
    std::string eb2=info.substr(p1,4);
    if (eb1!=eb2) { if (debug) { std::cout << " eb mismatch in " << info << "\n"; } return s; }

    size_t wp80_1st_a=info.find("1stWP80");
    size_t wp80_1st_b=info.find("1stWP80",div);
    size_t wp80_a=info.find("WP80");
    size_t wp80_b=info.find("WP80",div);
    size_t hlt_1st_a=info.find("1stHLT");
    size_t hlt_1st_b=info.find("1stHLT",div);
    size_t hlt2_a=info.find("HLT2");
    size_t hlt2_b=info.find("HLT2",div);

    int ok=0;
    if (PosOk(wp80_a) && PosOk(wp80_1st_b) && (wp80_a!=wp80_b)) {
      if (PosOk(hlt_1st_a) && PosOk(hlt_1st_b) && (hlt_1st_a!=hlt_1st_b)) {
	s="T&P eff.WP80 ";
	ok=1;
      }
    }
    else if (PosOk(wp80_a) && PosOk(wp80_b) && !PosOk(wp80_1st_a) && (wp80_a!=wp80_b)) {
      if (PosOk(hlt_1st_a) && PosOk(hlt_1st_b) && (hlt_1st_a==hlt_1st_b)) {
	s="T&P eff.HLT ";
	ok=1;
      }
    }

    if (pass) {
      if (!PosOk(wp80_1st_a) && !PosOk(wp80_1st_b) && PosOk(wp80_a) && PosOk(wp80_b)) {
	if (PosOk(hlt_1st_a) && PosOk(hlt_1st_b)) {
	  s="T&P eff.WP80 ";
	  ok=1;
	}
	else if (PosOk(hlt2_a) && PosOk(hlt2_b)) {
	  s="T&P eff.HLT ";
	  ok=1;
	}
      }
    }

    if (!ok) {
      if (debug) { std::cout << " failed to find expected pattern in " << info << "\n"; }
      return s;
    }
  }

  if (!ratio) {
    if (PosOk(info.find("1stWP80"))) s="tag ";
    else if (PosOk(info.find("WP80"))) s="T&P80 ";
    else if (debug) { std::cout << "failed to find WP80 part of " << info << "\n"; }
    if (debug_loc) std::cout << "(not ratio) starting string s=<" << s << ">\n";
  }
  if (PosOk(info.find("xEXEX"))) s+="all ";
  else if (PosOk(info.find("xEXEB"))) s+="probeB ";
  else if (PosOk(info.find("xEXEE"))) s+="probeE ";
  else if (PosOk(info.find("xEBEX"))) s+="tagB ";
  else if (PosOk(info.find("xEEEX"))) s+="tagE ";
  else if (PosOk(info.find("xEBEB"))) s+="B-B ";
  else if (PosOk(info.find("xEBEE"))) s+="B-E ";
  else if (PosOk(info.find("xEEEB"))) s+="E-B ";
  else if (PosOk(info.find("xEEEE"))) s+="E-E ";
  else if (PosOk(info.find("sEBEE"))) s+="B&E ";
  else if (debug) { std::cout << "failed to find EB/EE part of " << info << "\n"; }
  if (!ratio) {
    if (PosOk(info.find("1stWP80"))) {
      s+="tagWP80";
    }
    else if (PosOk(info.find("1stHLT"))) {
      s+="tagHLT";
    }
    else if (PosOk(info.find("HLT2"))) s+="bothHLT";
  }
  if (debug_loc) std::cout << "simplifying info=<" << info << ">, extracted name=<" << s << ">\n";
  if (pass!=0) {
    switch(pass) {
    case -1: s+="fail"; break;
    case 0: break;
    case 1: s+="pass"; break;
    case 3: s+="tot"; break;
    default: s+="???";
    }
  }
  return s;
}
};
*/
// ------------------------------------------------------

namespace helpers {
std::string SimplifyDielectronFlagInfo(const std::string &info_orig, int debug) {
  const int debug_loc=0;
  std::string info=info_orig;
  if (debug_loc) std::cout << "simplifying info=<" << info << ">\n";
  int pass=0;
  for (int ip=0; ip<3; ++ip) {
    const char *t=(ip==2) ? "tot" : ((ip) ? "pass":"fail");
    const int tsize=strlen(t);
    size_t p=info.find(t);
    if (PosOk(p)) {
      info.erase(p,tsize);
      p=info.find(t);
      if (PosOk(p)) info.erase(p,tsize);
      pass=(2*ip-1);
    }
  }
  if ((pass!=0) && debug_loc) std::cout << "modified info=<" << info << ">, pass=" << pass << "\n";

  std::string eb;
  if (PosOk(info.find("xEXEX"))) {
    size_t p=info.find("_Eta");
    if (!p) eb="allEB ";
    else {
      size_t t=info.find('_',p+2);
      if (!PosOk(t)) t=info.size();
      eb=info.substr(p+1,t-p-1);
      eb+=std::string(" ");
      //std::cout << " eb=<" << eb << ">\n";
    }
  }
  else if (PosOk(info.find("xEXEB"))) eb="probeB ";
  else if (PosOk(info.find("xEXEE"))) eb="probeE ";
  else if (PosOk(info.find("xEBEX"))) eb="tagB ";
  else if (PosOk(info.find("xEEEX"))) eb="tagE ";
  else if (PosOk(info.find("xEBEB"))) eb="BB ";
  else if (PosOk(info.find("xEBEE"))) eb="BE ";
  else if (PosOk(info.find("xEEEB"))) eb="EB ";
  else if (PosOk(info.find("xEEEE"))) eb="EE ";
  else if (PosOk(info.find("sEBEE"))) eb="B&E ";
  else if (debug) { std::cout << "failed to find EB/EE part of " << info << "\n"; }

  std::vector<std::string> targets;
  targets.reserve(10);
  targets.push_back("xEXEX"); targets.push_back("xEXEB");
  targets.push_back("xEXEE"); targets.push_back("xEBEX");
  targets.push_back("xEEEX"); targets.push_back("xEBEB");
  targets.push_back("xEBEE"); targets.push_back("xEEEB");
  targets.push_back("xEEEE"); targets.push_back("sEBEE");
  for (unsigned int i=0; i<targets.size(); ++i) {
    size_t p=info.find(targets[i]);
    while (PosOk(p)) {
      info.erase(p,targets[i].size());
      p=info.find(targets[i]);
    }
  }

  std::string plotname;
  {
    size_t p=info.find("HRef");
    if (PosOk(p)) {
      size_t t=info.find(' ',p);
      //std::cout << "href: p=" << p << ", t=" << t << "\n";
      plotname=std::string("Cnt. ") + info.substr(p+4,t-3);
      info.erase(p,t-p);
    }
    else if (1) {
      p=info.find("HRatio");
      if (PosOk(p)) {
	size_t t=info.find('_',p);
	plotname=std::string("Eff. ") + info.substr(p+6,t-5);
	info.erase(p,t-p);
      }
    }
  }
  if (debug_loc) std::cout << "info=<" << info << ">, plotname=" << plotname << "\n";

  std::string s;
  //int ratio=0;
  size_t div=info.find('_');
  if (!PosOk(div)) div=0;

  if (debug_loc) std::cout << "info=<" << info << ">, div=" << div << "\n";
  if (!div) {
    if (PosOk(info.find("1stWP801stHLT"))) s="Tag";
    else if (PosOk(info.find("WP801stHLT"))) s="ProbeWP80";
    else if (PosOk(info.find("WP80HLT2"))) s="ProbeHLT";
    else s=info;
  }
  else {
    if (PosOk(info.find("1stWP801stHLT_WP801stHLT"))) s="WP80eff";
    else if (PosOk(info.find("WP801stHLT_WP80HLT2"))) s="HLTeff";
    else s=info;
  }

  if (debug_loc) std::cout << "simplifying info=<" << info << ">, extracted name=<" << s << ">\n";
  if (pass!=0) {
    switch(pass) {
    case -1: s+="fail"; break;
    case 0: break;
    case 1: s+="pass"; break;
    case 3: s+="tot"; break;
    default: s+="???";
    }
  }

  std::string ans=plotname + eb + s;
  if (debug_loc) std::cout << " final answer is <" << ans << ">\n";
  return ans;
}
};

// ------------------------------------------------
// ------------------------------------------------

// --------------------------------------------------------------

void PrintLines(std::ostream &out, const std::vector<std::string> &lines) {
  for (unsigned int i=0; i<lines.size(); ++i) {
    out << lines[i] << "\n";
  }
}

// --------------------------------------------------------------

void PrintLines(std::ostream &out, const std::vector<std::vector<std::string>*> &lines) {
  for (unsigned int i=0; i<lines.size(); ++i) {
    out << "\n\nlines from #" << i << "\n";
    PrintLines(out,*lines[i]);
  }
}

// --------------------------------------------------------------


void PrintLinesHTML(std::vector<std::string> &html, const std::vector<std::string> &lines, const char *description, int start_table, int end_table, int start_doc, int end_doc) {
  //if (lines.size()<2) {
  //  std::cout << "PrintLinesHTML: built-in restriction is lines.size>=2\n";
  //  return;
  //}

  std::cout << "PrintLinesHTML: got lines: \n"; PrintLines(lines);

  const int buflen=300;
  char buf[buflen];
  const char *title=(description) ? description : "tables";
  if (start_table || start_doc) {
    if (start_doc) {
      std::cout << "pushing the Doc type\n";
      html.push_back("<!DOCTYPE html\n    PUBLIC \"-//W3C//DTD HTML 3.2//EN\">\n");
      sprintf(buf,"<html>\n<head><title>%s</title></head>\n",title); html.push_back(buf);
      sprintf(buf,"<body bgcolor=\"EEEEEE\">\n<h3 style=\"text-align:left; color:DD6600;\">%s</h3>\n",title); html.push_back(buf);
    }
    html.push_back("<table border=\"0\" cellspacing=\"2\" width=\"100\%\">\n");
  }
  
  if (description) {
    sprintf(buf,"<tr><h4 style=\"text-align:left; color:DD6600;\">%s</h4></tr>\n",description); html.push_back(buf);
  }

  // determine widths
  //int count=0;
  size_t pos;
  std::string s;
  std::vector<double> wfracs;
  {
    unsigned int i0=0, l0=0;
    unsigned int i1=1, l1=0;
    for (unsigned int i=0; i<lines.size(); ++i) {
      unsigned int sz=lines[i].size();
      if (sz>l0) {
	l1=l0; i1=i0;
	l0=sz; i0=i;
      }
      else if (sz>l1) {
	l1=sz; i1=i;
      }
    }

    std::vector<size_t> sizes0, sizes1;
    std::vector<int> widths;
    for (int i=0; i<2; ++i) {
      std::vector<size_t> *sz=(i==0) ? &sizes0 : &sizes1;
      s=(i==0) ? lines[i0] : lines[i1];
      pos=0; sz->push_back(0);
      while (PosOk(pos)) {
	pos=s.find('|',pos+1);
	if (PosOk(pos)) sz->push_back(pos);
      }
    }

    if (sizes0.size()!=sizes1.size()) {
      const std::vector<size_t> *sz=(sizes1.size()>sizes0.size()) ? &sizes1:&sizes0;
      widths.reserve(sz->size()-1);
      for (unsigned int i=1; i<sz->size(); ++i) {
	widths.push_back((*sz)[i]-(*sz)[i-1]);
      }
    }
    else {
      widths.reserve(sizes0.size()-1);
      for (unsigned int i=1; i<sizes0.size(); ++i) {
	int w0=sizes0[i]-sizes0[i-1];
	int w1=sizes1[i]-sizes1[i-1];
	widths.push_back((w0>w1) ? w0:w1);
      }
    }
    wfracs.reserve(widths.size());
    pos=(lines[0].size()>lines[1].size()) ? lines[0].size() : lines[1].size();
    for (unsigned int i=0; i<widths.size(); ++i) {
      //std::cout << "widths[" << i << "]=" << widths[i] << "\n";
      wfracs.push_back(trunc(widths[i]*100/double(pos)));
    }
  }

  for (unsigned int i=0; i<lines.size(); ++i) {
    s= lines[i];
    sprintf(buf,"<tr>\n"); html.push_back(buf);
    pos=s.find('|');
    if (!PosOk(pos)) {
      std::cout << " ! plain print " << lines[i] << "\n";
      sprintf(buf,"%s",lines[i].c_str()); html.push_back(buf);
    }
    else {
      std::string sb;
      size_t pos1;
      int corrected=0;
      int idx=0;
      while (PosOk(pos) && !corrected) {
	pos1=s.find('|',pos+1);
	if (PosOk(pos1)) { //{ pos1=s.size(); corrected=1; }
	  pos++;
	  //pos=s.find_first_not_of(' ',pos);
	  sb=s.substr(pos,pos1-pos-1);
	  //pos=sb.find_last_not_of(' ');
	  //sb.erase(pos+1);
	  sprintf(buf,"<td width=\"%4.1lf%%\">%s</td>\n",wfracs[idx],sb.c_str());
	  html.push_back(buf);
	}
	pos=pos1;
	idx++;
      }
    }
    html.push_back("</tr>\n");
  }

  if (end_table || end_doc) {
    html.push_back("</table>\n");
    html.push_back("<hr />\n");
    if (end_doc) {
      html.push_back("</body>\n");
      html.push_back("</html>\n");
    }
  }

  std::cout << "PrintLinesHTML: produced lines: \n"; PrintLines(lines);

  return;
}

// --------------------------------------------------------------


void PrintLinesHTML(std::ostream &htmlfile, const std::vector<std::string> &lines, const char *description, int start_table, int end_table, int start_doc, int end_doc) {
  //if (lines.size()<2) {
  //  std::cout << "PrintLinesHTML: built-in restriction is lines.size>=2\n";
  //  return;
  //}
  if (start_table || start_doc) {
    if (start_doc) {
      htmlfile << "<!DOCTYPE html" << std::endl;
      htmlfile << "    PUBLIC \"-//W3C//DTD HTML 3.2//EN\">" << std::endl;
      htmlfile << "<html>" << std::endl;
      const char *title=(description) ? description : "tables";
      htmlfile << "<head><title>" << title << "</title></head>" << std::endl;
      htmlfile << "<body bgcolor=\"EEEEEE\">" << std::endl;
      htmlfile << "<h3 style=\"text-align:left; color:DD6600;\">" << title << "</h3>" << std::endl;
      htmlfile << std::endl;
    }
    htmlfile << "<table border=\"0\" cellspacing=\"2\" width=\"100%\">" << std::endl;
  }
  
  if (description) {
    htmlfile << "<tr><h4 style=\"text-align:left; color:DD6600;\">" << description << "</h4></tr>" << std::endl;
  }

  // determine widths
  //int count=0;
  size_t pos;
  std::string s;
  std::vector<double> wfracs;
  {
    int i0=0, l0=0;
    int i1=1, l1=0;
    for (unsigned int i=0; i<lines.size(); ++i) {
      int sz=lines[i].size();
      if (sz>l0) {
	l1=l0; i1=i0;
	l0=sz; i0=i;
      }
      else if (sz>l1) {
	l1=sz; i1=i;
      }
    }

    std::vector<size_t> sizes0, sizes1;
    std::vector<int> widths;
    for (int i=0; i<2; ++i) {
      std::vector<size_t> *sz=(i==0) ? &sizes0 : &sizes1;
      s=(i==0) ? lines[i0] : lines[i1];
      pos=0; sz->push_back(0);
      while (PosOk(pos)) {
	pos=s.find('|',pos+1);
	if (PosOk(pos)) sz->push_back(pos);
      }
    }

    if (sizes0.size()!=sizes1.size()) {
      const std::vector<size_t> *sz=(sizes1.size()>sizes0.size()) ? &sizes1:&sizes0;
      widths.reserve(sz->size()-1);
      for (unsigned int i=1; i<sz->size(); ++i) {
	widths.push_back((*sz)[i]-(*sz)[i-1]);
      }
    }
    else {
      widths.reserve(sizes0.size()-1);
      for (unsigned int i=1; i<sizes0.size(); ++i) {
	int w0=sizes0[i]-sizes0[i-1];
	int w1=sizes1[i]-sizes1[i-1];
	widths.push_back((w0>w1) ? w0:w1);
      }
    }
    wfracs.reserve(widths.size());
    pos=(lines[0].size()>lines[1].size()) ? lines[0].size() : lines[1].size();
    for (unsigned int i=0; i<widths.size(); ++i) {
      //std::cout << "widths[" << i << "]=" << widths[i] << "\n";
      wfracs.push_back(trunc(widths[i]*100/double(pos)));
    }
  }

  for (unsigned int i=0; i<lines.size(); ++i) {
    s= lines[i];
    htmlfile << "<tr>" << std::endl;
    pos=s.find('|');
    if (!PosOk(pos)) htmlfile << lines[i] << std::endl;
    else {
      std::string sb;
      size_t pos1;
      int corrected=0;
      int idx=0;
      while (PosOk(pos) && !corrected) {
	pos1=s.find('|',pos+1);
	if (PosOk(pos1)) { //{ pos1=s.size(); corrected=1; }
	  pos++;
	  //pos=s.find_first_not_of(' ',pos);
	  sb=s.substr(pos,pos1-pos-1);
	  //pos=sb.find_last_not_of(' ');
	  //sb.erase(pos+1);
	  htmlfile << "<td width=\"" << wfracs[idx] << "%\">" << sb << "</td>" << std::endl;
	}
	pos=pos1;
	idx++;
      }
    }
    htmlfile << "</tr>" << std::endl;
  }

  if (end_table || end_doc) {
    htmlfile << "</table>" << std::endl;
    htmlfile << "<hr />" << std::endl;
    if (end_doc) {
      htmlfile << "</body>" << std::endl;
      htmlfile << "</html>" << std::endl;
    }
  }
  return;
}

// --------------------------------------------------------------

void PrintLinesHTML(std::ostream &out, const std::vector<std::vector<std::string>*> &lines, const std::vector<std::string> *descriptions) {
  for (unsigned int i=0; i<lines.size(); ++i) {
    const char *descr=(descriptions) ? (*descriptions)[i].c_str() : NULL;
    const int start_table=(i==0)?1:0;
    const int end_table=(i+1==lines.size()) ? 1:0;
    PrintLinesHTML(out,*lines[i],descr,start_table,end_table);
  }
}

// ------------------------------------------------
// ------------------------------------------------

int RemoveWrongStringFromFile(const char *fname, const char *find_str, const char *replace) {
  const char *fncname="RemoveWrongStringFromFile";
  std::ifstream fin;
  fin.open(fname);
  if (!fin) {
    std::cout << fncname << " failed to open <" << fname << ">\n";
    return 0;
  }
  std::vector<std::string> lines;
  lines.reserve(2000);
  std::string s;
  size_t pos=0;
  while (!fin.eof() && getline(fin,s)) {
    do {
      pos=s.find(find_str);
      if (pos!=std::string::npos) {
	std::string tmp;
	if (replace && (strlen(replace)>0)) tmp=s.substr(0,pos) + std::string(replace) + s.substr(pos+strlen(find_str));
	else tmp=s.substr(0,pos) + s.substr(pos+strlen(find_str));
	//std::cout << "got string {" << s << "}, modified to {" << tmp << "}\n";
	s=tmp;
      }
    } while (pos!=std::string::npos);
    lines.push_back(s);
  }
  fin.close();
  // recreate the file
  std::ofstream fout;
  fout.open(fname);
  if (!fout) {
    std::cout << fncname << " failed to recreate the file <" << fname << ">\n";
    return 0;
  }
  for (unsigned int i=0; i<lines.size(); ++i) {
    fout << lines[i] << "\n";
  }
  fout.close();
  return 1;
}

// ------------------------------------------------


#endif
