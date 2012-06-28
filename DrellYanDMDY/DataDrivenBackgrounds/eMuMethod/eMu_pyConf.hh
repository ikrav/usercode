#include "eMu.hh"

//BOOST PYTHON HEADERS
#include <boost/python/module.hpp>
#include <boost/python/def.hpp>
#include <boost/python/init.hpp>
#include <boost/python/class.hpp>

class eMu_If : public eMu {
  public:
  eMu_If();
  int help();
};

BOOST_PYTHON_MODULE(bkg)
{
  using namespace boost::python;

  class_<eMu>("eMu")
    .def(init<>())
    //    .def(init<int, char**>())
    .def("run", &eMu::run)
    .def("setMCreWeight", &eMu::setMCreWeight)
    .def("doDMDYanal", &eMu::doDMDYanal)
    .def("setSaveToRootFile", &eMu::setSaveToRootFile)
    .def("runPUreWeighting", &eMu::runPUreWeighting)
    .def_readwrite("emuNtupleDir", &eMu::emuNtupleDir)
    .def_readwrite("eeNtupleDir", &eMu::eeNtupleDir)    
    .def_readwrite("filePostfix", &eMu::filePostfix)
    .def_readwrite("lumiVal", &eMu::lumiVal)
    .def_readwrite("nBins", &eMu::nBins)
    .def_readwrite("subDir", &eMu::subDir)
    .def_readwrite("xmax", &eMu::xmax)
    .def_readwrite("xmin", &eMu::xmin)
    
    //    .def_readwrite("", &eMu::)
    ;

  class_<eMu_If, bases<eMu> >("eMu_If")
    .def(init<>())
    //    .def(init<int, char**>())
    .def("run", &eMu_If::run)
    .def("setMCreWeight", &eMu_If::setMCreWeight)
    .def("doDMDYanal", &eMu_If::doDMDYanal)
    .def("setSaveToRootFile", &eMu_If::setSaveToRootFile)
    .def("runPUreWeighting", &eMu_If::runPUreWeighting)
    .def_readwrite("emuNtupleDir", &eMu_If::emuNtupleDir)
    .def_readwrite("eeNtupleDir", &eMu_If::eeNtupleDir)    
    .def_readwrite("filePostfix", &eMu_If::filePostfix)
    .def("help", &eMu_If::help)
    .def_readwrite("lumiVal", &eMu_If::lumiVal)
    .def_readwrite("nBins", &eMu_If::nBins)
    .def_readwrite("subDir", &eMu_If::subDir)
    .def_readwrite("xmax", &eMu_If::xmax)
    .def_readwrite("xmin", &eMu_If::xmin)
    
    //    .def_readwrite("", &eMu_If::)
    ;
}


