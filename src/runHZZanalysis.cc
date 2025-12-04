#include <cstdlib>
#include <iostream>

#include <boost/algorithm/string.hpp>
#include <boost/program_options.hpp>

#include <DileptonTrees.h>
#include <TrileptonTrees.h>
#include <Logger.h>
#include <Looper.h>
#include <NrbAnalysis.h>
#include <NrbTrees.h>
#include <Options.h>
#include <PhotonTrees.h>
#include <Version.h>
#include <ZGammaTrees.h>
#include <ElectronTrees.h>
#include <EGammaFromMisid.h>

namespace po = boost::program_options;


enum class AnalysisType {
  DileptonTrees,
  TrileptonTrees,
  PhotonTrees,
  NRB,
  NrbTrees,
  ZGammaTrees,
  ElectronTrees,
  EGammaFromMisid
};


template<typename T>
void runAnalysis(int argc, char **argv,
                 po::options_description const &commonOptions) {
  Options options(
      argc, argv, {commonOptions, Looper<T>::OptionsDescription()});
  LOG_DEBUG << "Version: " << Version::Commit();
  Looper<T> looper{options};
  looper.Run();
}


int main(int argc, char **argv) {
  po::options_description analysisTypeOptions{"Analysis type"};
  analysisTypeOptions.add_options()
    ("analysis,a", po::value<std::string>()->default_value("DileptonTrees"),
     "Analysis to run; allowed values are \"DileptonTrees\", "
     "\"PhotonTrees\", \"NRB\", \"ZGammaTrees\", \"TrileptonTrees\", \"NrbTrees\", \"ElectronTrees\", \"EGammaFromMisid\"");

  // Command line options are checked twice. At the first pass only check the
  // analysis type and update the list of expected options accordingly.
  po::variables_map vm;
  po::store(po::command_line_parser(argc, argv).options(analysisTypeOptions)\
      .allow_unregistered().run(), vm);
  po::notify(vm);

  std::string analysisTypeArg{vm["analysis"].as<std::string>()};
  boost::to_lower(analysisTypeArg);
  AnalysisType analysisType;

  if (analysisTypeArg == "dileptontrees")
    analysisType = AnalysisType::DileptonTrees;
  else if (analysisTypeArg == "trileptontrees")
    analysisType = AnalysisType::TrileptonTrees;
  else if (analysisTypeArg == "photontrees")
    analysisType = AnalysisType::PhotonTrees;
  else if (analysisTypeArg == "nrb")
    analysisType = AnalysisType::NRB;
  else if (analysisTypeArg == "nrbtrees")
    analysisType = AnalysisType::NrbTrees;
  else if (analysisTypeArg == "zgammatrees")
    analysisType = AnalysisType::ZGammaTrees;
  else if (analysisTypeArg == "electrontrees")
    analysisType = AnalysisType::ElectronTrees;
  else if (analysisTypeArg == "egammafrommisid")
    analysisType = AnalysisType::EGammaFromMisid;
  else {
    LOG_ERROR << "Unknown analysis type \"" <<
      vm["analysis"].as<std::string>() << "\"";
    std::exit(EXIT_FAILURE);
  }


  // Perform the second pass over command line options and run the requested
  // analysis
  switch (analysisType) {
    case AnalysisType::DileptonTrees:
      runAnalysis<DileptonTrees>(argc, argv, analysisTypeOptions);
      break;
    
    case AnalysisType::TrileptonTrees:
      runAnalysis<TrileptonTrees>(argc, argv, analysisTypeOptions);
      break;

    case AnalysisType::PhotonTrees:
      runAnalysis<PhotonTrees>(argc, argv, analysisTypeOptions);
      break;

    case AnalysisType::NRB:
      runAnalysis<NrbAnalysis>(argc, argv, analysisTypeOptions);
      break;

    case AnalysisType::NrbTrees:
      runAnalysis<NrbTrees>(argc, argv, analysisTypeOptions);
      break;

    case AnalysisType::ZGammaTrees:
      runAnalysis<ZGammaTrees>(argc, argv, analysisTypeOptions);
      break;

    case AnalysisType::ElectronTrees:
      runAnalysis<ElectronTrees>(argc, argv, analysisTypeOptions);
      break;

    case AnalysisType::EGammaFromMisid:
      runAnalysis<EGammaFromMisid>(argc, argv, analysisTypeOptions);
      break;
  }
}
