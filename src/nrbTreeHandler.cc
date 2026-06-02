#include <cstdlib>
#include <iostream>
#include <vector>
#include <boost/algorithm/string.hpp>
#include <boost/program_options.hpp>

#include <Options.h>

#include <ROOT/RDataFrame.hxx>
#include <TBranch.h>
#include <TFile.h>
#include <TH1.h>
#include <TTree.h>
#include <TString.h>


namespace po = boost::program_options;


void AlphaHandler(std::string filename, std::string outputname) {
  ROOT::RDataFrame d("Vars", filename);
  auto N_ll = d.Filter("((ll_mass < 75. && ll_mass > 50.) || "
      "(ll_mass > 105 && ll_mass < 200)) && (btag_tight || "
      "(btag_tight_lowpt && jet_cat == 0)) && ll_pt > 55. && ptmiss > 80. "
      "&& lepton_cat != 2").Count();
  auto N_emu = d.Filter("((ll_mass < 75. && ll_mass > 50.) || "
      "(ll_mass > 105 && ll_mass < 200)) && (btag_tight || "
      "(btag_tight_lowpt && jet_cat == 0)) && ll_pt > 55. && ptmiss > 80. "
      "&& lepton_cat == 2").Count();
  double alpha = static_cast<double>( N_ll.GetValue())
    / static_cast<double>(N_emu.GetValue());
  std::cout << N_ll.GetValue() << " / " << N_emu.GetValue() << "== ";
  std::cout << "alpha:" << alpha << std::endl;
  TFile *oldFile = new TFile(filename.c_str());
  TTree *oldTree;
  oldFile->GetObject("Vars", oldTree);
  TFile *newFile = new TFile(outputname.c_str(), "recreate");
  oldTree->SetBranchStatus("tmp*", 0);
  auto newtree = oldTree->CloneTree(0);
  newtree->Branch("weight", &alpha, "weight/F");
  for (int i = 0; i < oldTree->GetEntries(); ++i)
    newtree->Fill();
  auto srtree = newtree->CopyTree("lepton_cat == 2 && "
      "ll_mass >76.2 && ll_mass < 106.2 && btag_loose == 0 && ll_pt > 55. ");
  srtree->Write();
  newFile->Close();
  oldFile->Close();
  delete newFile;
  delete oldFile;
}

void PtmissReweight(std::string filename, std::string outputname) {
  ROOT::RDataFrame d("Vars",filename);
  float ptmiss_binning[] = {50, 125, 140, 190, 250, 350, 1000};
  ROOT::RDF::TH1DModel histModel("","",
      sizeof(ptmiss_binning) / sizeof(float) - 1, ptmiss_binning);
  std::string baseCut_ll("((ll_mass < 76.2 && ll_mass > 50.) || "
      "(ll_mass > 106.2 && ll_mass < 201.2)) " 
      "&& ( lepton_cat == 0 || lepton_cat == 1)");
  std::string baseCut_emu("((ll_mass < 76.2 && ll_mass > 50.) || "
      "(ll_mass > 106.2 && ll_mass < 201.2))  && lepton_cat == 2");
  std::vector<std::string> extraCuts;
  extraCuts.emplace_back(" && jet_cat == 0 && btag_loose_lowpt ");
  extraCuts.emplace_back(" && jet_cat == 1 && btag_tight ");
  extraCuts.emplace_back(" && jet_cat == 2 && btag_tight ");
  std::vector<TH1F*> hists;
  hists.resize(3);	
  for (long unsigned int i = 0; i < extraCuts.size() ;i++) {
    auto d_ll = d.Filter(baseCut_ll + extraCuts[i]);
    auto d_emu = d.Filter(baseCut_emu + extraCuts[i]);
    auto proxy_ll = (TH1F *)(d_ll.Histo1D(histModel, 
        "ptmiss").GetValue().Clone());
    auto proxy_emu = (TH1F *)(d_emu.Histo1D(histModel,
        "ptmiss", "tmp_weight").GetValue().Clone());
    hists[i] = (TH1F *)(proxy_ll->Clone());
    hists[i]->Divide(proxy_emu);
  }
  TFile *oldFile = new TFile(filename.c_str());
  TTree *oldTree;
  oldFile->GetObject("Vars", oldTree);
  auto *brlist = oldTree->GetListOfBranches();
  std::vector<TString> tmpWeightBranches, weightBranches;
  for (int i = 0 ; i < brlist->GetEntries(); i++ ){ 
    TString name(((TBranch*)brlist->At(i))->GetName());
    if (name.Contains("tmp")) {
      tmpWeightBranches.push_back(name);
      weightBranches.push_back(name.ReplaceAll(TString("tmp_"),TString("")));
    }
  }
  std::vector<float> oldWeights, neWeights;
  float corr_up, corr_down, ptmiss;
	float f_corr_up, f_corr_down, f_corr;
  int jet_cat;
  oldWeights.resize(weightBranches.size());
  neWeights.resize(weightBranches.size());
  TFile *newFile = new TFile(outputname.c_str(), "recreate");
  auto newtree = oldTree->CloneTree(0);
  for (long unsigned int i = 0; i < weightBranches.size(); i++){
    oldTree->SetBranchAddress(tmpWeightBranches[i], &oldWeights[i]);
    newtree->Branch(weightBranches[i], &neWeights[i], weightBranches[i]+"/F");
  }
  oldTree->SetBranchAddress("ptmiss", &ptmiss);
  oldTree->SetBranchAddress("jet_cat", &jet_cat);
  newtree->Branch("weight_sidebandEff_up", &corr_up);
  newtree->Branch("weight_sidebandEff_down", &corr_down);
  newtree->Branch("fcorr", &f_corr_up);
  newtree->Branch("fcorr_up", &f_corr_up);
  newtree->Branch("fcorr_down", &f_corr_down);
  int Num = oldTree->GetEntries();
  for (int i = 0; i < Num ; ++i){
    oldTree->GetEntry(i);
    int bin = hists[jet_cat]->FindFixBin(ptmiss);
    f_corr = hists[jet_cat]->GetBinContent(bin);
    f_corr_up = f_corr + hists[jet_cat]->GetBinErrorUp(bin);
    f_corr_down = f_corr - hists[jet_cat]->GetBinErrorLow(bin);	
    float nominal = 1.;
    for (long unsigned int i = 0; i < weightBranches.size(); i++){
      if (tmpWeightBranches[i] == "tmp_weight")
        nominal = oldWeights[i];
      neWeights[i] = oldWeights[i] * f_corr;
    }
    corr_up = f_corr_up * nominal;
    corr_down = f_corr_down * nominal;
    newtree->Fill();
  }
  auto srtree = newtree->CopyTree("lepton_cat == 2 " 
      "&& ll_mass >76.2 && ll_mass < 106.2 && btag_loose == 0 && ll_pt > 55. ");
  srtree->Write();
  newFile->Close();
  oldFile->Close();
  delete newFile;
  delete oldFile;
}
void KMethodHandler(std::string filename, std::string outputname) {
  TFile *oldFile = new TFile(filename.c_str());
  TTree *oldTree;
  oldFile->GetObject("Vars", oldTree);
  auto *brlist = oldTree->GetListOfBranches();
  std::vector<TString> tmpWeightBranches, weightBranches;
  for (int i = 0 ; i < brlist->GetEntries(); i++ ){ 
    TString name(((TBranch*)brlist->At(i))->GetName());
    if (name.Contains("tmp")) {
      tmpWeightBranches.push_back(name);
      weightBranches.push_back(name.ReplaceAll(TString("tmp_"),TString("")));
    }
  }
  std::vector<float> oldWeights, neWeights;
  oldWeights.resize(weightBranches.size());
  TFile *newFile = new TFile(outputname.c_str(), "recreate");
  auto newtree = oldTree->CloneTree(0);
  for (long unsigned int i = 0; i < weightBranches.size(); i++){
    oldTree->SetBranchAddress(tmpWeightBranches[i], &oldWeights[i]);
    newtree->Branch(weightBranches[i], &oldWeights[i], weightBranches[i]+"/F");
  }
  int Num = oldTree->GetEntries();
  for (int i = 0; i < Num ; ++i){
    oldTree->GetEntry(i);
    newtree->Fill();
  }
  auto srtree = newtree->CopyTree("lepton_cat == 2 " 
      "&& ll_mass >76.2 && ll_mass < 106.2 && btag_loose == 0 && ll_pt > 55. ");
  srtree->Write();
  newFile->Close();
  oldFile->Close();
  delete newFile;
  delete oldFile;
}

int main(int argc, char **argv) {
  po::options_description analysisTypeOptions{"Analysis type"};
  analysisTypeOptions.add_options()
    ("input,i",po::value<std::string>()->default_value("Data.root"),
		 "input path: --input=/pathToFile or -i /pathToFile")
    ("output,o", po::value<std::string>()->default_value("NRB_weights.root"),
		 "output file name: --output=name.root or -o name.root")
    ("analysis,a", po::value<std::string>()->default_value("kmethod"),
     "Analysis to run; allowed values are \"kmethod\", "
     "\"alphamethod\",\"kcorrection\"");

  // Command line options are checked twice. At the first pass only check the
  // analysis type and update the list of expected options accordingly.
  po::variables_map vm;
  po::store(po::command_line_parser(argc, argv).options(analysisTypeOptions)\
      .allow_unregistered().run(), vm);
  po::notify(vm);

  std::string analysisTypeArg{vm["analysis"].as<std::string>()};
  std::string filename{vm["input"].as<std::string>()};
  std::string outputname{vm["output"].as<std::string>()};
  boost::to_lower(analysisTypeArg);

  if (analysisTypeArg == "kmethod")
    KMethodHandler(filename, outputname);
  else if (analysisTypeArg == "alphamethod")
    AlphaHandler(filename, outputname);
  else if (analysisTypeArg == "kcorrection")
    PtmissReweight(filename, outputname);
  else {
    LOG_ERROR << "Unknown analysis type \"" <<
      vm["analysis"].as<std::string>() << "\"";
    std::exit(EXIT_FAILURE);
  }
}
