#include <yaml-cpp/yaml.h>

YAML::Node config = YAML::LoadFile("config_2017.yaml");

std::string nvtx_weight_filepath = config["nvtx_weight_filepath"].as<std::string>();
std::string eta_weight_filepath = config["eta_weight_filepath"].as<std::string>();
std::string pt_weight_filepath = config["pt_weight_filepath"].as<std::string>();
TFile* nvtx_weight_file = TFile::Open(nvtx_weight_filepath.c_str());
TFile* eta_weight_file = TFile::Open(eta_weight_filepath.c_str());
TFile* pt_weight_file = TFile::Open(pt_weight_filepath.c_str());


Float_t GetContentFromTH1F(Float_t x, TH1F const * hist) {
  return hist->GetBinContent(hist->FindFixBin(x));
}

Float_t GetContentFromTH1F(Int_t x, TH1F const * hist) {
  return hist->GetBinContent(hist->FindFixBin(x));
}


Float_t GetNvtxWeight(Int_t nvtx, Int_t jet_cat) {

  R__ASSERT(nvtx_weight_file && "Error opening ROOT file.");

  Float_t result;

  if (jet_cat == 0) {
    auto const hist = (TH1F*) nvtx_weight_file->Get("WeightHisto_eq0jets_ll");
    result = GetContentFromTH1F(nvtx, hist);
  } else if (jet_cat == 1) {
    auto const hist = (TH1F*) nvtx_weight_file->Get("WeightHisto_eq1jets_ll");
    result = GetContentFromTH1F(nvtx, hist);
  } else {
    auto const hist = (TH1F*) nvtx_weight_file->Get("WeightHisto_geq2jets_ll");
    result = GetContentFromTH1F(nvtx, hist);
  }

  return result;
}

Float_t GetEtaWeight(Float_t eta, Int_t jet_cat) {

  R__ASSERT(eta_weight_file && "Error opening ROOT file.");

  Float_t result;

  if (jet_cat == 2) {
    Float_t const abs_eta = std::abs(eta);
    auto const hist = (TH1F*) eta_weight_file->Get("WeightHisto_geq2jets_ll");
    result = GetContentFromTH1F(abs_eta, hist);
  } else {
    result = 1.;
  }

  return result;
}

Float_t GetPtWeight(Float_t pt, Int_t jet_cat) {

  R__ASSERT(pt_weight_file && "Error opening ROOT file.");

  Float_t result;

  if (jet_cat == 0) {
    auto const hist = (TH1F*) pt_weight_file->Get("WeightHisto_eq0jets_ll");
    result = GetContentFromTH1F(pt, hist);
  } else if (jet_cat == 1) {
    auto const hist = (TH1F*) pt_weight_file->Get("WeightHisto_eq1jets_ll");
    result = GetContentFromTH1F(pt, hist);
  } else {
    auto const hist = (TH1F*) pt_weight_file->Get("WeightHisto_geq2jets_ll");
    result = GetContentFromTH1F(pt, hist);
  }

  return result;
}

void InstrMETReweighting() {
  std::cout << GetNvtxWeight(11, 1) << std::endl;
  std::cout << GetNvtxWeight(26, 0) << std::endl;

  std::cout << std::endl;
  std::cout << GetEtaWeight(0., 1) << std::endl;
  std::cout << GetEtaWeight(1., 1) << std::endl;
  std::cout << GetEtaWeight(2., 1) << std::endl;
  std::cout << GetEtaWeight(0., 2) << std::endl;
  std::cout << GetEtaWeight(0.5, 2) << std::endl;
  std::cout << GetEtaWeight(-1., 2) << std::endl;
  std::cout << GetEtaWeight(1.25, 2) << std::endl;
  std::cout << GetEtaWeight(2., 2) << std::endl;

  std::cout << std::endl;
  std::cout << GetPtWeight(55., 2) << std::endl;
  std::cout << GetPtWeight(100., 2) << std::endl;
}