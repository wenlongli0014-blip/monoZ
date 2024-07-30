#include <yaml-cpp/yaml.h>

YAML::Node config = YAML::LoadFile("config_2016HIPM.yaml");

std::string fake_ratio_filepath = config["fake_ratio_filepath"].as<std::string>();
std::string photon_trigger_efficiency_filepath = config["photon_trigger_efficiency_filepath"].as<std::string>();
std::string electron_trigger_efficiency_filepath = config["electron_trigger_efficiency_filepath"].as<std::string>();

TFile* fake_ratio_file = TFile::Open(fake_ratio_filepath.c_str());
TFile* photon_trigger_efficiency_file = TFile::Open(photon_trigger_efficiency_filepath.c_str());
TFile* electron_trigger_efficiency_file = TFile::Open(electron_trigger_efficiency_filepath.c_str());



Float_t GetContentFromTH2F(Float_t x, Float_t y, TH2F const * hist) {
  return hist->GetBinContent(hist->FindFixBin(x, y));
}

Float_t GetContentFromTH1F(Float_t x, TH1F const * hist) {
  return hist->GetBinContent(hist->FindFixBin(x));
}

Float_t GetFakeRatio(Float_t pt, Float_t eta) {

  R__ASSERT(fake_ratio_file && "Error opening ROOT file.");

  Float_t result;

  if (std::abs(eta) <= 1.4442) {
    auto const hist = (TH1F*) fake_ratio_file->Get("ratio_barrel");
    R__ASSERT(hist && "Histogram not found in the ROOT file.");

    result = GetContentFromTH1F(pt, hist);
  } else if (std::abs(eta) >= 1.5660) {
    auto const hist = (TH1F*) fake_ratio_file->Get("ratio_endcap");
    R__ASSERT(hist && "Histogram not found in the ROOT file.");

    result = GetContentFromTH1F(pt, hist);
  } else {
    result = 0.;
  }

  return result;
}

Float_t GetPhotonTrigEff(Float_t pt, Float_t eta) {
  R__ASSERT(photon_trigger_efficiency_file && "Error opening ROOT file.");

  auto const hist = (TH2F*) photon_trigger_efficiency_file->Get("h2_efficiency");
  Float_t const max_pt = hist->GetXaxis()->GetBinCenter(hist->GetNbinsX());
  Float_t const photon_trigger_efficiency = GetContentFromTH2F(std::min(pt, max_pt), std::abs(eta), hist);

  return photon_trigger_efficiency;
}

Float_t GetElectronTrigEff(Float_t pt, Float_t eta) {
  R__ASSERT(electron_trigger_efficiency_file && "Error opening ROOT file.");

  auto const hist = (TH2F*) electron_trigger_efficiency_file->Get("h2_efficiency");
  Float_t const max_pt = hist->GetXaxis()->GetBinCenter(hist->GetNbinsX());
  Float_t const electron_trigger_efficiency = GetContentFromTH2F(std::min(pt, max_pt), std::abs(eta), hist);

  return electron_trigger_efficiency;
}

Float_t GetTrigEffRatio(Float_t pt, Float_t eta) {
  Float_t const result = GetPhotonTrigEff(pt, eta) / GetElectronTrigEff(pt, eta);
  return result;
}


void WJets() {
  std::cout << GetFakeRatio(50., 1.0) << std::endl;
  std::cout << GetFakeRatio(60., -1.0) << std::endl;
  std::cout << GetFakeRatio(100., 1.0) << std::endl;
  std::cout << GetFakeRatio(400., 1.0) << std::endl;
  std::cout << GetFakeRatio(500., 1.0) << std::endl;
  std::cout << GetFakeRatio(60., 2.0) << std::endl;
  std::cout << GetFakeRatio(100., 2.0) << std::endl;
  std::cout << GetTrigEffRatio(100., 1.5) << std::endl;

  std::cout << "\n" << std::endl;
  std::cout << GetTrigEffRatio(50., 1.0) << std::endl;
  std::cout << GetTrigEffRatio(60., -1.0) << std::endl;
  std::cout << GetTrigEffRatio(100., 1.0) << std::endl;
  std::cout << GetTrigEffRatio(400., 1.0) << std::endl;
  std::cout << GetTrigEffRatio(500., 1.0) << std::endl;
  std::cout << GetTrigEffRatio(60., 2.0) << std::endl;
  std::cout << GetTrigEffRatio(100., 2.0) << std::endl;
  std::cout << GetTrigEffRatio(100., 1.5) << std::endl;
  // std::cout << GetTrigEffRatio(100., 2.6) << std::endl;

}