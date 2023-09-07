#include <PhotonWeight.h>

#include <sstream>
#include <stdexcept>

#include <TFile.h>

#include <FileInPath.h>
#include <Logger.h>

PhotonWeight::PhotonWeight(Dataset &, Options const &options,
                           PhotonBuilder const *photonBuilder)
    : photonBuilder_{photonBuilder} {
  for (auto const &photonPath : 
    Options::NodeAs<std::vector<std::string>>(
      options.GetConfig(), {"photon_efficiency", "photon"})) {
    photonTable_.emplace_back(ReadHistogram(photonPath));
  }
  LOG_WARN << "Photon trigger scale factors are missing";
}

double PhotonWeight::NominalWeight() const {
  double eff = 1.;
  for (auto &photon : photonBuilder_->Get()) {
    eff *= PhotonSF(photon);
  }
  return eff;
}

double PhotonWeight::PhotonSF(Photon const &photon) const {
  double eff = 1.;
  
  double pt = photon.p4.Pt();
  double eta = photon.p4.Eta();  // Use eta instead of etaSC
    
  for (auto &photonTable : photonTable_) {

    int photonEtaBin;
    if (photonTable->GetYaxis()->GetXmin() >= 0) {
        eta = fabs(eta);
    }
    photonEtaBin = photonTable->GetXaxis()->FindFixBin(eta);
    int photonPtBin = photonTable->GetYaxis()->FindFixBin(pt);
    if (pt > photonTable->GetYaxis()->GetXmax()) 
      photonPtBin = photonTable->GetNbinsY();

    eff *= photonTable->GetBinContent(photonEtaBin, photonPtBin);
  }
  return eff;
}

std::unique_ptr<TH2> PhotonWeight::ReadHistogram(
    std::string const &pathWithName) {
  auto const pos = pathWithName.find_last_of(':');

  if (pos == std::string::npos)
    throw HZZException{"Histogram path does not contain ':'."};

  std::filesystem::path path = FileInPath::Resolve(pathWithName.substr(0, pos));
  std::string name = pathWithName.substr(pos + 1);

  TFile inputFile{path.c_str()};

  if (inputFile.IsZombie()) {
    HZZException exception;
    exception << "Could not open file \"" << path << "\".";
    throw exception;
  }

  std::unique_ptr<TH2> hist;
  hist.reset(dynamic_cast<TH2 *>(inputFile.Get(name.c_str())));

  if (not hist) {
    HZZException exception;
    exception << "File " << path << " does not contain required histogram \""
        << name << "\".";
    throw exception;
  }

  hist->SetDirectory(nullptr);
  inputFile.Close();

  return hist;
}
