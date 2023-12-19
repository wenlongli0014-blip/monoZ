#include <ElectronBuilder.h>

#include <cmath>
#include <cstdlib>
#include <algorithm>


ElectronBuilder::ElectronBuilder(Dataset &dataset, Options const &options)
    : CollectionBuilder{dataset.Reader()},
      minPtLoose_{10.}, minPtTight_{15.},
      // maxRelIsoLoose_{0.4}, maxRelIsoTight_{0.1},
      srcPt_{dataset.Reader(), "Electron_pt"},
      srcEta_{dataset.Reader(), "Electron_eta"},
      srcPhi_{dataset.Reader(), "Electron_phi"},
      srcMass_{dataset.Reader(), "Electron_mass"},
      srcDeltaEtaSc_{dataset.Reader(), "Electron_deltaEtaSC"},
      // srcIsolation_{dataset.Reader(), "Electron_pfRelIso03_all"},
      srcCharge_{dataset.Reader(), "Electron_charge"},
      srcIdLoose_{dataset.Reader(), "Electron_mvaFall17V2Iso_WPL"},
      srcIdTight_{dataset.Reader(), "Electron_mvaFall17V2Iso_WP90"},
      srcECorr_{dataset.Reader(), "Electron_eCorr"},
      srcEnergyErr_{dataset.Reader(), "Electron_energyErr"} {
      // syst_{Syst::None} {

      systLabel = options.GetAs<std::string>("syst");
        // if (systLabel == "electronEnergy_up")
        //   syst_ = Syst::UnclEnergyUp;
        // else if (systLabel == "electronEnergy_down")
        //   syst_ = Syst::UnclEnergyDown;
        // else
        //   syst_ = Syst::None;
      }


std::vector<Electron> const &ElectronBuilder::GetLoose() const {
  Update();
  return looseElectrons_;
}


std::vector<Electron> const &ElectronBuilder::GetTight() const {
  Update();
  return tightElectrons_;
}


void ElectronBuilder::Build() const {

  looseElectrons_.clear();
  tightElectrons_.clear();

  for (unsigned i = 0; i < srcPt_.GetSize(); ++i) {
    double const eta = srcEta_[i];
    double const etaSc = srcDeltaEtaSc_[i] + eta;
    double const absEtaSc = std::abs(etaSc);
    bool const passIdLoose = srcIdLoose_[i];
    bool const passIdTight = srcIdTight_[i];

    double pt = srcPt_[i];

    // std::cout << "nominal pt " << pt << std::endl;
    // std::cout << systLabel << std::endl;
    if (systLabel == "electronEnergy_up") {
      pt += srcEnergyErr_[i] / std::cosh(eta);
    }
    if ((systLabel == "electronEnergy_down")) {
      pt -= srcEnergyErr_[i] / std::cosh(eta);
    }
    // std::cout << "varied pt " << pt << std::endl;

    if (not passIdLoose)
      continue;

    if (not (pt > minPtLoose_))
      continue;

    if (not (absEtaSc < 2.5))
      continue;

    Electron electron;
    electron.p4.SetPtEtaPhiM(pt, srcEta_[i], srcPhi_[i], srcMass_[i]);
    electron.charge = srcCharge_[i];
    electron.etaSc = etaSc;

    if (IsDuplicate(electron.p4, 0.1))
      continue;

    looseElectrons_.emplace_back(electron);

    // Propagate corrections to momenta of loose electrons to ptmiss
    TLorentzVector const uncorrP4 = electron.p4 * (1. / srcECorr_[i]);
    AddMomentumShift(uncorrP4, electron.p4);

    if (not passIdTight)
      continue;

    if (not (pt > minPtTight_))
      continue;

    if (absEtaSc > 1.4442 and absEtaSc < 1.5660)  // EB-EE gap
      continue;

    tightElectrons_.emplace_back(electron);


  }

  // Make sure the collections are ordered in pt
  std::sort(looseElectrons_.begin(), looseElectrons_.end(), PtOrdered);
  std::sort(tightElectrons_.begin(), tightElectrons_.end(), PtOrdered);
}

