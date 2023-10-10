#include <TauBuilder.h>

#include <cmath>
#include <cstdlib>
#include <algorithm>

#include <Utils.h>


TauBuilder::TauBuilder(Dataset &dataset, Options const &)
    : CollectionBuilder{dataset.Reader()},
      minLepPt_{20.},
      srcPt_{dataset.Reader(), "Tau_pt"},
      srcEta_{dataset.Reader(), "Tau_eta"},
      srcPhi_{dataset.Reader(), "Tau_phi"},
      srcDecayMode_{dataset.Reader(), "Tau_decayMode"},
      srcIdDeepTau2017v2p1VSe_{dataset.Reader(), "Tau_idDeepTau2017v2p1VSe"},
      srcIdDeepTau2017v2p1VSmu_{dataset.Reader(), "Tau_idDeepTau2017v2p1VSmu"},
      srcIdDeepTau2017v2p1VSjet_{dataset.Reader(), "Tau_idDeepTau2017v2p1VSjet"}
{}


std::vector<Tau> const &TauBuilder::Get() const {
  Update();
  return Taus_;
}


void TauBuilder::Build() const {
  Taus_.clear();

  for (unsigned i = 0; i < srcPt_.GetSize(); ++i) {
    if (srcDecayMode_[i] == 5 or srcDecayMode_[i] == 6) continue;  // decay modes 5 and 6 are rejected

    if (not (srcIdDeepTau2017v2p1VSe_[i] >= 2)) continue;  // VVLoose
    if (not (srcIdDeepTau2017v2p1VSmu_[i] >= 1)) continue;  // VVVLoose
    if (not (srcIdDeepTau2017v2p1VSjet_[i] >= 16)) continue;  // Medium

    float pt = srcPt_[i];
    if (not (pt > minLepPt_ )) continue; // lower pt treshold for taus

    if (not (std::abs(srcEta_[i]) < 2.3)) continue;

    Tau tau;
    tau.p4.SetPtEtaPhiM(srcPt_[i], srcEta_[i], srcPhi_[i], 0.);

    // Perform angular cleaning
    if (IsDuplicate(tau.p4, 0.4))
      continue;

    Taus_.emplace_back(tau);
  }

  // Make sure the collection is ordered in pt
  std::sort(Taus_.begin(), Taus_.end(), PtOrdered);
}
