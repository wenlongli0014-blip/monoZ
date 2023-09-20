#ifndef HZZ2L2NU_INCLUDE_EGAMMAFROMMISID_H_
#define HZZ2L2NU_INCLUDE_EGAMMAFROMMISID_H_

#include <array>
#include <optional>
#include <tuple>
#include <vector>
#include <variant>

#include <boost/program_options.hpp>
#include <TFile.h>
#include <TTree.h>
#include <TTreeReaderValue.h>

#include <Dataset.h>
// #include <EventNumberFilter.h>
#include <EventTrees.h>
#include <GenPhotonBuilder.h>
#include <Options.h>
#include <PhotonBuilder.h>
#include <PhotonWeight.h>
// #include <TriggerFilter.h>


class EGammaFromMisid final : public EventTrees {
 public:
  EGammaFromMisid(Options const &options, Dataset &dataset);

  /// Constructs descriptions for command line options
  static boost::program_options::options_description OptionsDescription();

  /// Performs the event selection and fills the output tree
  bool ProcessEvent();

 private:
  enum class EventCat : int {
    kEE = 0,
    kEGamma = 1
  };


  /**
   * \brief Performs selection on leptons
   *
   * \return If the current event does not pass the selection, the returned
   * optional is empty. Otherwise it contains a tuple of the determined lepton
   * category and pointers to the two leptons, which are ordered by pt.
   */
  // std::optional<std::tuple<EventCat, Lepton const *, Lepton const *>>
  // CheckElectrons() const;

  bool CheckProbe(std::variant<Electron const *, Photon const *>);

  /// Fills additional variables, mostly lepton and jet momenta
  void FillMoreVariables();

  static double constexpr kNominalMZ_ = 91.1876;

  /// Indicates that additional variables should be stored
  bool storeMoreVariables_;

  // TriggerFilter triggerFilter_;

  std::optional<GenPhotonBuilder> genPhotonBuilder_;

  PhotonBuilder photonBuilder_;

  PhotonWeight photonWeight_;

  // EventNumberFilter photonFilter_;

  std::optional<Int_t> datasetMinPtG_;
  std::optional<Int_t> datasetMaxPtG_;
  std::optional<Float_t> datasetLHEVptUpperLimitInc_;

  TTreeReaderValue<UInt_t> srcRun_;
  TTreeReaderValue<UInt_t> srcLumi_;
  TTreeReaderValue<ULong64_t> srcEvent_;

  mutable std::unique_ptr<TTreeReaderValue<Float_t>> srcLHEVpt_;
  mutable std::unique_ptr<TTreeReaderValue<UInt_t>> numGenPart_;
  mutable std::unique_ptr<TTreeReaderArray<Int_t>> genPartPdgId_;
  mutable std::unique_ptr<TTreeReaderArray<Float_t>> genPartPt_;
  mutable std::unique_ptr<TTreeReaderArray<Int_t>> genPartStatus_;
  mutable std::unique_ptr<TTreeReaderArray<Int_t>> genPartStatusFlags_;

  std::string datasetName_ = "";

  Int_t eventCat_;
  Int_t numPVGood_;
  Float_t probePt_, probeEta_, probePhi_, probeMass_;
  Float_t probeR9_;
  Float_t tagPt_, tagEta_, tagPhi_, tagMass_;
  Float_t totPt_, totEta_, totPhi_, totMass_;
  Float_t missPt_, missPhi_;
  // Float_t mT_;

  TTreeReaderValue<int> srcNumPVGood_;

  UInt_t run_, lumi_;
  ULong64_t event_;
  // Int_t leptonCharge_[2];
  // Float_t leptonPt_[2], leptonEta_[2], leptonPhi_[2], leptonMass_[2];
  // static int const maxSize_ = 32;
  Int_t jetSize_;
  // Float_t jetPt_[maxSize_], jetEta_[maxSize_], jetPhi_[maxSize_],
  //         jetMass_[maxSize_];

  bool isWJetsToLNu_;
};

#endif  // HZZ2L2NU_INCLUDE_EGAMMAFROMMISID_H_
