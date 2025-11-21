#ifndef HZZ2L2NU_INCLUDE_TRILEPTONTREES_H_
#define HZZ2L2NU_INCLUDE_TRLEPTONTREES_H_

#include <array>
#include <optional>
#include <tuple>
#include <vector>

#include <boost/program_options.hpp>
#include <TFile.h>
#include <TTree.h>
#include <TTreeReaderValue.h>

#include <EventTrees.h>
#include <Dataset.h>
#include <GenZZBuilder.h>
#include <Options.h>
#include <TriggerFilter.h>


/**
 * \brief Implements an analysis in the dilepton channel
 *
 * This class applies a slightly looser version of the standard analysis
 * selection in the dilepton channel: the leptons are not required to be of the
 * same flavour, and the cut on ptmiss is loosened. For each selected event, a
 * few observables, together with the event weight, are stored in a ROOT tree.
 * In addition, momenta and other properties of jets and leptons can be stored
 * if flag --more-vars is provided.
 *
 * Event weights are saved with the help of the base class EventTrees.
 */
class TrileptonTrees final : public EventTrees {
 public:
  TrileptonTrees(Options const &options, Dataset &dataset);

  /// Constructs descriptions for command line options
  static boost::program_options::options_description OptionsDescription();

  /// Performs the event selection and fills the output tree
  bool ProcessEvent();

 private:
  enum class LeptonCat : int {
    kEEE = 0,
    kMuMuMu = 1,
    kEEMu = 2,
    kMUMUE = 3
  };

  enum class JetCat : int {
    kEq0J,
    kEq1J,
    kGEq2J
  };

  /**
   * \brief Performs selection on leptons
   *
   * \return If the current event does not pass the selection, the returned
   * optional is empty. Otherwise it contains a tuple of the determined lepton
   * category and pointers to the two leptons, which are ordered by pt.
   */
  std::optional<std::tuple<LeptonCat, Lepton const *, Lepton const *, Lepton const *>>
  CheckLeptons() const;

  /// Fills additional variables, mostly lepton and jet momenta
  void FillMoreVariables(std::array<Lepton, 3> const &leptons,
      std::vector<Jet> const &jets);

  static double constexpr kNominalMZ_ = 91.1876;

  /// Indicates that additional variables should be stored
  bool storeMoreVariables_;

  /// Specifies the cut on ptmiss. Default is 80.
  double ptMissCut_;
  
  // Specific for Trilepton tree
  double mTwCut_ = 60;

  mutable std::unique_ptr<TTreeReaderValue<Float_t>> srcLHEVpt_;

  std::optional<Float_t> datasetLHEVptUpperLimitInc_;

  /**
   * \brief An object to reconstruct generator-level ZZ system
   *
   * Only created for datasets for ZZ production with decays to 2l2nu.
   */
  std::optional<GenZZBuilder> genZZBuilder_;

  TriggerFilter triggerFilter_;

  TTreeReaderValue<UInt_t> srcRun_;
  TTreeReaderValue<UInt_t> srcLumi_;
  TTreeReaderValue<ULong64_t> srcEvent_;

  Int_t leptonCat_, jetCat_, numPVGood_;
  Float_t llPt_, llEta_, llPhi_, llMass_;
  Float_t l3Pt_, l3Eta_, l3Phi_, l3Mass_;
  Float_t missPt_, missPhi_, missSignificance_, missSignificanceCorrected_;
  Float_t mTw_;
  Float_t dPhiVisiblesPtmiss_;
  Float_t mT_;

  TTreeReaderValue<int> srcNumPVGood_;

  UInt_t run_, lumi_;
  ULong64_t event_;
  Float_t genMZZ_;
  Int_t leptonCharge_[3];
  Float_t leptonPt_[3], leptonEta_[3], leptonPhi_[3], leptonMass_[3];
  static int const maxSize_ = 32;
  Int_t jetSize_;
  Float_t jetPt_[maxSize_], jetEta_[maxSize_], jetPhi_[maxSize_],
          jetMass_[maxSize_];

  Float_t dijetMass_;
};

#endif  // HZZ2L2NU_INCLUDE_TRILEPTONTREES_H_
