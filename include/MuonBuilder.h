#ifndef MUONBUILDER_H_
#define MUONBUILDER_H_

#include <memory>
#include <optional>
#include <vector>

#include <TTreeReaderArray.h>

#include <CollectionBuilder.h>
#include <Dataset.h>
#include <Options.h>
#include <PhysicsObjects.h>
#include <RoccoR.h>
#include <TabulatedRandomGenerator.h>


/**
 * \brief Lazily builds collections of reconstructed muons
 *
 * For each event two collections of muons are constructed: tight and loose.
 * They differ in the minimal pt cut as well as identification requirements. The
 * tight collection is a subset of the loose one.
 *
 * Rochester corrections for muon momenta are applied. The changes in momenta of
 * loose muons are aggregated for \ref GetSumMomentumShift.
 */
class MuonBuilder : public CollectionBuilder<Muon> {
 public:
  /**
   * \brief Constructor
   *
   * \param[in] dataset    Dataset that will be processed.
   * \param[in] options    Configuration options for the job.
   * \param[in] rngEngine  Engine to construct TabulatedRandomGenerator.
   */
  MuonBuilder(Dataset &dataset, Options const &options,
              TabulatedRngEngine &rngEngine);

  /// Alias for \ref GetTight
  std::vector<Muon> const &Get() const override;

  /// Returns collection of loose muons
  std::vector<Muon> const &GetLoose() const;

  /// Returns collection of tight muons
  std::vector<Muon> const &GetTight() const;

 private:
  /**
   * \brief Applies Rochester correction to momentum of the muon
   *
   * \param[in] index     Index of the muon to choose channel for tabulatedRng_.
   * \param[in,out] muon  Muon to be corrected.
   * \param[in] trackerLayers  Number of tracker layers with measurements for
   *   the given muon.
   * 
   * The <a href="https://twiki.cern.ch/twiki/bin/view/CMS/RochcorMuon">
   * Rochester correction</a> is applied in place.
   */
  void ApplyRochesterCorrection(int index, Muon *muon, int trackerLayers) const;

  /// Constructs muons for the current event
  void Build() const override;
  
  /**
   * \brief Finds matching generator-level muon using (eta, phi) metric
   *
   * \param[in] muon  Recontructed muon for which a generator-level match needs
   *   to be found.
   * \param[in] maxDR  Maximal allowed distance in the (eta, phi) metric.
   * \return  An optional that contains the matching GenParicle. If no match is
   *   found in the specified cone, the optional is empty.
   */
  std::optional<GenParticle> FindGenMatch(Muon const &muon, double maxDR) const;

  /// Minimal pt for loose muons, GeV
  double minPtLoose_;

  /// Minimal pt for tight muons, GeV
  double minPtTight_;

  /// Maximal rel iso for loose muons
  double maxRelIsoLoose_;

  /// Maximal rel iso for tight muons
  double maxRelIsoTight_;

  /// Collection of muons passing loose selection
  mutable std::vector<Muon> looseMuons_;

  /// Collection of muons passing tight selection
  mutable std::vector<Muon> tightMuons_;

  /// Indicates whether running on simulation or data
  bool isSim_;

  std::string systLabel;

  /// Object to compute Rochester correction to muon pt
  std::unique_ptr<RoccoR> rochesterCorrection_;

  /// Random number generator
  TabulatedRandomGenerator tabulatedRng_;

  mutable TTreeReaderArray<float> srcPt_, srcEta_, srcPhi_, srcMass_;
  mutable TTreeReaderArray<int> srcCharge_;
  mutable TTreeReaderArray<float> srcIsolation_;
  mutable TTreeReaderArray<bool> srcIsPfMuon_, srcIsGlobalMuon_;
  mutable TTreeReaderArray<bool> srcIsTrackerMuon_;
  mutable TTreeReaderArray<bool> srcIdLoose_, srcIdTight_;
  mutable TTreeReaderArray<float> srcdxy_, srcdz_;
  mutable TTreeReaderArray<int> srcTrackerLayers_;
  mutable std::unique_ptr<TTreeReaderArray<int>> genPartId_;
  mutable std::unique_ptr<TTreeReaderArray<float>> genPartPt_, genPartEta_;
  mutable std::unique_ptr<TTreeReaderArray<float>> genPartPhi_;
};


inline std::vector<Muon> const &MuonBuilder::Get() const {
  return GetTight();
}

#endif  // MUONBUILDER_H_

