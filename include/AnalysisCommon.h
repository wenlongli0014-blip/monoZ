#ifndef HZZ2L2NU_INCLUDE_ANALYSISCOMMON_H_
#define HZZ2L2NU_INCLUDE_ANALYSISCOMMON_H_

#include <optional>

#include <boost/program_options.hpp>

#include <BTagger.h>
#include <BTagWeight.h>
#include <Dataset.h>
#include <ElectronBuilder.h>
#include <EWCorrectionWeight.h>
#include <GenJetBuilder.h>
#include <GenWeight.h>
#include <IsoTrackBuilder.h>
#include <JetBuilder.h>
#include <JetGeometricVeto.h>
#include <KFactorCorrection.h>
#include <L1TPrefiringWeight.h>
#include <LeptonWeight.h>
#include <MeKinFilter.h>
#include <MetFilters.h>
#include <MuonBuilder.h>
#include <Options.h>
#include <PileUpIdFilter.h>
#include <PileUpIdWeight.h>
#include <PileUpWeight.h>
#include <PSWeight.h>
#include <PtMissBuilder.h>
#include <RunSampler.h>
#include <TabulatedRandomGenerator.h>
#include <TauBuilder.h>
#include <TriggerWeight.h>
#include <WeightCollector.h>


/**
 * \brief Implements common aspects of analyses
 *
 * This class creates a number of blocks used in multiple analyses, thus
 * avoiding code duplication. Provided blocks include
 * - builders for leptons, jets, and missing pt,
 * - common event filters,
 * - common reweighting objects,
 * and some others. They are provided to the inheriting classes as protected
 * data members.
 */
class AnalysisCommon {
 public:
  AnalysisCommon(Options const &options, Dataset &dataset);

  /// Constructs descriptions for common command line options
  static boost::program_options::options_description OptionsDescription();

 protected:
  /**
   * \brief Checks if event passes common event-level filters
   *
   * Method ProcessEvent of a derived class should normally call this.
   */
  bool ApplyCommonFilters() const;

  /// Computes the absolute value of phi between ptmiss and the system of
  /// leptons (or photon, for the CR) and jets.
  double DPhiPtMiss(
    const std::initializer_list<CollectionBuilderBase const *> &builders);

  double DPhiPtMiss2(const TLorentzVector &p4Miss,
    const std::initializer_list<CollectionBuilderBase const *> &builders);
  /// Integrated luminosity, 1/pb
  double intLumi_;

  /// Selection cuts values
  // double zMassWindow_, minPtLL_, minDphiLLPtMiss_, minDphiJetPtMiss_,
        //  minDphiLeptonsJetsPtMiss_;
  double zMassWindow_, minPtLL_, minDphiLLPtMiss_, minDphiJetPtMiss_;

  /// Indicates whether this is simulation or real data
  bool isSim_;

  /// Common random number generator engine
  TabulatedRngEngine tabulatedRngEngine_;

  RunSampler runSampler_;

  BTagger bTagger_;
  PileUpIdFilter pileUpIdFilter_;

  ElectronBuilder electronBuilder_;
  MuonBuilder muonBuilder_;
  IsoTrackBuilder isotrkBuilder_;
  TauBuilder tauBuilder_;
  std::optional<GenJetBuilder> genJetBuilder_;
  JetBuilder jetBuilder_;
  PtMissBuilder ptMissBuilder_;

  LeptonWeight leptonWeight_;
  TriggerWeight triggerWeight_;
  std::optional<GenWeight> genWeight_;
  std::optional<PSWeight> psWeight_;
  std::optional<KFactorCorrection> kFactorCorrection_;
  std::optional<EWCorrectionWeight> ewCorrectionWeight_;
  std::optional<PileUpWeight> pileUpWeight_;
  std::optional<L1TPrefiringWeight> l1tPrefiringWeight_;
  std::optional<BTagWeight> bTagWeight_;
  std::optional<PileUpIdWeight> pileUpIdWeight_;
  WeightCollector weightCollector_;


 private:
  MeKinFilter meKinFilter_;
  MetFilters metFilters_;
  JetGeometricVeto jetGeometricVeto_;
};

#endif  // HZZ2L2NU_INCLUDE_ANALYSISCOMMON_H_
