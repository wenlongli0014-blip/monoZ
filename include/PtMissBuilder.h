#ifndef PTMISSBUILDER_H_
#define PTMISSBUILDER_H_

#include <initializer_list>
#include <algorithm>
#include <optional>
#include <string>
#include <utility>
#include <vector>

#include <TTreeReaderArray.h>

#include <CollectionBuilder.h>
#include <Dataset.h>
#include <EventCache.h>
#include <Options.h>
#include <PhysicsObjects.h>


/**
 * \brief Lazily builds reconstructed missing pt
 *
 * Starts from raw missing pt. Type 1 correction, i.e. changes caused by
 * corrections applied to other objects, such as jets, can be included using
 * method \ref PullCalibration. Variations in "unclustered" momentum are applied
 * if requested via the \c syst option.
 *
 * If the master configuration contains field \c fix_ee_2017 in section
 * \c ptmiss and it is set to true, applies the
 * <a href="https://twiki.cern.ch/twiki/bin/viewauth/CMS/MissingETUncertaintyPrescription?rev=91#Instructions_for_2017_data_with">EE noise mitigation</a>.
 */
class PtMissBuilder {
 public:
  PtMissBuilder(Dataset &dataset, Options const &options);

  /// Returns missing pt in the current event
  PtMiss const &Get() const;

  /**
   * \brief Requests that changes in the total momentum produced by the given
   * builders are propagated into missing pt
   *
   * For example, changes in jet momenta due to variations in their corrections
   * can be included this way. If this method is called multiple times, new
   * builders are added to the list of already registered builders.
   */
  void PullCalibration(
    std::initializer_list<CollectionBuilderBase const *> builders);

 private:
  /// Supported systematic variations
  enum class Syst {
    None,
    UnclEnergyUp,
    UnclEnergyDown
  };

  /// Constructs ptmiss in the current event
  void Build() const;

  /// Systematic variation to be applied
  Syst syst_;

  /// Whether to apply the EE noise mitigation
  bool applyEeNoiseMitigation_;

  /// An object to facilitate caching
  EventCache cache_;

  /// Objects registered with \ref PullCalibration
  std::vector<CollectionBuilderBase const *> calibratingBuilders_;

  /// Object representing ptmiss in the current event
  mutable PtMiss ptMiss_;

  bool isSim_, applyXYCorrections_;

  mutable TTreeReaderValue<int> srcNumPV_;
  mutable TTreeReaderValue<float> srcPt_, srcPhi_;
  mutable TTreeReaderValue<UInt_t> srcRun_;
  mutable TTreeReaderValue<float> srcCOVxx_, srcCOVxy_, srcCOVyy_;
  mutable std::optional<TTreeReaderValue<float>> srcSignificance_;
  mutable std::optional<TTreeReaderValue<float>> srcUnclEnergyUpDeltaX_,
      srcUnclEnergyUpDeltaY_;

  /// Year for XY corrections. Needs to be 2016, 2016APV, 2016nonAPV, 2017 or 2018 to apply corr.
  std::string metXYCorrectionYear_;

  bool isUL_;

  /**
   * \brief Type 1 corrected default ptmiss and ptmiss with EE noise mitigation
   *
   * Used when \ref applyEeNoiseMitigation_ is set.
   */
  mutable std::optional<TTreeReaderValue<float>> srcDefaultPt_, srcDefaultPhi_,
      srcFixedPt_, srcFixedPhi_;
 };

#endif  // PTMISSBUILDER_H_
