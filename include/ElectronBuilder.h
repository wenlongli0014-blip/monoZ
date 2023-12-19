#ifndef ELECTRONBUILDER_H_
#define ELECTRONBUILDER_H_

#include <vector>

#include <TTreeReaderArray.h>

#include <CollectionBuilder.h>
#include <Dataset.h>
#include <PhysicsObjects.h>
#include <Options.h>


/**
 * \brief Lazily builds collections of reconstructed electrons
 *
 * For each event two collections of electrons are constructed: tight and loose.
 * They differ in the minimal pt cut as well as identification requirements. The
 * tight collection is a subset of the loose one.
 *
 * Residual scale corrections in momenta of loose electrons are aggregated for
 * GetSumMomentumShift.
 */
class ElectronBuilder : public CollectionBuilder<Electron> {
 public:
  /// Constructor
  ElectronBuilder(Dataset &dataset, Options const &options);

  /// Alias for \ref GetTight
  std::vector<Electron> const &Get() const override;

  /// Returns collection of loose electrons
  std::vector<Electron> const &GetLoose() const;

  /// Returns collection of tight electrons
  std::vector<Electron> const &GetTight() const;

 private:
  // enum class Syst {
  //   None,
  // };

  /// Constructs electrons for the current event
  void Build() const override;

  // /// Systematic variation to be applied
  // Syst syst_;

  std::string systLabel;

  /// Minimal pt for loose electrons, GeV
  double minPtLoose_;

  /// Minimal pt for tight electrons, GeV
  double minPtTight_;

  /// Collection of electrons passing loose selection
  mutable std::vector<Electron> looseElectrons_;

  /// Collection of electrons passing tight selection
  mutable std::vector<Electron> tightElectrons_;

  mutable TTreeReaderArray<float> srcPt_, srcEta_, srcPhi_, srcMass_, srcDeltaEtaSc_;
  // mutable TTreeReaderArray<float> srcIsolation_;
  mutable TTreeReaderArray<int> srcCharge_;
  mutable TTreeReaderArray<bool> srcIdLoose_, srcIdTight_;
  mutable TTreeReaderArray<float> srcECorr_;

  mutable TTreeReaderArray<float> srcEnergyErr_;
};


inline std::vector<Electron> const &ElectronBuilder::Get() const {
  return GetTight();
}

#endif  // ELECTRONBUILDER_H_

