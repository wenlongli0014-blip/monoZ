#ifndef TAUBUILDER_H_
#define TAUBUILDER_H_

#include <initializer_list>
#include <vector>

#include <TTreeReaderArray.h>
#include <TTreeReaderValue.h>

#include <CollectionBuilder.h>
#include <Dataset.h>
#include <Options.h>
#include <PhysicsObjects.h>


/**
 * \brief Lazily builds a collection of reconstructed Taus
 *
 * The Tau veto is very effective to reject taus
 *
 */
class TauBuilder : public CollectionBuilder<Tau> {
 public:
  TauBuilder(Dataset &dataset, Options const &options);

  /// Returns collection of Taus
  std::vector<Tau> const &Get() const override;

 private:
  /// Constructs Taus for the current event
  void Build() const override;

  /// Minimal pt for Taus to select, GeV
  double minLepPt_;

  /// Collection of Taus
  mutable std::vector<Tau> Taus_;

  mutable TTreeReaderArray<float> srcPt_, srcEta_, srcPhi_;
  mutable TTreeReaderArray<int> srcDecayMode_;
  mutable TTreeReaderArray<UChar_t> srcIdDeepTau2017v2p1VSe_, srcIdDeepTau2017v2p1VSmu_, srcIdDeepTau2017v2p1VSjet_;

};

#endif  // TAUBUILDER_H_


