#ifndef HZZ2L2NU_INCLUDE_JETBUILDER_H_
#define HZZ2L2NU_INCLUDE_JETBUILDER_H_

#include <initializer_list>
#include <optional>
#include <vector>

#include <TTreeReaderArray.h>
#include <TTreeReaderValue.h>

#include <CollectionBuilder.h>
#include <Dataset.h>
#include <GenJetBuilder.h>
#include <JetCorrector.h>
#include <Options.h>
#include <PhysicsObjects.h>
#include <PileUpIdFilter.h>


/**
 * \brief Lazily builds collection of reconstructed jets that pass analysis
 * selection
 *
 * Reads parameters from sections \c jets and \c ptmiss of the master
 * configuration. The latter specify how jets contribute to the type 1
 * correction of missing pt.
 *
 * Jets are checked against the kinematical selection specified in the
 * configuration, as well as jet ID. Selection on pileup ID is applied if the
 * corresponding filter object is provided to the constructor.
 *
 * Systematic variations in jets are implemented with the help of JetCorrector,
 * which also applies JER smearing. To follow the standard smearing algorithm,
 * this builder needs to be made aware of generator-level jets via method
 * \ref SetGenJetBuilder.
 *
 * Jet with corrected pt > 15 GeV (the exact meaning depends on the
 * configuration) are aggregated for \ref GetSumMomentumShift to be used for the
 * type 1 correction of missing pt. The details are controlled by section
 * \c ptmiss of the master configuration.  The following fields are supported:
 * - \c jer  Boolean indicating whether JER smearing should be propagated into
 *     missing pt.
 * - \c pog_jets  Optional boolean; defaults to false. Requests that jets with a
 *     high EM energy fraction are skipped from the computation and constituent
 *     muons are subtracted from the remaining jets.  This is how the standard
 *     type 1 correction is computed.
 * - \c fix_ee_2017 Optional boolean; defaults to false. Indicates whether the
 *   <a href="https://twiki.cern.ch/twiki/bin/viewauth/CMS/MissingETUncertaintyPrescription?rev=91#Instructions_for_2017_data_with">EE noise mitigation</a>
 *   should be applied.
 * Jets that geometrically overlap with other objects (as checked with
 * IsDuplicate) are never propagated into the missing pt.
 */
class JetBuilder : public CollectionBuilder<Jet> {
 public:
  JetBuilder(
      Dataset &dataset, Options const &options, TabulatedRngEngine &rngEngine,
      PileUpIdFilter const *pileUpIdFilter = nullptr);

  /// Returns collection of jets
  std::vector<Jet> const &Get() const override;
  
	std::vector<Jet> const &GetLowPt() const;
  /**
   * \brief Returns collection of jets that have been rejected by pileup ID but
   * satisfy other requirements.
   *
   * Filled only if a PileUpIdFilter has been provided.
   */
  std::vector<Jet> const &GetRejected() const;

  /**
   * \brief Specifies an object that provides generator-level jets
   *
   * Generator-level jets are used for JER smearing. However, if no such object
   * is given, stochastic version of the smearing will be performed. When given,
   * the object must have an approriate life time.
   */
  void SetGenJetBuilder(GenJetBuilder const *genJetBuilder);

 private:
  /**
   * \brief Adds a contribution to the type 1 correction of missing pt
   *
   * \param[in] rawP4  Raw four-momentum of the jet.
   * \param[in] area  Area of the jet.
   * \param[in] jecOrig  The full JEC applied during production of NanoAOD.
   * \param[in] jecNew  The full JEC applied in the analysis, including the JEC
   *   systematic variation.
   * \param[in] jerFactor  Scale factor representing JER smearing on top of
   *   the full JEC.
   * \param[in] emFraction  Total (neutral and charged) electromagnetic energy
   *   fraction in a jet. A non-positive value disables the corresponding cut.
   * \param[in] muonFraction  Fraction of (raw) jet energy carried by muons.
   *
   * Implements the type 1 correction according to the full - L1 scheme.
   * The full correction includes JEC systematic variations (so that they are
   * propagated into missing pt) and, depending on the value of \ref ptMissJer_,
   * JER smearing.
   *
   * If \ref ptMissPogJets_ is \c true, jets with a high EM energy fraction are
   * skipped and constituent muons are removed. This mimics the standard
   * procedure to compute the type 1 correction.
   *
   * If the EE noise mitigation is enabled, soft jets in the EE are excluded
   * from this procedure. At the same time these jets are used to undo some of
   * the contributions in the starting ptmiss.
   */
  void AddType1Correction(
      TLorentzVector const &rawP4, double area,
      double jecOrig, double jecNew, double jerFactor,
      double emFraction, double muonFraction) const;

  /**
   * \brief Constructs jets in the current event
   *
   * Actual work is delegated to ProcessJets and ProcessSoftJets.
   */
  void Build() const override;

  /**
   * \brief Finds matching generator-level jet using JERC definition
   *
   * Returns a nullptr if no match is found within the allowed cone.
   */
  GenJet const *FindGenMatch(TLorentzVector const &p4,
                             double ptResolution) const;

  /**
   * \brief Computes momentum scale factor that accounts for JER smearing
   *
   * \param[in] corrP4  Corrected four-momentum of the jet.
   * \param[in] rngChannel  Channel for the random number generator. Should be
   *   set to the index of the jet in the current event.
   */
  double GetJerFactor(TLorentzVector const &corrP4, int rngChannel) const;

  /// Constructs collection of jets in the current event
  void ProcessJets() const;

  /**
   * \brief Processes soft jets in the corrent event
   *
   * These are jets from branches <tt>CorrT1METJet_*</tt>.
   */
  void ProcessSoftJets() const;

  /**
   * \brief Deals with pileup ID for given jet
   *
   * Sets data members in \c jet related to pileup ID and applies the filtering
   * on pileup ID if it is enabled.
   *
   * \param[in,out] jet  Jet to be checked and updated.
   * \param[in] index  Index of the jet in the source branches.
   * \return  Boolean indicating of the jet passes the filtering on pileup ID.
   *   If the filtering is disabled, always returns true.
   */
  bool SetPileUpInfo(Jet &jet, int index) const;

  /**
   * \brief Non-owning pointer to an object that produces generator-level jets
   *
   * May be nullptr.
   */
  GenJetBuilder const *genJetBuilder_;

  /**
   * \brief Non-owning pointer to an object that implements jet filtering based
   * on pileup ID
   *
   * May be nullptr.
   */
  PileUpIdFilter const *pileUpIdFilter_;

  /// Minimal pt for jets, GeV
  double minPt_;

  /// Maximal |eta| for jets
  double maxAbsEta_;

  /**
   * \brief Minimal (corrected) pt for jets to be used in type 1 correction of
   * missing pt
   */
  double minPtType1Corr_;

  /// Whether to apply the EE noise mitigation in missing pt
  bool ptMissEeNoise_;

  /// Whether to propagate JER smearing into missing pt
  bool ptMissJer_;

  /**
   * \brief Whether to compute the type 1 correction for missing pt skipping
   * jets with a high EM energy fraction and removing leptons from them
   *
   * This is how the standard type 1 correction is computed.
   */
  bool ptMissPogJets_;

  /// Range of pt, in GeV, where pileup ID is applicable
  double pileUpIdMinPt_, pileUpIdMaxPt_;
  int pileUpIdLooseBit_, pileUpIdMediumBit_, pileUpIdTightBit_;
  /// Collection of jets
  mutable std::vector<Jet> jets_;

  mutable std::vector<Jet> lowptJets_;
  /// Collection of jets rejected by pileup ID
  mutable std::vector<Jet> rejectedJets_;

  /// Indicates whether running on simulation or data
  bool isSim_;

  /// Index of the bit with jet ID decision
  int jetIdBit_;

  /// Object that computes JEC
  JetCorrector jetCorrector_;

  mutable TTreeReaderArray<float> srcPt_, srcEta_, srcPhi_, srcMass_;
  mutable TTreeReaderArray<float> srcArea_, srcRawFactor_;
  mutable TTreeReaderArray<float> srcChEmEF_, srcNeEmEF_, srcMuonFraction_;
  mutable TTreeReaderArray<float> srcBTag_;
  mutable TTreeReaderArray<int> srcId_, srcPileUpId_;
  mutable TTreeReaderValue<float> puRho_;
  mutable std::optional<TTreeReaderArray<int>> srcHadronFlavour_,
      srcPartonFlavour_, srcGenJetIdx_;

  // Properties of soft jets, which are used in the type 1 correction of ptmiss
  mutable TTreeReaderArray<float> softRawPt_, softEta_, softPhi_, softArea_;
  mutable TTreeReaderArray<float> softMuonFraction_;
};

#endif  // HZZ2L2NU_INCLUDE_JETBUILDER_H_

