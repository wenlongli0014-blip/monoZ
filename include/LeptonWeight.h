#ifndef HZZ2L2NU_INCLUDE_LEPTONWEIGHT_H_
#define HZZ2L2NU_INCLUDE_LEPTONWEIGHT_H_

#include <WeightBase.h>

#include <Dataset.h>
#include <ElectronBuilder.h>
#include <MuonBuilder.h>
#include <Options.h>
#include <PhysicsObjects.h>


class PtEtaHistogram;


/**
 * \brief Applies lepton identification efficiency scale factors
 *
 * The computation is done using tight leptons provided by the builders given to
 * the constructor. The scale factors are read from files specified in the
 * master configuration. Multiple scale factors for the same lepton flavour are
 * multiplied together. See documentation for class PtEtaHistogram for a
 * description of the configuration for individual scale factors.
 *
 * Trigger scale factors are missing.
 *
 * Systematic variations are not provided.
 */
class LeptonWeight : public WeightBase {
 public:
  /// Constructor
  LeptonWeight(Dataset &dataset, Options const &options,
               ElectronBuilder const *electronBuilder,
               MuonBuilder const *muonBuilder, int efficiencyType = 0);

  ~LeptonWeight();

  /// Computes the total lepton efficiency weight for the current event

  double NominalWeight() const override {
    if (cache_.IsUpdated())
      Update();
    return weights_[0];
  }


  double operator()() const override {
    if (cache_.IsUpdated())
      Update();
    return weights_[defaultWeightIndex_];
  }
  
	double RelWeight(int variation) const override {

    if (cache_.IsUpdated())
      Update();
    return weights_[variation + 1] / weights_[defaultWeightIndex_];
	}
	
	int NumVariations() const override {
    return 4;
  }
  /// Computes single lepton efficiency with given pt & eta & variation

  double ElectronEff(double pt, double eta, int syst = 0) const;
  
  double MuonEff(double pt, double eta, int syst = 0) const;
	
  /// Computes transfer factor for emu reweighting
  double TransferFactor(Lepton const *l1, Lepton const *l2, int syst = 0) const;

	std::string_view VariationName(int variation) const override;
 
 private:
  /// Computes all weights for the current event
  void Update() const;
  
	/**
   * \brief Gives scale factor for an electron
   *
   * Uses pseudorapidity of the ECAL supercluster to look up the scale factor.
   */
  double ElectronSF(Electron const &electron, int syst) const;

  /**
   * \brief Gives scale factor for a muon
   *
   * Uses uncorrected momentum to look up the scale factor.
   */
  double MuonSF(Muon const &muon, int syst) const;
  
	/**
   * \brief Cached weights
   *
   * They are given in the order nominal, 
	 * muonEff_syst_up, muonEff_syst_down, 
	 * muonEff_stat_up, muonEff_stat_down, 
	 * muonEff_syst_up, muonEff_syst_down, 
	 * muonEff_stat_up, muonEff_stat_down, 
	 * electronEff_syst_up, electronEff_syst_down,
	 * electronEff_stat_up, electronEff_syst_down.
	 * electronEff_syst_up, electronEff_syst_down,
	 * electronEff_stat_up, electronEff_syst_down.
   */
  mutable std::array<double, 17> weights_;
  
	EventCache cache_;
  
	/**
   * \brief Index of the default weight to be returned by operator()
   *
   * Corresponds to array \ref weights_.
   */
  int defaultWeightIndex_;
  
	/// Individual scale factors for muons and electrons
  std::array<std::vector<PtEtaHistogram>, 9> 
		muonScaleFactors_, electronScaleFactors_;

  /// Non-owning pointer to object that provides collection of electrons
  ElectronBuilder const *electronBuilder_;

  /// Non-owning pointer to object that provides collection of muons
  MuonBuilder const *muonBuilder_;
};

#endif  // HZZ2L2NU_INCLUDE_LEPTONWEIGHT_H_

