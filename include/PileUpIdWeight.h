#ifndef HZZ2L2NU_INCLUDE_PILEUPIDWEIGHT_H_
#define HZZ2L2NU_INCLUDE_PILEUPIDWEIGHT_H_

#include <WeightBase.h>

#include <array>
#include <memory>
#include <optional>
#include <string_view>

#include <TH2.h>
#include <TTreeReaderValue.h>

#include <Dataset.h>
#include <EventCache.h>
#include <JetBuilder.h>
#include <Options.h>
#include <PhysicsObjects.h>
#include <PileUpIdFilter.h>
// #include <XGBoostPredictor.h>


/**
 * \brief Computes event weights that account for pileup ID scale factors
 *
 * A multivariate parameterization is used for the efficiency of the pileup ID
 * in simulation. Computed weights are cached on per-event basis.
 *
 * The behaviour is controlled by the section \c pileup_id of the master
 * configuration. Some information from it is extracted via PileUpIdFilter --
 * see its documentation on how the jet selection should be specified. In class
 * PileUpIdWeight it is assumed that section \c pileup_id exists; otherwise an
 * object of this class should not be constructed.
 *
 * Systematic variations in the scale factors are provided. The scale factors
 * for pileup and matched jets are varied simulatneously within each category
 * and independently between the categories.
 */
class PileUpIdWeight : public WeightBase {
 public:
  /**
   * \brief Constructor
   *
   * \param[in] dataset  Dataset that will be processed.
   * \param[in] options  Configuration options.
   * \param[in] pileUpIdFilter  Pointer to a valid PileUpIdFilter, which
   *   describes jet selection based on pileup ID.
   * \param[in] jetBuilder  Pointer to a valid JetBuilder.
   */
  PileUpIdWeight(
      Dataset &dataset, Options const &options,
      PileUpIdFilter const *pileUpIdFilter,
      JetBuilder const *jetBuilder);
  double NominalWeight() const override {
    if (cache_.IsUpdated())
      Update();
    return weights_[0];
  }

  // int NumVariations() const override {
  //   return 4;
  // }

  double operator()() const override {
    if (cache_.IsUpdated())
      Update();
    return weights_[int(defaultVariation_)];
  }

  // double RelWeight(int variation) const override {
  //   if (cache_.IsUpdated())
  //     Update();
  //   return weights_[variation + 1] / weights_[0];
  // }

  // std::string_view VariationName(int variation) const override;

 private:
  /**
   * \brief Supported systematic variations
   *
   * Cast to int, a variable of this type can index array \ref weights_.
   */
  enum class Variation : int {
    kNominal = 0,
    // kTagUp = 1,
    // kTagDown = 2,
    // kMistagUp = 3,
    // kMistagDown = 4
  };

  /// Auxiliary structure that provides details specific to each |eta| range
  struct Context {
    Context(Jet::PileUpId wp)
        : workingPoint{wp} {}

    /// Working point for pileup ID used in this range
    Jet::PileUpId workingPoint;

    /**
     * \brief Histograms with pileup ID scale factors and their uncertainties
     * for matched and pileup jets
     */
    std::shared_ptr<TH2> sf, eff, sf_unc;
  };

  /**
   * \brief Finds which |eta| range the jet falls into and returns the
   * corresponding context object
   */
  Context const &FindContext(Jet const &jet) const;

  /// Computes pileup ID efficiency in simulation for given jet
  double GetEfficiency(Context const &context, Jet const &jet) const;

  /// Finds pileup ID scale factor for given jet and variation
  double GetScaleFactor(Context const &context, Jet const &jet) const;

  /**
   * \brief Reads histograms with scale factors for all used working points
   *
   * \param[in] config  Master configuration.
   * \param[in] year  Year of data taking the scale factors for which should be
   *   read.
   */
  void LoadScaleFactors(YAML::Node const config, std::string year, bool isUL);

  /// Computes event weights for all systematic variations
  void Update() const;

  PileUpIdFilter const *pileUpIdFilter_;
  JetBuilder const *jetBuilder_;
  /**
   * \brief Edges between different |eta| regions
   *
   * Copied from pileUpIdFilter_.
   */
  std::vector<double> absEtaEdges_;

  /// \ref Context for each |eta| region
  std::vector<Context> contexts_;

  // /**
  //  * \brief Array with input features for the computation of the pileup ID
  //  * efficiency in simulation
  //  *
  //  * It will be reused for each computation. Features that don't change are set
  //  * in the constructor.
  //  */
  // mutable std::array<float, 13> effFeatures_;

  // /// Model that parameterizes the pileup ID efficiency in simulation
  // std::optional<XGBoostPredictor> effCalc_;

  /// Interface to read the expected number of pileup interactions
  mutable TTreeReaderValue<float> expPileUp_;

  /// Requested systematic variation
  Variation defaultVariation_;

  EventCache cache_;
  std::string systLabel_;
  /**
   * \brief Cached weights for all systematic variations
   *
   * The order is the same as in enum \ref Variation.
   */
  // mutable std::array<double, 5> weights_;
  mutable std::array<double, 1> weights_;
};

#endif  // HZZ2L2NU_INCLUDE_PILEUPIDWEIGHT_H_

