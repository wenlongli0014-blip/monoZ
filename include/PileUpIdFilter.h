#ifndef HZZ2L2NU_INCLUDE_PILEUPIDFILTER_H_
#define HZZ2L2NU_INCLUDE_PILEUPIDFILTER_H_

#include <utility>
#include <vector>

#include <Options.h>
#include <PhysicsObjects.h>


/**
 * \brief Applies pileup ID selection to jets
 *
 * The selection is specified in the section \c pileup_id of the master
 * configuration. If this section is missing, the filter is disabled and accepts
 * all jets. The following fields are read from this section:
 * - \c pt_range Jets whose pt is outside of this range are accepted
 *   automatically.
 * - \c abs_eta_edges Edges between bins in |eta|, for which different working
 *   points of pileup ID are to be applied. The left edge of the first bin and
 *   the right edge of the last bin are not included (e.g. if there are two bins
 *   in total, the sequence must include a single number). If the same working
 *   point is to be applied regardless of |eta|, the sequence must be empty. In
 *   this case it's optional and may be omitted from the configuration.
 * - \c working_points Working points for pileup ID to be used in each |eta|
 *   bin. Must be encoded with one-letter strings "N" (for no selection), "L",
 *   "M", or "T" (for loose, medium, or tight working points).
 */
class PileUpIdFilter {
 public:
  PileUpIdFilter(Options const &options);

  /// Return edges between bins in |eta| read from the configuration
  std::vector<double> const &GetAbsEtaEdges() const {
    return absEtaEdges_;
  }

  /// Returns range of pt where pileup ID is applicable
  std::pair<double, double> GetPtRange() const {
    return {minPt_, maxPt_};
  }

  /// Return working points read from the configuration
  std::vector<Jet::PileUpId> const &GetWorkingPoints() const {
    return workingPoints_;
  }

  /// Indicates whether pileup ID is applicable for given jet
  bool IsTaggable(Jet const &jet) const {
    return jet.p4.Pt() >= minPt_ and jet.p4.Pt() < maxPt_;
  }

  /**
   * \brief Checks whether the given jet passes the pileup ID selection
   * specified in the configuration
   *
   * If this filter is disabled or if the jet's pt falls outside of the range of
   * validity for pileup ID, the jet is automatically accepted.
   */
  bool operator()(Jet const &jet) const;

 private:
  /// Indicates whether the filter is enabled or works as a trivial pass-through
  bool enabled_;

  /// Range in pt where pileup ID is applicable
  double minPt_, maxPt_;

  /**
   * \brief Edges between diffent bins in |eta|
   *
   * Values are sorted in a strictly increasing order. The vector is empty if a
   * single working point is used for all pseudorapidities.
   */
  std::vector<double> absEtaEdges_;

  /**
   * \brief Working points for all bins in |eta|
   *
   * The size of this vector is exactly one unit larger than the size of
   * absEtaEdges_. The elements represent the loosest values of Jet::pileUpId
   * for the jet to be accepted.
   */
  std::vector<Jet::PileUpId> workingPoints_;
};

#endif  // HZZ2L2NU_INCLUDE_PILEUPIDFILTER_H_

