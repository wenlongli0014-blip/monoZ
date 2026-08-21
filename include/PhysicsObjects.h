#ifndef HZZ2L2NU_INCLUDE_PHYSICSOBJECTS_H_
#define HZZ2L2NU_INCLUDE_PHYSICSOBJECTS_H_

#include <cmath>
#include <cstdlib>
#include <limits>

#include <TLorentzVector.h>


/// Non-polymorphic base class for particle-like physics objects
struct Particle {

  /// Four-momentum, in GeV
  TLorentzVector p4;
};


/**
 * \brief Comparison function to sort collections of particle-like objects in pt
 *
 * \param[in] p1, p2  Particle-like objects.
 * \return True if the pair of given particle-like objects is ordered in pt,
 *   meaning that pt of the first object is greater than pt of the second one.
 */
inline bool PtOrdered(Particle const &p1, Particle const &p2) {
  return (p1.p4.Pt() > p2.p4.Pt());
}


/// Generator-level particle
struct GenParticle : public Particle {

  /// Constructor from PDG ID
  GenParticle(int pdgId) noexcept;

  /// PDG ID code
  int pdgId;
};


inline GenParticle::GenParticle(int pdgId_) noexcept
    : Particle{}, pdgId{pdgId_} {}


/// Generator-level jet
struct GenJet : public Particle {};


/// Reconstructed jet
struct Jet : public Particle {
  /// Working points of pileup jet ID
  enum class PileUpId : int {
    None = 0,
    Loose,
    Medium,
    Tight,
    /// Used for high-pt jets for which pileup ID is not applicable
    PassThrough
  };

  /// Default constructor
  Jet() noexcept;

  /// Sets jet flavours
  void SetFlavours(int hadronFlavour, int partonFlavour);

  /**
   * \brief Value of CSVv2 b-tagging discriminator
   *
   * NaN if not set.
   */
  double bTag;

  /**
   * \brief Hadron flavour of the jet
   *
   * Possible values are 0, 4, and 5.
   */
  int hadronFlavour;

  /**
   * \brief Combination of hadron and parton flavours
   *
   * Equals to the hadron flavour if it's non-zero, otherwise the absolute value
   * of the parton flavour.
   */
  int combFlavour;

  /**
   * \brief Tightest working point of pileup jet ID passed by this jet
   *
   * PileUpId::PassThrough if not set.
   */
  PileUpId pileUpId;

  /**
   * \brief Boolean indicating whether the jet is a true pileup jet
   *
   * The definition should be consistent with what is used for pileup ID. The
   * value is always false in real data.
   */
  bool isPileUp;
};


inline Jet::Jet() noexcept
    : Particle{}, bTag{std::numeric_limits<double>::quiet_NaN()},
      hadronFlavour{0}, combFlavour{0},
      pileUpId{PileUpId::PassThrough}, isPileUp{false} {}


inline void Jet::SetFlavours(int hadronFlavour_, int partonFlavour) {
  hadronFlavour = hadronFlavour_;
  combFlavour =
      (hadronFlavour != 0) ? hadronFlavour : std::abs(partonFlavour);
}


/// Reconstructed missing pt
struct PtMiss : public Particle {

  /// Default constructor
  PtMiss() noexcept;

  /**
   * \brief Value of ptmiss significance
   *
   * NaN if not set.
   */
  double significance;
  double significance_corrected;
};


inline PtMiss::PtMiss() noexcept
    : Particle{}, significance{std::numeric_limits<double>::quiet_NaN()}, significance_corrected{std::numeric_limits<double>::quiet_NaN()} {}


/// Reconstructed photon
struct Photon : public Particle {

   /// Gen photon origin, for matched gen-level particle
  enum class Origin {
    PromptPhoton,
    PromptElectron,
    Unmatched
  };

  /// Photon flavour (prompt photon, prompt electron, unknown/unmatched)
  Origin flavour = Origin::Unmatched;

  /// Four-momentum of the matched gen-level particle, in GeV
  TLorentzVector genP4;

  /// Photon variables
  bool passElecVeto;
  bool isEB;
  float r9;
  float sieie;

  /// Default constructor
  Photon() noexcept;
};

inline Photon::Photon() noexcept
    : Particle{} {}


struct IsoTrack : public Particle {

  /// Default constructor
  IsoTrack() noexcept;
};

inline IsoTrack::IsoTrack() noexcept
    : Particle{} {}

/// Reconstructed charged lepton
struct Lepton : public Particle {

  /// Lepton flavour
  enum class Flavour {
    Electron,
    Muon,
    Tau
  };

  /// Constructor from flavour
  Lepton(Flavour flavour) noexcept;

  /// Flavour
  Flavour flavour;

  /**
   * \brief Electric charge
   *
   * Allowed values are +-1 and 0, the latter meaning that the charge has not
   * been specified.
   */
  int charge;

  /**
   * \brief Uncorrected momentum
   *
   * If momentum of the lepton is corrected, its old value is expected to be
   * stored in this field. This is useful to access lepton scale factors using
   * the old momentum. If the momentum has not been corrected, set to a null
   * vector.
   */
  TLorentzVector uncorrP4;
};


inline Lepton::Lepton(Flavour flavour_) noexcept
    : Particle{}, flavour{flavour_}, charge{0}, uncorrP4{} {}


/// Reconstructed electron
struct Electron : public Lepton {

  /// Default constructor
  Electron() noexcept;

  /**
   * \brief Pseudorapidity of associated ECAL supercluster
   *
   * NaN if not set.
   */
  double etaSc;
};


inline Electron::Electron() noexcept
    : Lepton{Lepton::Flavour::Electron},
      etaSc{std::numeric_limits<double>::quiet_NaN()} {}


/// Reconstructed muon
struct Muon : public Lepton {

  /// Default constructor
  Muon() noexcept;
};


inline Muon::Muon() noexcept
    : Lepton{Lepton::Flavour::Muon} {}

/// Reconstructed tau
struct Tau : public Lepton {

  /// Default constructor
  Tau() noexcept;
};


inline Tau::Tau() noexcept
    : Lepton{Lepton::Flavour::Tau} {}
#endif  // HZZ2L2NU_INCLUDE_PHYSICSOBJECTS_H_
