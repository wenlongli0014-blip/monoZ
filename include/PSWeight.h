#ifndef HZZ2L2NU_INCLUDE_PSWEIGHT_H_
#define HZZ2L2NU_INCLUDE_PSWEIGHT_H_

#include <WeightBase.h>

#include <filesystem>
#include <memory>
#include <string>

#include <TH2.h>

#include <Dataset.h>
#include <Options.h>
#include <PhotonBuilder.h>
#include <PhysicsObjects.h>


class PSWeight : public WeightBase {
 public:
  /// Constructor
  PSWeight(Dataset &dataset, Options const &options);

  virtual double NominalWeight() const override;

  int NumVariations() const override {
    return 0;
  }

//  private:

//   void Update() const;

};

#endif  // HZZ2L2NU_INCLUDE_PSWEIGHT_H_
