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

  // int NumVariations() const override {
  //   return 2;
  // }

  double operator()() const override {
    // std::cout << systLabel << std::endl;
    if (not systLabel.empty()) {
      // std::cout << RelWeight(systLabel) << std::endl;
      return NominalWeight() * RelWeight(systLabel);
    } else
      return NominalWeight();
  }

  double RelWeight(std::string const systLabel) const;

 private:
  
  //   void Update() const;
  
  std::string systLabel;

  mutable TTreeReaderArray<float> srcPSWeights_;
};

#endif  // HZZ2L2NU_INCLUDE_PSWEIGHT_H_
