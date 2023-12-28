#include <PSWeight.h>

#include <sstream>
#include <stdexcept>

#include <Logger.h>

PSWeight::PSWeight(Dataset &dataset, Options const &options) 
  : srcPSWeights_{dataset.Reader(), "PSWeight"} {
    systLabel = options.GetAs<std::string>("syst");
  }

double PSWeight::NominalWeight() const {
  return 1.;
}

double PSWeight::RelWeight(std::string const systLabel) const {
  // std::cout << systLabel << std::endl;
  if (systLabel == "PS_isr_up")
    return srcPSWeights_[0];
  else if (systLabel == "PS_isr_down")
    return srcPSWeights_[2];
  else if (systLabel == "PS_fsr_up")
    return srcPSWeights_[1];
  else if (systLabel == "PS_fsr_down")
    return srcPSWeights_[3];
  else
    return 1.;
}