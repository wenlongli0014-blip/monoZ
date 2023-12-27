#include <PSWeight.h>

#include <sstream>
#include <stdexcept>

#include <Logger.h>

PSWeight::PSWeight(Dataset &, Options const &) {}

double PSWeight::NominalWeight() const {
  return 1.;
}
