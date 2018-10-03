
#include "XolotlExecutioner.h"

template<>
InputParameters validParams<XolotlExecutioner>()
{
  InputParameters params = validParams<Transient>();
  return params;
}

XolotlExecutioner::XolotlExecutioner(const InputParameters & parameters) :
    Transient(parameters)
{
}

// This method is called only one time at the start of the entire simulation.
void XolotlExecutioner::init() {
  Transient::init();
}
