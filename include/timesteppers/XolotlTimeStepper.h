#ifndef XOLOTLTIMESTEPPER_H
#define XOLOTLTIMESTEPPER_H

#include "TimeStepper.h"

class XolotlTimeStepper;

template<>
InputParameters validParams<XolotlTimeStepper>();

class XolotlTimeStepper : public TimeStepper
{
public:
  XolotlTimeStepper(const InputParameters & parameters);

protected:
  virtual Real computeInitialDT() override;
  virtual Real computeDT() override;
  virtual void step() override;
  virtual void postExecute() override;
  virtual bool converged() override;

private:
  Real _dt;
  std::string _ext_lib_path_name;
  std::string _func_name;
};

#endif //XOLOTLTIMESTEPPER_H
