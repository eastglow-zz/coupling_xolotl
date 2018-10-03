
#ifndef XOLOTLEXECUTIONER_H
#define XOLOTLEXECUTIONER_H

#include "Transient.h"

class XolotlExecutioner;

template<>
InputParameters validParams<XolotlExecutioner>();

class XolotlExecutioner : public Transient
{
 public:
   XolotlExecutioner(const InputParameters & parameters);

   virtual void init();
};
#endif //XOLOTLEXECUTIONER_H
