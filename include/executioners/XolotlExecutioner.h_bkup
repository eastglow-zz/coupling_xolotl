
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
 private:
   std::string _ext_lib_path_name;
   std::string _xolotl_input_path_name;
};
#endif //XOLOTLEXECUTIONER_H
