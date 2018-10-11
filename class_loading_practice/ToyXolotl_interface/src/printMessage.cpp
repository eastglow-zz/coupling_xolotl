#include "printMyStr.h"                                                       
#include <string>                                                             
                                                                             
class printMessage : public printMyStr {                                            
public:                                                                      
    virtual std::string message() {                                            
        return "This message was brought from the external dylib\n";                    
    }                                                                        
};                                                                           
                                                                             
                                                                             
// the class factories                                                       
                                                                             
extern "C" printMyStr* create() {                                               
    return new printMessage;                                                     
}                                                                            
                                                                             
extern "C" void destroy(printMyStr* p) {                                        
    delete p;                                                                
}           
