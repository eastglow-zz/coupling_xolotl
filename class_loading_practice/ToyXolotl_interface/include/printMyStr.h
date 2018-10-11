#ifndef PRINTMYSTR_H                                                          
#define PRINTMYSTR_H   

#include <string>                                                       
                                                                             
class printMyStr {                                                              
public:                                                                      
    printMyStr() {
    }                                                 
    virtual ~printMyStr() {
    }                                                    
                                                                             
    virtual std::string  message() = 0; 
};                                                                           
                                                                             
// the types of the class factories                                          
//typedef printMyStr* create_t();                                                 
//typedef void destroy_t(printMyStr*);                                            
                                                                             
#endif                                                                       

