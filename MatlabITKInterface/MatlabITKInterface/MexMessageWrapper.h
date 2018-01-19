// Copyright 2017 Jean-Marie Mirebeau, University Paris-Sud, CNRS, University Paris-Saclay
// Distributed WITHOUT ANY WARRANTY. Licensed under the Apache License, Version 2.0, see http://www.apache.org/licenses/LICENSE-2.0

#ifndef MatlabITKFastMarching_MexMessageWrapper_h
#define MatlabITKFastMarching_MexMessageWrapper_h


#include <sstream>
#include "mex.h"

template<int MessageGrade=0>
struct MexMessageWrapper {
    enum MexMessageGrade {Normal,Warning,Error};
    
    template<typename DataType>
    MexMessageWrapper & operator << (DataType data){oss << data; return (*this);}
    
    ~MexMessageWrapper(){
        std::string str = oss.str();
        switch (MessageGrade) {
            case Warning:
                mexWarnMsgTxt(str.c_str());
                break;
                
            case Error:
                mexErrMsgTxt(str.c_str());
                break;
                
            default:
            { // Solve an issue for displaying the % sign.
                size_t start_pos = 0; const std::string from="%", to ="%%";
                while((start_pos = str.find(from, start_pos)) != std::string::npos) {
                    str.replace(start_pos, from.length(), to);
                    start_pos += to.length(); // Handles case where 'to' is a substring of 'from'
                }
            }
                printf(str.c_str());
                break;
        };
        mexEvalString("drawnow");
        mexEvalString("pause(0.0001)");
    };
    
    std::ostringstream oss;
};
        
typedef MexMessageWrapper<MexMessageWrapper<>::Warning> MexWarnMsg;
typedef MexMessageWrapper<MexMessageWrapper<>::Normal>  MexMsg;
        
        
        

#endif
