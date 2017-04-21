// Copyright 2017 Jean-Marie Mirebeau, University Paris-Sud, CNRS, University Paris-Saclay
// Distributed WITHOUT ANY WARRANTY. Licensed under the Apache License, Version 2.0, see http://www.apache.org/licenses/LICENSE-2.0

#ifndef FileIO_h
#define FileIO_h

#include <fstream>
#include <sstream>
#include "BaseIO.h"
/* 
 A simple interface for importing/exporting strings and numeric arrays (with double type) to files.
 
 Mostly used as an interface with Mathematica(R) and Python.
 */

        
struct FileIO : BaseIO {
    FileIO(std::string inputPrefix, std::string outputPrefix_);
    ~FileIO();
    
    /*
    bool HasField(KeyCRef) const;
    std::string GetString(KeyCRef) const;
    void SetString(KeyCRef, std::string);
    */
protected:
    const std::string outputPrefix;

//    template<typename T> std::pair<std::vector<DiscreteType>,const T*> GetDimsPtr(KeyCRef) const;

    /*
    std::vector<ScalarType> inputData;
    
    std::ostringstream outputFormat;
    std::vector<ScalarType> outputData;
     
    
    const InputFormatElement & GetInputFormat(KeyCRef) const;
    void SetExported(KeyCRef);
     */
};
    
#include "FileIO.hxx"

#endif /* FileIO_h */
