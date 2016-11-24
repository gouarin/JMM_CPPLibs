//
//  FileIO.h
//  LiftedFastMarching
//
//  Created by Jean-Marie Mirebeau on 14/09/2016.
//
//

#ifndef FileIO_h
#define FileIO_h

#include <map>
#include <fstream>
#include <sstream>
#include "IO.h"
/* 
 A simple interface for importing/exporting strings and numeric arrays (with double type) to files.
 
 Mostly used as an interface with Mathematica(R) and Pyhton.
 */

namespace File {
        
struct BaseIO : virtual TraitsIO {
    typedef double ScalarType;

    template<bool warn> struct Msg_ {
        Msg_(){assert(output!=nullptr); if(warn) (*output) << "----- ! Warning ! -----\n";}
        template<typename T> Msg_ & operator<<(const T & t){(*output) << t; return *this;}
        ~Msg_(){if(warn) (*output) << "-----------------------\n";}
        static std::ostream * output;
    };
    
    BaseIO(std::string inputPrefix, std::string outputPrefix_);
    BaseIO(const BaseIO &) = delete;
    ~BaseIO();
    
    bool HasField(KeyCRef) const;
    std::string GetString(KeyCRef) const;
    void SetString(KeyCRef, std::string);
    
protected:
    template<typename T> std::pair<std::vector<DiscreteType>,const T*> GetDimsPtr(KeyCRef) const;
    template<typename T, size_t d, typename F> void Set(KeyCRef, DimType<d>, const F &);
    mutable std::set<KeyType> unused, defaulted, exported;

    struct InputFormatElement;
    std::map<KeyType,InputFormatElement> inputFormat;
    std::vector<ScalarType> inputData;
    
    std::ostringstream outputFormat;
    std::vector<ScalarType> outputData;
    const std::string outputPrefix;
    
    const InputFormatElement & GetInputFormat(KeyCRef) const;
    void SetExported(KeyCRef);
};
    
    template<bool warn> std::ostream * BaseIO::Msg_<warn>::output = &std::cout;
    
#include "FileIO.hxx"

}
#endif /* FileIO_h */
