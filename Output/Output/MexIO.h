//
//  MexIO.h
//  LiftedFastMarching_MatlabBuildTest
//
//  Created by Jean-Marie Mirebeau on 17/09/2016.
//
//

#ifndef MexIO_h
#define MexIO_h

#include "mex.h"
#include <sstream>

#define __IgnoreTrailingSingletonDimensions 1
#include "IO.h"

namespace Mex {
    
    struct BaseIO : TraitsIO {
        typedef double ScalarType;
        
        template<bool warn> struct Msg_ {
            std::ostringstream oss;
            template<typename T> Msg_ & operator << (const T & t){oss << t; return (*this);}
            ~Msg_();
        };
        
        bool HasField(KeyCRef) const;
        std::string GetString(KeyCRef) const;
        void SetString(KeyCRef, std::string);
        
        BaseIO(const mxArray *, mxArray **);
        BaseIO(const BaseIO &) = delete;
    protected:
        template<typename T> std::pair<std::vector<DiscreteType>,const T*> GetDimsPtr(KeyCRef) const;
        template<typename T, size_t d, typename F> void Set(KeyCRef, DimType<d>, const F &);
        mutable std::set<KeyType> unused, defaulted, exported;
        void SetField(KeyCRef, mxArray *);
        
        const mxArray * mxInput;
        mxArray ** mxOutput;
    };
    
#include "MexIO.hxx"
}


#endif /* MexIO_h */
