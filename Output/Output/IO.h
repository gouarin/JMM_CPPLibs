//
//  IO.h
//  ExternLFM
//
//  Created by Jean-Marie Mirebeau on 21/09/2016.
//
//

#ifndef IO_h
#define IO_h

/** A simple import/export interface. 
 Note that most of the complexity is due to Matlab(R)'s habits of tampering with trailing singleton dimensions, and permuting the first two coordinate axes in images.
 */

#include <string>
#include <set>
#include "Output/ExceptionMacro.h"
#include "LinearAlgebra/ArrayType.h"

#ifndef __IgnoreTrailingSingletonDimensions
#define __IgnoreTrailingSingletonDimensions 0
#endif

struct TraitsIO {
    template<typename TC, size_t VD> using Array = LinearAlgebra::Array<TC,VD>;
    typedef std::string KeyType;
    typedef std::string const & KeyCRef;
    typedef Array<double,1>::DiscreteType DiscreteType;
protected:
    template<size_t d> using DimType = LinearAlgebra::Point<DiscreteType,d>;
};

/*
 Base must provide the following
 
 public:
 template<bool warn> struct Msg_;
 bool HasField(KeyCRef) const;
 std::string GetString(KeyCRef) const;
 void SetString(KeyCRef, std::string);

 protected:
 template<typename T> std::pair<std::vector<DiscreteType>,const T*> GetDimsPtr(KeyCRef) const;
 template<typename T, size_t d, typename F> void Set(KeyCRef, DimType<d>, const F &);
 mutable std::set<KeyType> unused, defaulted;
 */

template<typename Base> struct IO_ : Base, TraitsIO {
    using Base::Base;
    IO_(const IO_ &) = delete;
    ~IO_();
    
    typedef typename Base::template Msg_<true> WarnMsg;
    typedef typename Base::template Msg_<false> Msg;
    using Base::GetString;
    std::string GetString(KeyCRef, const std::string &) const;
    
    template<typename T> T Get(KeyCRef) const;
    template<typename T> T Get(KeyCRef, const T &) const;
    template<typename T> std::vector<T> GetVector(KeyCRef key) const {return GetArray<T, 1>(key);}
    template<typename T, size_t d> Array<T, d> GetArray(KeyCRef) const;
    template<typename T> std::vector<DiscreteType> GetDimensions(KeyCRef) const;
    
    // Output
    template<typename T> void Set(KeyCRef, const T &);
    template<typename T> void SetVector(KeyCRef, const std::vector<T> &);
    template<typename T, size_t d> void SetArray(KeyCRef, const Array<T, d> &);
    
    int verbosity=1;  bool transposeFirstTwoCoordinates=false;
protected:
    template<typename V> static V TransposeDims(const V &);
    template<typename T, size_t d> struct TransposeVals;
    template<typename T, size_t d> void Set(KeyCRef, DimType<d>, const T*);
};

#include "IO.hxx"

#endif /* IO_h */
