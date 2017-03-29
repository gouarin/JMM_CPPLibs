// Copyright 2017 Jean-Marie Mirebeau, University Paris-Sud, CNRS, University Paris-Saclay
// Distributed WITHOUT ANY WARRANTY. Licensed under the Apache License, Version 2.0, see http://www.apache.org/licenses/LICENSE-2.0

#ifndef IO_h
#define IO_h

/** A simple import/export interface. 
 Note that most of the complexity is due to Matlab(R)'s habit of tampering with trailing singleton dimensions, and of permuting the first two coordinate axes in images.
 
 The code WILL throw an exception, with a descriptive error message, whenever the user input is invalid. Make sure that your user input code is exception safe.
 */

#include <string>
#include <set>
#include "Output/ExceptionMacro.h"
#include "LinearAlgebra/ArrayType.h"
#include "Output/EnumToString.h"


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

enum class ArrayOrdering {Default, Reversed, Transposed};
template<> char const * enumStrings<ArrayOrdering>::data[] = {
    "Default", "Reversed", "Transposed"};

template<typename Base> struct IO_ : Base, virtual TraitsIO {
    using Base::Base;
    IO_(const IO_ &) = delete;
    ~IO_();
    
    typedef typename Base::template Msg_<true> WarnMsg;
    typedef typename Base::template Msg_<false> Msg;
    using Base::GetString;
    std::string GetString(KeyCRef, const std::string &, int=1) const;
    
    template<typename T> T Get(KeyCRef) const;
    template<typename T> T Get(KeyCRef, const T &, int=1) const;
    template<typename T> std::vector<T> GetVector(KeyCRef key) const;
    template<typename T, size_t d> Array<T, d> GetArray(KeyCRef) const;
    template<typename T> std::vector<DiscreteType> GetDimensions(KeyCRef) const;
    
    // Output
    template<typename T> void Set(KeyCRef, const T &);
    template<typename T> void SetVector(KeyCRef, const std::vector<T> &);
    template<typename T, size_t d> void SetArray(KeyCRef, const Array<T, d> &);
    
    int verbosity=1;
    ArrayOrdering arrayOrdering = ArrayOrdering::Default;
protected:
    template<typename V> static V TransposeDims(const V &);
    template<typename V> static V ReverseDims(const V &);
    template<typename T, size_t d> struct TransposeVals;
    template<typename T, size_t d> struct ReverseVals;
    template<typename T, size_t d> void Set(KeyCRef, DimType<d>, const T*);
};

#include "IO.hxx"

#endif /* IO_h */
