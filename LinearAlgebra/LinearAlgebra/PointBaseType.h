// Copyright 2017 Jean-Marie Mirebeau, University Paris-Sud, CNRS, University Paris-Saclay
// Distributed WITHOUT ANY WARRANTY. Licensed under the Apache License, Version 2.0, see http://www.apache.org/licenses/LICENSE-2.0

#ifndef AmongConvex2_PointBaseType_h
#define AmongConvex2_PointBaseType_h

#include <iostream>
#include <limits>
#include <numeric>
#include <array>
#include <cassert>
#include <cmath>

namespace LinearAlgebra {
    
template<typename TComponent, size_t VDimension>
struct PointBase :
public std::array<TComponent,VDimension>
{
    typedef TComponent ComponentType;
    static const int Dimension = VDimension;
    typedef std::array<TComponent,VDimension> DataType;
    
    bool AreAllCoordinatesNonNegative() const;
    bool AreAllCoordinatesPositive() const;
    static bool IsInRange(int pos) {return 0<=pos && pos<Dimension;}
    bool IsFinite() const;
    
    ComponentType SumOfCoordinates() const;
    ComponentType ProductOfCoordinates() const;
    static PointBase Constant(ComponentType c){PointBase p; p.fill(c); return p;}
    
    struct LexicographicCompare;
/*    template<typename ...Components> PointBase(Components... comps) : DataType{std::forward<Components>(comps)...}{}*/
protected:
//    constexpr PointBase(const DataType & data):DataType(data){};
    template<typename ...T,typename dummy = typename std::enable_if<sizeof...(T)==Dimension>::type >
    constexpr PointBase(T... t):DataType{{ComponentType(t)...}}{};
    PointBase(){};
    
    template<typename Component2>
    static PointBase CastCoordinates(const PointBase<Component2,Dimension> & u){
        PointBase v;
        for(size_t i=0; i<Dimension; ++i){
            if (std::numeric_limits<ComponentType>::is_integer // constexpr
                && !std::numeric_limits<Component2>::is_integer
                && std::numeric_limits<Component2>::is_specialized)
                 v[i] = ComponentType(round(u[i]));
            else v[i] = ComponentType(u[i]);
        }
        return v;
    }    
};


    
template<typename TC, size_t VD>
std::ostream & operator << (std::ostream & f, const PointBase<TC,VD> & p)
{
    f<<"{";
    for(int i=0; i<VD; ++i){
        if(i>0) f<<",";
        f<< (std::is_same<TC, char>::value ? (int) p[i] : p[i]);
    }
    f<<"}";
    return f;
}
    
    
#include "Implementation/PointBaseType.hxx"
}

#endif
