//
//  PointBaseType.h
//  AmongConvex2
//
//  Created by Jean-Marie Mirebeau on 24/07/13.
//  Copyright (c) 2013 Jean-Marie Mirebeau. All rights reserved.
//

#ifndef AmongConvex2_PointBaseType_h
#define AmongConvex2_PointBaseType_h

#include <iostream>
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

    PointBase(std::initializer_list<TComponent> c){
        assert(c.size()==Dimension);
        std::copy(c.begin(), c.end(), this->begin());
    }
    
    PointBase(){};
    
    bool AreAllCoordinatesNonNegative() const;
    bool AreAllCoordinatesPositive() const;
    static bool IsInRange(int pos) {return 0<=pos && pos<Dimension;}
    bool IsFinite() const;
    
    ComponentType SumOfCoordinates() const;
    ComponentType ProductOfCoordinates() const;
    static PointBase Constant(ComponentType c){PointBase p; p.fill(c); return p;}
    
    struct LexicographicCompare;
    
protected:
    template<typename Component2>
    static PointBase CastCoordinates(const PointBase<Component2,Dimension> & u){
        PointBase v;
        for(size_t i=0; i<Dimension; ++i) v[i] = ComponentType(u[i]);
        return v;
    }    
};


    
template<typename TC, size_t VD>
std::ostream & operator << (std::ostream & f, const PointBase<TC,VD> & p)
{
    f<<"{";
    for(int i=0; i<VD; ++i){
        if(i>0) f<<",";
        f<<p[i];
    }
    f<<"}";
    return f;
}
    
template<typename TC, size_t VD>
bool
PointBase<TC,VD>::AreAllCoordinatesNonNegative() const {
    for(int i=0; i<Dimension; ++i)
        if(this->operator[](i)<0) return false;
    return true;
}

template<typename TC, size_t VD>
bool
PointBase<TC,VD>::AreAllCoordinatesPositive() const {
    for(int i=0; i<Dimension; ++i)
        if(this->operator[](i)<=0) return false;
    return true;
}

template<typename TC, size_t VD>
bool
PointBase<TC,VD>::IsFinite() const {
    for(int i=0; i<Dimension; ++i)
        if(!(std::fabs(this->operator[](i))<std::numeric_limits<ComponentType>::infinity()))
            return false;
    return true;
}
    
template<typename TC, size_t VD>
typename PointBase<TC,VD>::ComponentType
PointBase<TC,VD>::SumOfCoordinates() const {
    ComponentType sum(0);
    for(int i=0; i<Dimension; ++i)
        sum += this->operator[](i);
    return sum;
}

template<typename TC, size_t VD>
typename PointBase<TC,VD>::ComponentType
PointBase<TC,VD>::ProductOfCoordinates() const {
    ComponentType prod(1);
    for(int i=0; i<Dimension; ++i)
        prod *= this->operator[](i);
    return prod;
}

    template<typename TC, size_t VD>
    struct PointBase<TC, VD>::LexicographicCompare {
        bool operator()(const PointBase & p, const PointBase & q) const {
            for(auto pi=p.begin(), qi=q.begin(); pi!=p.end(); ++pi, ++qi){
                if(*pi < *qi) return true;
                if(*pi > *qi) return false;
            }
            return false; // strict ordering
        }
    };
    
    
}

#endif
