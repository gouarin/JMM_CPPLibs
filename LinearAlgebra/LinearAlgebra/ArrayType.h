//
//  ArrayType.h
//  LiftedFastMarching
//
//  Created by Jean-Marie Mirebeau on 14/09/2016.
//
//

#ifndef ArrayType_h
#define ArrayType_h

/*
 A very basic array type. The fast index is the first index.
 */

#include <vector>
#include "PointType.h"

namespace LinearAlgebra {

template<typename TComponent, size_t VDim, typename TDiscrete=int>
struct Array : public std::vector<TComponent> {
    typedef TComponent ComponentType;
    typedef TDiscrete DiscreteType;
    static const size_t Dimension = VDim;
    typedef Point<DiscreteType, Dimension> IndexType;
    
    typedef std::vector<ComponentType> Superclass;
    Array(Superclass && super_, IndexType dims_):Superclass(std::move(super_)),dims(dims_){};
    Array(const Superclass & super_,IndexType dims_):Superclass(super_),dims(dims_){};
    Array(){};
    
    IndexType dims;
    bool CheckDims() const {return dims.ProductOfCoordinates()==this->size();}
    bool InRange(const IndexType & index) const {
        for(size_t i=0; i<Dimension; ++i) if(index[i]<0 || index[i]>=dims[i]) return false;
        return true;}
    
    auto operator() (const IndexType & index) ->decltype((*this)[0]) {assert(CheckDims());
        return this->operator[](Convert(index));}
    auto operator() (const IndexType &) const ->decltype((*this)[0]) {assert(CheckDims());
        return this->operator[](Convert(index));}
    
    void PrintSelf(std::ostream & os) const;
    friend std::ostream & operator << (const Array & a, std::ostream & os){a.PrintSelf(os); return os;}

    DiscreteType Convert(const IndexType & index) const {
        assert(InRange(index));
        DiscreteType result = index[Dimension-1];
        for(int i=(int)Dimension-2; i>=0; --i)
            result = result*dims[i] + index[i];
        return result;
    }
    
    IndexType Convert(DiscreteType n) const {
        assert(0<=n && n<this->size());
        IndexType result;
        for(int i=0; i<Dimension-1; ++i){
            result[i]=n%dims[i];
            n/=dims[i];
        }
        result[Dimension-1]=n;
        return result;
    }
};
    
template<typename TC, size_t VD, typename TD> void
Array<TC,VD,TD>::PrintSelf(std::ostream & os) const {
    if(!CheckDims()) {
        os << "{\"Dimensions " << dims << "inconsistent with size " << this->size() << "\"}";
        return;
    }
    IndexType mods;
    mods[0]=dims[0];
    for(int i=1; i<Dimension;++i)
        mods[i]=dims[i]*mods[i-1];
    
    for(DiscreteType n=0; n<this->size();){
        for(int i=0; i<Dimension; ++i){
            if(n%mods[i]==0) os << "{";
            else break;
        }
        os << this->operator[](n);
        ++n;
        for(int i=0; i<Dimension; ++i){
            if(n%mods[i]==0) os << "}";
            else break;
        }
        if(n<this->size()) os << ",";
        else break;
    }
}

    
}

#endif /* ArrayType_h */
