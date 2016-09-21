//
//  DifferentiationType.h
//  OptimalTransport
//
//  Created by Jean-Marie Mirebeau on 17/01/2014.
//  Copyright (c) 2014 Jean-Marie Mirebeau. All rights reserved.
//

#ifndef OptimalTransport_DifferentiationType_h
#define OptimalTransport_DifferentiationType_h

// A class for automatic differentiation

#include "FriendOperators.h"
#include "math.h"
#include <cassert>

template<typename TScalar, typename TVector>
struct DifferentiationType :
LinearAlgebra::vector_space< DifferentiationType<TScalar,TVector>, TScalar>,
LinearAlgebra::totally_ordered<DifferentiationType<TScalar, TVector> >,
LinearAlgebra::equality_comparable<DifferentiationType<TScalar, TVector> >,
LinearAlgebra::multiplicative<DifferentiationType<TScalar, TVector> >
{
    typedef TScalar ScalarType;
    typedef TVector VectorType;
    typedef DifferentiationType DT;
    
    ScalarType s;
    VectorType v;
    
    DifferentiationType(){};
    DifferentiationType(ScalarType t, VectorType w):s(t),v(w){};
    explicit DifferentiationType(ScalarType t):s(t){v.fill(0.);}
    explicit DifferentiationType(ScalarType t, size_t i):s(t){v.fill(0.); assert(i<v.size()); v[i]=1;}
    
    const DT & operator *= (ScalarType t) {s*=t; v*=t; return *this;}
    const DT &  operator /= (ScalarType t) {s/=t; v/=t; return *this;}
    const DT &  operator += (ScalarType t) {s+=t; return *this;}
    const DT &  operator -= (ScalarType t) {s-=t; return *this;}
    bool operator < (ScalarType t) {return s<t;}
    
    const DT &  operator *= (const DT & y) {v= v*y.s + y.v*s; s*=y.s; return *this;}
    const DT &  operator /= (const DT & y) {v= v/y.s - y.v*s/(y.s*y.s); s/=y.s; return *this;}
    const DT &  operator += (const DT & y) {s+=y.s; v+=y.v; return *this;}
    const DT &  operator -= (const DT & y) {s-=y.s; v-=y.v; return *this;}
    
    friend DT operator / (ScalarType s, const DT & y){return DT(s)/y;}
    DT operator -() const {return DT(-s,-v);}

    // Note: only base scalars are compared.
    bool operator < (const DT & y) const  {return s<y.s;}
    bool operator == (const DT & y) const {return s==y.s;}
    
    friend DT min(const DT & x, ScalarType s){return x.s < s ? x : DT(s);}
    friend DT sqrt(const DT & x){const ScalarType ss = sqrt(x.s); return ss>0 ? DT(ss,x.v/(2.*ss)) : DT(0.);};
    friend DT fabs(const DT & x){return x.s>=0 ? x : DT(-x.s,-x.v);}
    friend DT pow(const DT &x, ScalarType p){const ScalarType xp = pow(x.s,p-1); return DT(xp*x.s, p*xp*x.v);}
    //friend std::ostream & operator << (std::ostream & os, const ScalarType & y){return os << "{" << y.s << "," << y.v << "}";}
};



#endif
