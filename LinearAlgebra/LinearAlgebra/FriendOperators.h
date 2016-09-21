//
//  FriendOperators.h
//  AmongConvex2
//
//  Created by Jean-Marie Mirebeau on 24/07/13.
//  Copyright (c) 2013 Jean-Marie Mirebeau. All rights reserved.
//

#ifndef AmongConvex2_FriendOperators_h
#define AmongConvex2_FriendOperators_h

// Operators for a vector space V on a field K. (or a module on a ring, but be careful with divisions.)
// K elements are passed by value, since we assume that it is a machine type (e.g. double, int)

namespace LinearAlgebra {
    
template<class K> class multiplicative {
    friend K operator * (K k, const K & l){k*=l; return k;}
    friend K operator / (K k, const K & l){k/=l; return k;}
};

// Define u+=v, u-=v, u*=k, u/=k
template <class V, class K> class vector_space {
    friend V operator + (V u, const V &v) {u+=v; return u;}
    friend V operator - (V u, const V &v) {u-=v; return u;}
    friend V operator * (V u, K k)  {u*=k; return u;}
    friend V operator * (K k, V u)  {u*=k; return u;}
    friend V operator / (V u, K k)  {u/=k; return u;}
};

// define p+=v, p-=v
template <class P, class V> class offsettable {
    friend P operator +(P p, const V & v){p+=v; return p;}
    friend P operator +(const V & v, P p){p+=v; return p;}
    friend P operator -(P p, const V & v){p-=v; return p;}
};

template <class P, class V, class K> class affine_space :
offsettable< P, V>,
vector_space<V, K>
{};

// define a<b, a==b
template <class A> class totally_ordered {
    friend bool operator > (const A &a, const A &b){return (b<a);};
    friend bool operator <= (const A &a, const A &b){return !(b<a);}
    friend bool operator >= (const A &a, const A &b){return !(a<b);}
};

// define a==b
template <class A> class equality_comparable {
    friend bool operator != (const A &a, const A &b){return !(a==b);}
};
    
}

#endif
