//
//  SquareCube.h
//  AmongConvex2
//
//  Created by Jean-Marie Mirebeau on 25/07/13.
//  Copyright (c) 2013 Jean-Marie Mirebeau. All rights reserved.
//

#ifndef AmongConvex2_SquareCube_h
#define AmongConvex2_SquareCube_h


template<int i, typename T>
inline T pow(T t){
    if(i<0) return T(1)/pow<-i,T>(t);
    if(i==0) return T(1);
    if(i==1) return t;
    const int j=i%2;
    const T s = pow< (i>>1),T>(t);
    if(j==0) return s*s;
    else return t*s*s;
}


template<typename T>
inline T square(T t){return t*=t;}

template<typename T>
inline T cube(T t){return t*=square(t);}

#endif
