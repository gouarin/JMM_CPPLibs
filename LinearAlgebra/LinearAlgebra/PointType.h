
//
//  PointType.h
//  AmongConvex2
//
//  Created by Jean-Marie Mirebeau on 23/07/13.
//  Copyright (c) 2013 Jean-Marie Mirebeau. All rights reserved.
//

#ifndef AmongConvex2_PointType_h
#define AmongConvex2_PointType_h

#include "PointBaseType.h"
#include <array>
#include <numeric>
#include <math.h>

namespace LinearAlgebra {
    
template<typename TComponent, size_t VDimension>
struct Point :
PointBase<TComponent,VDimension>
{
    typedef TComponent ComponentType;
    static const int Dimension = VDimension;
    
    typedef PointBase<TComponent,VDimension> PointBaseType;
    
    Point(std::initializer_list<TComponent> c):PointBaseType(c){};
    Point(){};
    
    static Point FromOrigin(const PointBaseType & p){return Point(p);}
    static Point Constant(ComponentType c){return PointBaseType::Constant(c);}
    template<typename ComponentType2>
    static Point CastCoordinates(const Point<ComponentType2,Dimension> & p){
        return Point(PointBaseType::CastCoordinates(p));}
    
    template<size_t n>
    static Point
    Barycenter(const std::array<Point,n> & points, const std::array<TComponent,n> & weights) {
        assert(fabs(std::accumulate(weights.begin(), weights.end(), ComponentType(0)) - 1) < 1e-6);
        Point p;
        for(int i=0; i<Dimension; ++i){
            p[i]=0;
            for(int j=0; j<n; ++j)
                p[i]+=weights[j]*points[j][i];
        }
        return p;
    }
    
    static Point Barycenter(const Point & p, const Point & q, ComponentType t) {
        assert(0<=t && t<=1);
        Point b;
        for(int i=0; i<Dimension; ++i)
            b[i]=(1-t)*p[i]+t*q[i];
        return b;
    }
protected:
    Point(const PointBaseType & p): PointBaseType(p){};
};
    
}
#endif
