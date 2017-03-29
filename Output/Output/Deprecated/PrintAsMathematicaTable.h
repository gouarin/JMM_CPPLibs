//
//  PrintAsMathematicaTable.h
//  OptimalTransport
//
//  Created by Jean-Marie Mirebeau on 05/02/2014.
//  Copyright (c) 2014 Jean-Marie Mirebeau. All rights reserved.
//

#ifndef OptimalTransport_PrintAsMathematicaTable_h
#define OptimalTransport_PrintAsMathematicaTable_h

#include <sstream>
#include "PrintContainers.h"
#include <cassert>

template <class ForwardIterator, class DimensionsType>
void PrintAsMathematicaTable(std::ostream & f,
                             ForwardIterator first,
                             ForwardIterator last,
                             DimensionsType N) {
    const int Dimension = N.size();
    int prod[Dimension];
    prod[0]=N[0];
    for(int k=1; k<Dimension; ++k)
        prod[k]=N[k]*prod[k-1];
    
    for(int k=0; k<Dimension; ++k) f<<"{";
    if(first!=last) f<<*first++;
    int i=1;
    for(; first!=last; ++first, ++i){
        for(int k=0; k<Dimension; ++k) if(i%(prod[k])==0) f<<"}";
        f<<",";
        for(int k=0; k<Dimension; ++k) if(i%(prod[k])==0) f<<"{";
        f << *first;
    }
    for(int k=0; k<Dimension; ++k) f<<"}";
    assert(i==prod[Dimension-1]); //size and format must match
}

// Definition of n intervals needs (n+1) integers (interval bounds).
template<typename IndicesContainer, typename DataContainer>
void PrintArrayOfArrays(std::ostream &os, const IndicesContainer & ind, const DataContainer & data){
    os << "{";
    assert(ind.size()>0);
    auto it = ind.begin();
    while(true){
        assert(0<=*it && *it<=data.size());
        const auto beginIt = &data[*it];
        const bool atBegin = it==ind.begin();
        if(++it==ind.end()) break;
        if(!atBegin) os << ",";
        assert(0<=*it && *it<=data.size());
        const auto endIt = &data[*it];
        print_range(os, beginIt, endIt);
    }
    os << "}";
}

#endif
