//
//  PrintContainers.h
//  AmongConvex2
//
//  Created by Jean-Marie Mirebeau on 25/07/13.
//  Copyright (c) 2013 Jean-Marie Mirebeau. All rights reserved.
//

#ifndef AmongConvex2_PrintContainers_h
#define AmongConvex2_PrintContainers_h

#include <iostream>
#include <map>


template<typename A,typename B>
std::ostream & operator << (std::ostream & f, std::pair<A,B> p){
    return f << "{" << p.first << "," << p.second << "}";
}

template <typename ForwardIterator>
void print_range(std::ostream & f, ForwardIterator first, ForwardIterator last) {
    f<<"{";
    if(first!=last) f<<*first++;
    while(first!=last) f << "," << *first++;
    f<<"}";
}

template <typename ContainerType>
void print_container(std::ostream & f, const ContainerType & u) {
    print_range(f,begin(u),end(u));
}

/*
template<typename KeyType, typename ValueType>
std::ostream & operator << (std::ostream & f, std::multimap<KeyType,ValueType> & u){
    if(u.empty()) {f << "{}"; return f;}
    bool first=true;
    KeyType oldKey = u.begin()->first;
    for(const auto & keyVal : u){
        if(first){
            first=false;
            f << "{{" << keyVal.first << ",{" << keyVal.second;
        } else if(keyVal.first != oldKey){
            oldKey = keyVal.first;
            f << "}},{" << keyVal.first << "'{" << keyVal.second;
        } else {
            f << "," << keyVal.second;
        }
    }
    f << "}}";
}
 */
#endif
