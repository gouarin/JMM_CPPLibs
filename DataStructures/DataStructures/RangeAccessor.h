//
//  RangeAccessor.h
//  LiftedFastMarching
//
//  Created by Jean-Marie Mirebeau on 30/08/2016.
//
//

#ifndef RangeAccessor_h
#define RangeAccessor_h

template<typename Pointer> class RangeAccessor {
    Pointer _begin, _end;
public:
    typedef Pointer pointer;
    RangeAccessor(pointer __begin, pointer __end):_begin(__begin),_end(__end){};
    
    pointer begin() const {return _begin;}
    pointer end() const {return _end;}
    size_t  size() const {return _end - _begin;}
    
    typedef decltype(*_begin) value_type;
    value_type & operator[](size_t i) {assert(i<size()); return _begin[i];}
    const value_type & operator[](size_t i) const {assert(i<size()); return _begin[i];}
    value_type & back() {assert(size()>0); pointer ptr(_end); return *(--ptr);}
    const value_type & back() const {assert(size()>0); pointer ptr(_end); return *(--ptr);}
    value_type & front() {assert(size()>0); return *_begin;}
    const value_type & front() const {assert(size()>0); return *_begin;}
};

#endif /* RangeAccessor_h */
