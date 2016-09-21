//
//  VirtualConstIterator.h
//  LiftedFastMarching
//
//  Created by Jean-Marie Mirebeau on 30/08/2016.
//
//

#ifndef VirtualConstIterator_h
#define VirtualConstIterator_h

template<typename TElement>
struct VirtualConstIterator {
    typedef TElement ElementType;
    virtual ~VirtualConstIterator(){};
    
    virtual void operator++() = 0;
    bool IsAtEnd() const {return status==AtEnd;}
    const ElementType & operator*() const {return current;}
    VirtualConstIterator(int _status, const ElementType & _current):status(_status),current(_current){}
protected:
    enum {AtEnd=0};
    int status;
    ElementType current;
};

#endif /* VirtualConstIterator_h */
