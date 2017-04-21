// Copyright 2017 Jean-Marie Mirebeau, University Paris-Sud, CNRS, University Paris-Saclay
// Distributed WITHOUT ANY WARRANTY. Licensed under the Apache License, Version 2.0, see http://www.apache.org/licenses/LICENSE-2.0

#ifndef MathematicaIO_h
#define MathematicaIO_h

#include "BaseIO.h"

/*
A simple interface for making Python data available to the c++ code.
*/

struct PythonIO : BaseIO {
    // An adaptation of the file IO:
    PythonIO(){};
    typedef boost::python::numpy::ndarray ndarray;
    ndarray PyGetArray(KeyCRef) const;
    void PySetArray(KeyCRef, ndarray);
    std::string GetComputedKeys() const;
protected:
    typedef boost::python::numpy::dtype dtype;
};


#include "PythonIO.hxx"

#endif /* MathematicaIO_h */
