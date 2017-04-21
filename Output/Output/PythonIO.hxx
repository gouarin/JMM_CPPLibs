// Copyright 2017 Jean-Marie Mirebeau, University Paris-Sud, CNRS, University Paris-Saclay
// Distributed WITHOUT ANY WARRANTY. Licensed under the Apache License, Version 2.0, see http://www.apache.org/licenses/LICENSE-2.0

#ifndef PythonIO_hxx
#define PythonIO_hxx

void PythonIO::PySetArray(KeyCRef key, ndarray arr) {
    if(arr.get_dtype() != dtype::get_builtin<double>())
        ExceptionMacro("PythonIO error : ndarray " << key << " is not made of doubles");
    
    SetDefined(key);
    RawElement & raw = rawElems.find(key)->second;
    raw.Clear(currentSetter);
    const int ndims = arr.get_nd();
    raw.dims.resize(ndims);
    DiscreteType stride = sizeof(ScalarType);
    
    for(int i=0; i<ndims; ++i){
        const int j=ndims-i-1;
        raw.dims[i] = arr.shape(j);
        
        //TODO : FIX ME.  Support arbitrary strides.
        if(arr.strides(j)!=stride){
            ExceptionMacro("PythonIO error : Sorry, only canonical strides are supported.");}
        stride*=raw.dims[i];
    }
    const DiscreteType flattenedLength = raw.FlattenedLength();
    raw.data.resize(flattenedLength);
    
    ScalarType * beg = reinterpret_cast<ScalarType*>(arr.get_data());
    std::copy(beg,beg+flattenedLength, raw.data.begin());
}

auto PythonIO::PyGetArray(KeyCRef key) const -> ndarray {
    const RawElement & raw = GetRaw(key);
    if(raw.IsString()) {ExceptionMacro("PythonIO error : field " << key << " is a string, not an array");}
    if(raw.dims.empty()){ExceptionMacro("PythonIO error : field " << key << " is a scalar, not an array");}
    boost::python::tuple shape;
    for(DiscreteType dim : raw.dims){
        shape = boost::python::extract<boost::python::tuple>(boost::python::make_tuple(dim)+shape);}
    ndarray result = boost::python::numpy::zeros(shape, dtype::get_builtin<double>());
    
    namespace p=boost::python;
    std::cout << p::extract<char const *>(p::str(shape)) << std::endl;
    std::cout << p::extract<char const *>(p::str(result)) << std::endl;

    assert(raw.FlattenedLength()==raw.data.size());
    ScalarType * beg = reinterpret_cast<ScalarType*>(result.get_data());
    std::copy(raw.data.begin(),raw.data.end(),beg);
    return result;
}

std::string PythonIO::GetComputedKeys() const {
    std::ostringstream oss;
    oss << "(";
    for(const auto & a : rawElems){
        if(a.second.setter == SetterTag::Compute) {
            oss << "(" << a.first << "," <<
            (a.second.IsString() ? "String" : a.second.IsScalar() ? "float" : "array")
            << "),";}
    }
    oss << ")";
    return oss.str();
}
#endif /* PythonIO_hxx */
