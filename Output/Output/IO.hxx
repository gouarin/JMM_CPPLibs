// Copyright 2017 Jean-Marie Mirebeau, University Paris-Sud, CNRS, University Paris-Saclay
// Distributed WITHOUT ANY WARRANTY. Licensed under the Apache License, Version 2.0, see http://www.apache.org/licenses/LICENSE-2.0

#ifndef IO_hxx
#define IO_hxx

// Input

template<typename Base>
std::string IO_<Base>::GetString(KeyCRef key, const std::string & val, int verb) const {
    if(Base::HasField(key)) return Base::GetString(key);
    if(this->verbosity>=verb) Msg() << "Field " << key << " defaults to " << val << "\n";
    this->visitedUnset.erase(key); this->defaulted.insert(key);
    return val;
}

template<typename Base> template<typename T>
T IO_<Base>::Get(KeyCRef key) const {
    const auto dimsPtr = this->template GetDimsPtr<T>(key);
    if( !dimsPtr.first.empty() )
        ExceptionMacro("IO::Get error : Field " << key << " has incorrect dimensions for this object.");
    return *dimsPtr.second;
}

template<typename Base> template<typename T>
T IO_<Base>::Get(KeyCRef key, const T & val, int verb) const {
    if(this->HasField(key)) return Get<T>(key);
    if(this->verbosity>=verb) Msg() << "Field " << key << " defaults to " << val << "\n";
    this->visitedUnset.erase(key); this->defaulted.insert(key);
    return val;
}

template<typename Base> template<typename T>
std::vector<T> IO_<Base>::GetVector(KeyCRef key) const {
    return GetArray<T, 1>(key);
}

template<typename Base> template<typename T, size_t d>
typename Base::template Array<T, d> IO_<Base>::GetArray(KeyCRef key) const {
    auto dimsPtr = this->template GetDimsPtr<T>(key);
    auto & dims = dimsPtr.first;
    
#if __IgnoreTrailingSingletonDimensions
    if(dims.size()<d) dims.resize(d,1);
#endif
    
    if(dims.size()!=d)
        ExceptionMacro("IO::GetArray error : Field " << key << " has depth " << dims.size() <<
                       " instead of expected " << d << ".");
    
    Array<T,d> result;
    std::copy(dims.begin(), dims.end(), result.dims.begin());
    result.resize(result.dims.ProductOfCoordinates());
    
    switch (this->arrayOrdering) {
        case ArrayOrdering::Default: {
            std::copy(dimsPtr.second, dimsPtr.second+result.size(), result.begin());
            break;}
        case ArrayOrdering::Transposed: {
            const TransposeVals<T,d> TrVals(result.dims, dimsPtr.second);
            result.dims = TransposeDims(result.dims);
            for(size_t i=0; i<result.size(); ++i) result[i] = TrVals((DiscreteType)i);
            break;}
        case ArrayOrdering::Reversed: {
            const ReverseVals<T,d> RevVals(result.dims, dimsPtr.second);
            result.dims = ReverseDims(result.dims);
            for(size_t i=0; i<result.size(); ++i) result[i] = RevVals((DiscreteType)i);
            break;}
    }
    return result;
}

template<typename Base> template<typename T>
auto IO_<Base>::GetDimensions(KeyCRef key) const -> std::vector<DiscreteType> {
    auto dims = this->template GetDimsPtr<T>(key).first;
    switch (this->arrayOrdering) {
        case ArrayOrdering::Default: return  dims;
        case ArrayOrdering::Transposed: return TransposeDims(dims);
        case ArrayOrdering::Reversed: return ReverseDims(dims);
    }
}

// Output

template<typename Base> template<typename T, size_t d>
void IO_<Base>::Set(KeyCRef key, DimType<d> dims, const T* pVals){
    struct {const T *pVal; const T & operator()(DiscreteType i) const {return pVal[i];}} fVal{pVals};
    Base::template Set<T,d>(key, dims, fVal);
    
}

template<typename Base> template<typename T>
void IO_<Base>::Set(KeyCRef key, const T & val) {
    this->template Set<T,0>(key,{},&val);
}

template<typename Base> template<typename T>
void IO_<Base>::SetVector(KeyCRef key, const std::vector<T> & val) {
    Set<T,1>(key,{(DiscreteType)val.size()},&val[0]);
}

template<typename Base> template<typename T, size_t d>
void IO_<Base>::SetArray(KeyCRef key, const Array<T, d> & val) {
    switch (this->arrayOrdering) {
        case ArrayOrdering::Default: return Set<T,d>(key,val.dims,&val[0]);
        case ArrayOrdering::Transposed:
            return Base::template Set<T,d>(key, TransposeDims(val.dims),
                                           TransposeVals<T,d>(val.dims,&val[0]) );
        case ArrayOrdering::Reversed:
            return Base::template Set<T,d>(key, ReverseDims(val.dims),
                                           ReverseVals<T,d>(val.dims,&val[0]) );
    }
}

// Transposition

template<typename Base> template<typename V>
V IO_<Base>::TransposeDims(const V & dims) {
    V result=dims;
    if(dims.size()>=2)
        std::swap(result[0], result[1]);
    return result;
}

template<typename Base> template<typename T, size_t d>
struct IO_<Base>::TransposeVals {
    const T * const pVals;
    Array<T,d> a, ta;
    TransposeVals(const DimType<d> & dims, const T * pVals_)
    :pVals(pVals_){a.dims=dims; ta.dims=IO_::TransposeDims(dims);}
    const T & operator()(DiscreteType index) const {
        return pVals[a.Convert(TransposeDims(ta.Convert(index)))];}
};

// Reversal

template<typename Base> template<typename V>
V IO_<Base>::ReverseDims(const V & dims) {
    V result=dims;
    const int n=(int)dims.size();
    for(int i=0; i<n/2; ++i)
        std::swap(result[i],result[n-1-i]);
    return result;
}

template<typename Base> template<typename T, size_t d>
struct IO_<Base>::ReverseVals {
    const T * const pVals;
    Array<T,d> a, ta;
    ReverseVals(const DimType<d> & dims, const T * pVals_)
    :pVals(pVals_){a.dims=dims; ta.dims=IO_::ReverseDims(dims);}
    const T & operator()(DiscreteType index) const {
        return pVals[a.Convert(ReverseDims(ta.Convert(index)))];}
};

// Destruction

template<typename Base>
IO_<Base>::~IO_(){
    this->UsageReport();
}
#endif /* IO_hxx */
