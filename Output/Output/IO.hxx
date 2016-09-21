//
//  IO.hxx
//  ExternLFM
//
//  Created by Jean-Marie Mirebeau on 21/09/2016.
//
//

#ifndef IO_hxx
#define IO_hxx

// Input

template<typename Base>
std::string IO_<Base>::GetString(KeyCRef key, const std::string & val) const {
    if(Base::HasField(key)) return Base::GetString(key);
    if(verbosity) Msg() << "Field " << key << " defaults to " << val << "\n";
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
T IO_<Base>::Get(KeyCRef key, const T & val) const {
    if(this->HasField(key)) return Get<T>(key);
    if(verbosity) Msg() << "Field " << key << " defaults to " << val << "\n";
    return val;
}

template<typename Base> template<typename T, size_t d>
TraitsIO::Array<T, d> IO_<Base>::GetArray(KeyCRef key) const {
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
    
    if(transposeFirstTwoCoordinates){
        const TransposeVals<T,d> TrVals(result.dims, dimsPtr.second);
        result.dims = TransposeDims(result.dims);
        for(size_t i=0; i<result.size(); ++i) result[i] = TrVals(i);
    } else std::copy(dimsPtr.second, dimsPtr.second+result.size(), result.begin());
    
    return result;
}

template<typename Base> template<typename T>
std::vector<TraitsIO::DiscreteType> IO_<Base>::GetDimensions(KeyCRef key) const {
    auto dims = this->template GetDimsPtr<T>(key).first;
    if(transposeFirstTwoCoordinates) return TransposeDims(dims);
    else return dims;
}

// Output

template<typename Base> template<typename T, size_t d>
void IO_<Base>::Set(KeyCRef key, DimType<d> dims, const T* pVals){
    struct {const T *pVal; const T & operator()(DiscreteType i) const {return pVal[i];}} fVal{pVals};
    Base::template Set<T,d>(key, dims, fVal);
    
}

template<typename Base> template<typename T>
void IO_<Base>::Set(KeyCRef key, const T & val) {
    this-> template Set<T,0>(key,{},&val);
}

template<typename Base> template<typename T>
void IO_<Base>::SetVector(KeyCRef key, const std::vector<T> & val) {
    Set<T,1>(key,{(DiscreteType)val.size()},&val[0]);
}

template<typename Base> template<typename T, size_t d>
void IO_<Base>::SetArray(KeyCRef key, const Array<T, d> & val) {
    if(transposeFirstTwoCoordinates)
        Base::template Set<T,d>(key, TransposeDims(val.dims), TransposeVals<T,d>(val.dims,&val[0]) );
    else Set<T,d>(key,val.dims,&val[0]);
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

// Destruction

template<typename Base>
IO_<Base>::~IO_(){
    if(!this->unused.empty()){
        std::ostringstream oss;
        for(KeyCRef key : this->unused) oss << key << " ";
        this->SetString("unused", oss.str());
        Msg() << "Unused fields : " << oss.str() << "\n";
    }
    if(!this->defaulted.empty()){
        std::ostringstream oss;
        for(KeyCRef key : this->defaulted) oss << key << " ";
        this->SetString("defaulted", oss.str());
    }
}
#endif /* IO_hxx */
