// Copyright 2017 Jean-Marie Mirebeau, University Paris-Sud, CNRS, University Paris-Saclay
// Distributed WITHOUT ANY WARRANTY. Licensed under the Apache License, Version 2.0, see http://www.apache.org/licenses/LICENSE-2.0


#ifndef BaseIO_hxx
#define BaseIO_hxx

// ---- Data container ----

struct BaseIO::RawElement {
    // Case of a string field
    std::string str;
    
    // case of a numerical field
    std::vector<DiscreteType> dims;
    std::vector<ScalarType> data;
    
    SetterTag setter = SetterTag::Unknown;
    
    bool IsString() const {return !str.empty();}
    bool IsScalar() const {return !IsString() && dims.size()==0;}
    void Clear(SetterTag tag) {str.clear(); dims.clear(); data.clear();setter=tag;}
    RawElement(SetterTag tag):setter(tag){};
    DiscreteType FlattenedLength() const {
        if(IsString()) {ExceptionMacro("BaseIO::RawElement::FlattenedLength error, field is a string");}
        DiscreteType res=1;
        for(auto dim : dims) res*=dim;
        return res;
    }
};


// ----- Public fields -----

bool BaseIO::HasField(KeyCRef key) const {
    if (rawElems.find(key) != rawElems.end()) return true;
    visitedUnset.insert(key);
    return false;
}

bool BaseIO::EraseField(KeyCRef key) {
    unused.erase(key);
    return rawElems.erase(key);
}

std::string BaseIO::GetString(KeyCRef key) const {
    const auto & val = GetRaw(key);
    if (!val.data.empty())
        ExceptionMacro("BaseIO import error : field " << key << " is not a string.");
    unused.erase(key);
    return val.str;
}

void BaseIO::SetString(KeyCRef key, const std::string & val) {
    SetDefined(key);
    RawElement & raw = rawElems.find(key)->second;
    raw.str = val;
}

// ----- Protected fields -----

void BaseIO::UsageReport(){
    if(!unused.empty()){
        std::ostringstream ossUser, ossCompute;
        for(KeyCRef key : unused) {
            switch(GetSetter(key)){
                case SetterTag::User: ossUser << key << " "; break;
                case SetterTag::Compute: ossCompute << key << " "; break;
                case SetterTag::Unknown: ExceptionMacro("BaseIO Error : unknown setter for key " << key);
            }
        }
        const std::string & fromUser = ossUser.str(), fromCompute = ossCompute.str();
        SetString("unusedFromUser", fromUser);
        SetString("unusedFromCompute", fromCompute);
        
        if(verbosity>=1 && !fromUser.empty())    _Msg<true,BaseIO>(this) << "Unused fields from user: " << fromUser << "\n";
        if(verbosity>=2 && !fromCompute.empty()) _Msg<false,BaseIO>(this) << "Unused fields from compute : " << fromCompute << "\n";
    }
    if(!defaulted.empty()){
        std::ostringstream oss;
        for(KeyCRef key : defaulted) oss << key << " ";
        SetString("defaulted", oss.str());
        if(verbosity>=2) _Msg<false,BaseIO>(this) << "Defaulted fields : " << oss.str() << "\n";
    }
    if(!visitedUnset.empty()){
        std::ostringstream oss;
        for(KeyCRef key : visitedUnset) oss << key << " ";
        SetString("visitedUnset", oss.str());
        if(verbosity>=3) _Msg<false,BaseIO>(this) << "Visited but unset fields : " << oss.str() << "\n";
    }
}

auto BaseIO::GetSetter(KeyCRef key) const -> SetterTag {
    return GetRaw(key,false).setter;
}

const BaseIO::RawElement & BaseIO::GetRaw(KeyCRef key, bool forUse) const {
    const auto it = rawElems.find(key);
    if (it == rawElems.end()) ExceptionMacro("BaseIO import error : field " << key << " not found.");
    if(forUse) unused.erase(key);
    return it->second;
}


void BaseIO::SetDefined(KeyCRef key) {
    auto it = rawElems.find(key);
    if (it != rawElems.end() && verbosity>=1){
        _Msg<true,BaseIO>(this) << "BaseIO: redefining field " << key << ".\n";}
    rawElems.insert(std::pair<std::string,RawElement>(key,RawElement(currentSetter)));
    unused.insert(key);
}

template<typename T> std::pair<std::vector<BaseIO::DiscreteType>, const T*> BaseIO::GetDimsPtr(KeyCRef key) const {
    const auto & raw = GetRaw(key);
    if(raw.IsString()) {ExceptionMacro("IO error : key " << key << " is a string, not a numerical field.");}
    auto dims = raw.dims;
    static_assert(sizeof(T) % sizeof(ScalarType) == 0, "Field is not made of scalars.");
    const DiscreteType sizeRatio = sizeof(T) / sizeof(ScalarType);
    if (!std::is_same<T, ScalarType>::value) {
        if (dims.empty() || dims[0] != sizeRatio)
            ExceptionMacro("PythonIO input error first dimension " << (dims.empty() ? 0 : dims[0]) << " of field "
                           << key << " does not match expected value " << sizeRatio << ".");
        dims.erase(dims.begin());
    }
    return {dims, reinterpret_cast<const T*>(&raw.data[0]) };
}

template<typename T, size_t d, typename F> void BaseIO::Set(KeyCRef key, DimType<d> dims, const F & vals) {
    SetDefined(key);
    
    static_assert(sizeof(T) % sizeof(ScalarType) == 0, "Type is not built of scalars.");
    const ScalarType sizeRatio = sizeof(T) / sizeof(ScalarType);
    const DiscreteType size = dims.ProductOfCoordinates();
    
    RawElement & raw = rawElems.find(key)->second;
    raw.data.resize(sizeRatio*size);
    
    T* input = reinterpret_cast<T*>(&raw.data[0]);
    for (DiscreteType i = 0; i<size; ++i)
        input[i] = vals(i);
    
    if(!std::is_same<T, ScalarType>::value) raw.dims.push_back(sizeRatio);
    for(auto dim : dims) raw.dims.push_back(dim);
}




#endif /* BaseIO_h */
