//
//  MexIO.hxx
//  LiftedFastMarching_MatlabBuildTest
//
//  Created by Jean-Marie Mirebeau on 17/09/2016.
//
//

#ifndef MexIO_hxx
#define MexIO_hxx

// ----- Public -----

template<bool warn> BaseIO::Msg_<warn>::~Msg_(){
    std::string str = oss.str();
    if(warn) mexWarnMsgTxt(str.c_str());
        else {
            size_t start_pos = 0; const std::string from="%", to ="%%";
            while((start_pos = str.find(from, start_pos)) != std::string::npos) {
                str.replace(start_pos, from.length(), to);
                start_pos += to.length(); // Handles case where 'to' is a substring of 'from'
                
            }
            printf(str.c_str());
        }
    mexEvalString("drawnow");
    mexEvalString("pause(0.0001)");
};


bool BaseIO::HasField(KeyCRef key) const {
    if( mxGetFieldNumber(mxInput, key.c_str()) >= 0 ) return true;
    defaulted.insert(key);
    return false;
}

std::string BaseIO::GetString(KeyCRef key) const {
    if(!HasField(key)) ExceptionMacro("GetString error : field" << key << " does not exist.");
    unused.erase(key);
    
    const mxArray * field = mxGetField(mxInput,0,key.c_str());
    if(field==NULL) ExceptionMacro("GetString error : field " << key << " could not be imported.");
    
    const char * ans = mxArrayToString(field);
    if(ans==NULL) ExceptionMacro("GetString error : field " << key << " is not a string.");
    
    return ans;
}

void BaseIO::SetString(KeyCRef key, std::string val){
    SetField(key, mxCreateString(val.c_str()));}

// ----- Protected -----

template<typename T> std::pair<std::vector<TraitsIO::DiscreteType>,const T*> BaseIO::GetDimsPtr(KeyCRef key) const {
    static_assert(sizeof(T)%sizeof(ScalarType)==0, "Error : element is not made of scalars.");
    if(!HasField(key)) ExceptionMacro("IO error : field " << key << " does not exist.");
    unused.erase(key);
    
    mxArray * p = mxGetField(mxInput,0,key.c_str());
    std::string errMsg = "GetScalarArrayRaw(" + key + ") error : ";
    
    if(!mxIsNumeric(p)) ExceptionMacro(errMsg << "Not a numeric array.");
    if(mxIsComplex(p)) ExceptionMacro(errMsg << "Complex arrays not supported.");
    if(mxGetClassID(p) != mxDOUBLE_CLASS) ExceptionMacro(errMsg << "Array elements must be of double type.");
    
    const mwSize * mxSize = mxGetDimensions(p);
    const mwSize ArrayDimension = mxGetNumberOfDimensions(p);
    
    std::vector<DiscreteType> dims(mxSize, mxSize+ArrayDimension);
    if(!std::is_same<T, ScalarType>::value){
        const DiscreteType sizeRatio = sizeof(T)/sizeof(ScalarType);
        if(dims.empty() || dims[0]!=sizeRatio)
            ExceptionMacro("IO input error : first dimension " << (dims.empty() ? 0 : dims[0]) << " of field " << key << " does not match expected value " << sizeRatio << ".");
        dims.erase(dims.begin());
    }

    while(!dims.empty() && dims.back()==1)
        dims.pop_back(); // Eliminate trailing singleton dimensions.
    
    return {dims,reinterpret_cast<const T*>(mxGetPr(p))};
}

template<typename T, size_t d, typename F> void BaseIO::Set(KeyCRef key, DimType<d> dims, const F & f){
    static_assert(sizeof(T)%sizeof(ScalarType)==0, "IO error : field must be made of scalars.");
    const int sizeRatio = sizeof(T)/sizeof(ScalarType);
    std::vector<mwSize> mxDims;
    if(!std::is_same<T, ScalarType>::value) mxDims.push_back(sizeRatio);
    for(DiscreteType dim : dims)            mxDims.push_back(dim);
    for(size_t i=mxDims.size(); i<2; ++i)   mxDims.push_back(1);
    
    mxArray * pArr = mxCreateNumericArray(mxDims.size(), &mxDims[0], mxDOUBLE_CLASS, mxREAL);
    
    if(pArr==NULL)
        ExceptionMacro("IO error : could not create mex array of dimensions" << dims
                       << " for elements of size " << sizeof(T)/sizeof(ScalarType) << ".");
    
    T * pVal = reinterpret_cast<T*>(mxGetPr(pArr));
    const DiscreteType size = dims.ProductOfCoordinates();
    for(DiscreteType i=0; i<size; ++i)
        pVal[i] = f(i);
    
    SetField(key, pArr);
}


void BaseIO::SetField(KeyCRef key, mxArray * pVal){
    if(exported.find(key)!=exported.end())
        ExceptionMacro("FileIO error : exporting two fields with identical names");
    exported.insert(key);

    if(pVal==NULL)  ExceptionMacro("IO error : could not create Matlab(R) representation of field" << key << ".");
    
    if(key.empty())     ExceptionMacro("IO output error : trying to export field with empty name.");
    
    const int fieldNum = mxAddField(*mxOutput, key.c_str());
    if(fieldNum==-1)    ExceptionMacro("IO error, could not add output field " << fieldNum << ".");
    
    mxSetFieldByNumber(*mxOutput, 0, fieldNum, pVal);
}

BaseIO::BaseIO(const mxArray * mxInput_, mxArray ** mxOutput_)
:mxInput(mxInput_), mxOutput(mxOutput_){
    
    // Examine input
    if(!mxIsStruct(mxInput))
        ExceptionMacro("Input not a struct");
    
    for(size_t i = 0; i<mxGetNumberOfDimensions(mxInput); ++i)
        if(mxGetDimensions(mxInput)[i] != 1)
            ExceptionMacro("Input structure is not flat");
    
    for(int i=0; i<mxGetNumberOfFields(mxInput); ++i)
        unused.insert(std::string(mxGetFieldNameByNumber(mxInput, i)));
    
    
    // Prepare output
    const mwSize nDims = 2;
    const mwSize dims[2]={1,1};
    *mxOutput = mxCreateStructArray(nDims,dims,0,nullptr);
}


#endif /* MexIO_hxx */
