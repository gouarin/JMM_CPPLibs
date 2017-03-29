// Copyright 2017 Jean-Marie Mirebeau, University Paris-Sud, CNRS, University Paris-Saclay
// Distributed WITHOUT ANY WARRANTY. Licensed under the Apache License, Version 2.0, see http://www.apache.org/licenses/LICENSE-2.0

#ifndef FileIO_hxx
#define FileIO_hxx

// Protected fields

struct BaseIO::InputFormatElement {
    // Case of a string field
    std::string str;
    
    // case of a numerical field
    std::vector<DiscreteType> dims;
    size_t pos = std::numeric_limits<size_t>::max();
};

const BaseIO::InputFormatElement & BaseIO::GetInputFormat(KeyCRef key) const {
    const auto it = inputFormat.find(key);
    if(it==inputFormat.end()) ExceptionMacro("FileIO import error : field " + key + " not found.");
    unused.erase(key);
    return it->second;
}


void BaseIO::SetExported(KeyCRef key){
    if(exported.find(key)!=exported.end())
        ExceptionMacro("FileIO error : exporting two fields with identical names");
    exported.insert(key);
}

template<typename T> std::pair<std::vector<BaseIO::DiscreteType>,const T*> BaseIO::GetDimsPtr(KeyCRef key) const {
    const auto & format = GetInputFormat(key);
    auto dims = format.dims;
    static_assert(sizeof(T)%sizeof(ScalarType)==0, "Field is not made of scalars.");
    const DiscreteType sizeRatio = sizeof(T)/sizeof(ScalarType);
    if(!std::is_same<T, ScalarType>::value){
        if(dims.empty() || dims[0]!=sizeRatio)
            ExceptionMacro("FileIO input error first dimension " << (dims.empty() ? 0 : dims[0]) << " of field " << key << " does not match expected value " << sizeRatio << ".");
        dims.erase(dims.begin());
    }
    return {dims, reinterpret_cast<const T*>(&inputData[format.pos])};
}

template<typename T, size_t d, typename F> void BaseIO::Set(KeyCRef key, DimType<d> dims, const F & vals) {
    SetExported(key);
    static_assert(sizeof(T)%sizeof(ScalarType)==0, "Type is not built of scalars.");
    const ScalarType sizeRatio = sizeof(T)/sizeof(ScalarType);
    
    const size_t pos = outputData.size();
    const DiscreteType size = dims.ProductOfCoordinates();
    outputData.resize(pos+sizeRatio*size);
    
    T* output = reinterpret_cast<T*>(&outputData[pos]);
    for(DiscreteType i=0; i<size; ++i)
        output[i] = vals(i);
    
    outputFormat << key << "\n";
    if(std::is_same<ScalarType, T>::value) outputFormat << d << "\n";
    else outputFormat << (d+1) << "\n" << sizeRatio << "\n";
    for(int i=0; i<d; ++i) outputFormat << dims[i] << "\n";
    outputFormat << "\n";
}


// Public fields

bool BaseIO::HasField(KeyCRef key) const {
    if( inputFormat.find(key)!=inputFormat.end() ) return true;
    defaulted.insert(key);
    return false;
}

std::string BaseIO::GetString(KeyCRef key) const {
    const auto & val = GetInputFormat(key);
    if(val.dims.size()>0 || val.pos!=std::numeric_limits<size_t>::max())
        ExceptionMacro("FileIO import error : field " << key << " is not a string.");
    return val.str;
}

void BaseIO::SetString(KeyCRef key, std::string val) {
    SetExported(key);
    outputFormat << key << "\n-1\n" << val << "\n\n";
}

// Construction and destruction


BaseIO::BaseIO(std::string inputPrefix, std::string outputPrefix_):outputPrefix(outputPrefix_){
    
    // ------- Read format file -------
    
    struct {
        void operator()(std::istream & is){
            is.ignore(std::numeric_limits<std::streamsize>::max(), is.widen('\n'));}
    } nextLine;
    
    size_t pos=0;
    std::ifstream ifs;
    ifs.open(inputPrefix+"_Format.txt");
    if(!ifs.is_open())
        ExceptionMacro("FileIO error : Input format file could not be opened "<<inputPrefix<<"_Format.txt.");
    
    while(!ifs.eof()){
        std::string key; ifs >> key; nextLine(ifs);
        if(key.empty()){
            if(ifs.eof()) break; // Allow for a final blank line
            else ExceptionMacro("FileIO Error empty field name");
        }
        if(ifs.eof() || ifs.fail()) ExceptionMacro("FileIO error : invalid input format file.");
        if(!unused.insert(key).second)
            ExceptionMacro("FileIO error : duplicate key " << key << " in format file.");
        
        int dim; ifs >> dim; nextLine(ifs);
        if(ifs.eof() || ifs.fail() || dim<-1)
            ExceptionMacro("FileIO error : invalid dimension value (" << std::to_string(dim) << ") for key " << key << "in input format file.");
        
        InputFormatElement format;
        if(dim==-1){
            ifs >> format.str; nextLine(ifs);
            if(ifs.eof() || ifs.fail())
                ExceptionMacro("FileIO error : invalid string value for key " << key << "in input format file.");
        } else {
            format.pos = pos;
            format.dims.resize(dim);
            DiscreteType size=1;
            for(int i=0; i<dim; ++i){
                DiscreteType & dimi = format.dims[i];
                ifs >> dimi; nextLine(ifs);
                if(ifs.eof() || ifs.fail() || dimi<0)
                    ExceptionMacro("FileIO error : invalid dimension for key " << key << " in input format file.");
                size*=dimi;
            }
            pos+=size;
        }
        inputFormat[key] = format;
        nextLine(ifs);
    }
    ifs.close();
    
    // ------- read data file --------
    if(pos==0) return;
    inputData.resize(pos);
    
    ifs.open(inputPrefix+"_Data.dat");
    if(!ifs.is_open())
        ExceptionMacro("FileIO error : input data file could not be opened " << inputPrefix << "_Data.dat");
    
    ifs.seekg(0, ifs.end);
    if(ifs.tellg()!=inputData.size()*sizeof(ScalarType))
        ExceptionMacro("FileIO error : input data file has incorrect size.");
    
    ifs.seekg(0, ifs.beg);
    ifs.read((char*)&inputData[0], inputData.size()*sizeof(ScalarType));
    ifs.close();
}

BaseIO::~BaseIO(){
    std::ofstream ofs;
    ofs.open(outputPrefix+"_Format.txt");
    ofs << outputFormat.str();
    ofs.close();
    
    ofs.open(outputPrefix+"_Data.dat",  std::ios::out | std::ios::binary);
    ofs.write((char*) &outputData[0], outputData.size()*sizeof(ScalarType));
}

#endif /* FileIO_hxx */
