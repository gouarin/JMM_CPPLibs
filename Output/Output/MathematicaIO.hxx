// Copyright 2017 Jean-Marie Mirebeau, University Paris-Sud, CNRS, University Paris-Saclay
// Copyright 2017 Erik Bekkers, Eindhoven University of Technology, the Netherlands
// Distributed WITHOUT ANY WARRANTY. Licensed under the Apache License, Version 2.0, see http://www.apache.org/licenses/LICENSE-2.0

#ifndef MathematicaIO_hxx
#define MathematicaIO_hxx

// Protected fields

struct BaseIO::InputFormatElement {
	// Case of a string field
	std::string str;

	// case of a numerical field
	std::vector<DiscreteType> dims;
	size_t pos = std::numeric_limits<size_t>::max();

	// in case of vector valued vectors:
	DiscreteType sizeRatio;
};

const BaseIO::InputFormatElement & BaseIO::GetInputFormat(KeyCRef key) const {
	const auto it = inputOutputFormat.find(key);
	if (it == inputOutputFormat.end()) ExceptionMacro("MathematicaIO import error : field " + key + " not found.");
	unused.erase(key);
	return it->second;
}


void BaseIO::SetDefined(KeyCRef key) {
    static bool hasWarned=false;
    
    if (defined.find(key) != defined.end()){
        Msg_<true>() << "MathematicaIO warning: redefining field " << key <<
        (!hasWarned ? ", but memory is not freed until library is unloaded.\n" : "\n");
        hasWarned=true;
    }
	defined.insert(key);
    unused.insert(key);
}


template<typename T> std::pair<std::vector<BaseIO::DiscreteType>, const T*> BaseIO::GetDimsPtr(KeyCRef key) const {
	const auto & format = GetInputFormat(key);
	auto dims = format.dims;
	static_assert(sizeof(T) % sizeof(ScalarType) == 0, "Field is not made of scalars.");
	const DiscreteType sizeRatio = sizeof(T) / sizeof(ScalarType);
	if (!std::is_same<T, ScalarType>::value) {
		if (dims.empty() || dims[0] != sizeRatio)
			ExceptionMacro("MathematicaIO input error first dimension " << (dims.empty() ? 0 : dims[0]) << " of field " << key << " does not match expected value " << sizeRatio << ".");
		dims.erase(dims.begin());
	}
	return { dims, reinterpret_cast<const T*>(&inputOutputData[format.pos]) };
}

template<typename T, size_t d, typename F> void BaseIO::Set(KeyCRef key, DimType<d> dims, const F & vals) {
	SetDefined(key);
	
	static_assert(sizeof(T) % sizeof(ScalarType) == 0, "Type is not built of scalars.");
	const ScalarType sizeRatio = sizeof(T) / sizeof(ScalarType);
	const DiscreteType size = dims.ProductOfCoordinates();
	
	const size_t posIn = inputOutputData.size();
	inputOutputData.resize(posIn + sizeRatio*size);

	T* input = reinterpret_cast<T*>(&inputOutputData[posIn]);
	for (DiscreteType i = 0; i<size; ++i)
		input[i] = vals(i);

	InputFormatElement format;
	format.pos = posIn;
	int dim = dims.Dimension;
	format.dims.resize(dim);
	for (int i = 0; i<dim; ++i) {
		format.dims[i] = dims[i];
	}
	format.sizeRatio = sizeRatio;
	inputOutputFormat[key] = format;
}


// Public fields

bool BaseIO::HasField(KeyCRef key) const {
	if (inputOutputFormat.find(key) != inputOutputFormat.end()) return true;
	defaulted.insert(key);
	return false;
}

std::string BaseIO::GetString(KeyCRef key) const {
	const auto & val = GetInputFormat(key);
	if (val.dims.size()>0 || val.pos != std::numeric_limits<size_t>::max())
		ExceptionMacro("MathematicaIO import error : field " << key << " is not a string.");
    unused.erase(key);
	return val.str;
}

void BaseIO::SetString(KeyCRef key, std::string val) {
	SetDefined(key);
	InputFormatElement format;
	format.str = val;
	inputOutputFormat[key] = format;
}

template<typename T> std::vector<T> BaseIO::MTensorToVector(MTensor vectorMath_T) {
	double* vectorMath = libData->MTensor_getRealData(vectorMath_T);
	int vectorMathLen = libData->MTensor_getFlattenedLength(vectorMath_T);
	std::vector<T> vectorOut(vectorMath, vectorMath + vectorMathLen);
	return vectorOut;
}

template<typename T, size_t d> TraitsIO::Array<T,d> BaseIO::MTensorToArray(MTensor arrayMath_T) {
    
	// Get the MTensor data
	double* arrayMath = libData->MTensor_getRealData(arrayMath_T);// The data
	int arrayMathLen = libData->MTensor_getFlattenedLength(arrayMath_T);// Flattened length
	mint const* dimsMath = libData->MTensor_getDimensions(arrayMath_T);// Dimensions

	// Copy dimensions to the array dimensions
	TraitsIO::Array<T, d> arrayOut;
	for (int i = 0; i < d; i++)
		arrayOut.dims[i] = dimsMath[d - 1 - i];

	// Copy the MTensor data to the array
	arrayOut.resize(arrayMathLen);
	memcpy(&arrayOut[0], &arrayMath[0], arrayMathLen * sizeof(double));

	// Return the array
	return arrayOut;
}

template<typename T> MTensor BaseIO::VectorToMTensor(const std::vector<T> & vector) {
    
	// Dimension (rank) of the vector is 1
	int d = 1;

	// Get size of the vector
	mint size[d];
	size[0] = vector.size();

	// Create empty MTensor
	MTensor vectorMath_T;
	libData->MTensor_new(MType_Real, d, size, &vectorMath_T);

	// Copy data to the MTensor
	double* vectorMath;
	vectorMath = libData->MTensor_getRealData(vectorMath_T);
    static_assert(std::is_same<T,double>::value,"MTensor only accepts double values");
	memcpy(&vectorMath[0], &vector[0], size[0] * sizeof(double));

	// Return the MTensor
	return vectorMath_T;
}

template<typename T, int d> MTensor BaseIO::ArrayToMTensor(const TraitsIO::Array<T, d> & arrayIn) {
	// Get size of the array
	mint size[d];
	mint flattenedLength = 0;
	for (int i = 0; i < d; i++) {
		size[i] = arrayIn.dims[d - 1 - i];
		flattenedLength += size[i];
	}

	// Create empty MTensor
	MTensor arrayMath_T;
	libData->MTensor_new(MType_Real, d, size, &arrayMath_T);

	// Copy data to the MTensor
	double* arrayMath;
	arrayMath = libData->MTensor_getRealData(arrayMath_T);
	memcpy(&arrayMath[0], &arrayIn[0], flattenedLength * sizeof(double));

	// Return the MTensor
	return arrayMath_T;
}

template<typename T, int d> MTensor BaseIO::GetMTensor(KeyCRef key) {
    
	// Get information about the formatting of the data
	InputFormatElement format = GetInputFormat(key);
	
	// Set size of the array
    const size_t nDims = format.dims.size() + (format.sizeRatio!=1);
    
    if(nDims!=d)
        ExceptionMacro("Mathematica Export error : Field " << key << " has total depth " << nDims
                                << " distinct from requested depth " << d);
    
    mint size[d];
	for (int i = 0; i < format.dims.size(); i++) {
		size[i] = format.dims[format.dims.size() - 1 - i];
	}
	// For vector valued vectors the 1D vector is actually a 2D vector
	if (d > format.dims.size()) {
		size[d-1] = format.sizeRatio;
	}

	// Create empty MTensor
	MTensor arrayMath_T;
	libData->MTensor_new(MType_Real, d, size, &arrayMath_T);
	mint flattenedLength = libData->MTensor_getFlattenedLength(arrayMath_T);

	// Copy data to MTensor
	double* arrayMath = libData->MTensor_getRealData(arrayMath_T);;
	memcpy(&arrayMath[0], &inputOutputData[format.pos], flattenedLength * sizeof(double));
		
	// Return the MTensor
	return arrayMath_T;
}



// Constructor and destructor
BaseIO::BaseIO(WolframLibraryData libDataIn) {
	
	libData = libDataIn;
	SetString("arrayOrdering", "Reversed");
	
}

BaseIO::~BaseIO() {
}

#endif /* MathematicaIO_hxx */
