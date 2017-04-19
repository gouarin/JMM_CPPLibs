// Copyright 2017 Jean-Marie Mirebeau, University Paris-Sud, CNRS, University Paris-Saclay
// Copyright 2017 Erik Bekkers, Eindhoven University of Technology, the Netherlands
// Distributed WITHOUT ANY WARRANTY. Licensed under the Apache License, Version 2.0, see http://www.apache.org/licenses/LICENSE-2.0

#ifndef MathematicaIO_h
#define MathematicaIO_h


#include "MathLink.h"
#include "WolframLibrary.h"
#include <map>
#include <fstream>
#include <sstream>
#include "Output/IO.h"

/*
A simple interface for making Mathematica data available to the c++ code.
 
 // TODO make sure MSG comes from the correct io...
*/

namespace Mathematica {

	struct BaseIO : TraitsIO {
		typedef double ScalarType;

		template<bool warn> struct Msg_ {
            std::ostringstream oss;
            template<typename T> Msg_ & operator << (const T & t){oss << t; return (*this);}
            Msg_(WolframLibraryData _libData):libData(_libData){oss << "WriteString[\"stdout\",\"";};
            Msg_():Msg_(BaseIO::latestLibData){};
            ~Msg_(){
                if(!libData) return;
                oss << "\"]";
                libData->evaluateExpression(libData,(char*)oss.str().c_str(),6,0,0);
            };
            WolframLibraryData libData;
		};

		// An adaptation of the file IO:		
		BaseIO(WolframLibraryData);
		BaseIO(const BaseIO &) = delete;
		~BaseIO();

		bool HasField(KeyCRef) const;
        bool EraseField(KeyCRef);
		std::string GetString(KeyCRef) const;
		void SetString(KeyCRef, std::string);
				
		// Conversion from Wolfram MTensors to the used IO types
		template<typename T> std::vector<T> MTensorToVector(MTensor);
		template<typename T, size_t d> TraitsIO::Array<T,d> MTensorToArray(MTensor);

		// Type conversion
		template<typename T> MTensor VectorToMTensor(const std::vector<T> &);
		template<typename T, int d> MTensor ArrayToMTensor(const TraitsIO::Array<T, d> &);

		// Directoy obtain key values in MTensor format
		template<typename T, int d> MTensor GetMTensor(KeyCRef);

		// Setting libData, which is necessary to deal with all MTensor objects
        void SetWolframLibraryData(WolframLibraryData libDataIn) { libData=libDataIn; latestLibData=libDataIn;};
	protected:
		WolframLibraryData libData = nullptr;// libData is necessary in the type conversions
        static WolframLibraryData latestLibData;
		template<typename T> std::pair<std::vector<DiscreteType>, const T*> GetDimsPtr(KeyCRef) const;
		template<typename T, size_t d, typename F> void Set(KeyCRef, DimType<d>, const F &);

		mutable std::set<KeyType> defined, unused, defaulted;

		struct InputFormatElement;
		std::map<KeyType, InputFormatElement> inputOutputFormat;
		
		const InputFormatElement & GetInputFormat(KeyCRef) const;
		void SetDefined(KeyCRef);
        void SetUsed(KeyCRef);
	};
    WolframLibraryData BaseIO::latestLibData = nullptr;


#include "MathematicaIO.hxx"

}
#endif /* MathematicaIO_h */
