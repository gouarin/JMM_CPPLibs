// Copyright 2017 Jean-Marie Mirebeau, University Paris-Sud, CNRS, University Paris-Saclay
// Copyright 2017 Erik Bekkers, Eindhoven University of Technology, the Netherlands
// Distributed WITHOUT ANY WARRANTY. Licensed under the Apache License, Version 2.0, see http://www.apache.org/licenses/LICENSE-2.0

#ifndef MathematicaIOExternC_h
#define MathematicaIOExternC_h


#define MathematicaTryCatch(code) \
try{code; \
return LIBRARY_NO_ERROR; \
} catch (const std::logic_error & e) { \
IO::WarnMsg() << "Exception caught. " << e.what(); \
return LIBRARY_FUNCTION_ERROR; \
}




// All functions that interact with Wolfram Mathematica are included in -extern "C"
extern "C"
{//Begin of Extern "C"

	  /* ***************************************************************************** */
	  /* ********************* All the set/get functions ***************************** */
	  /* ***************************************************************************** */

	  // GetScalar/SetScalar
	DLLEXPORT int GetScalar(WolframLibraryData libData, mint Argc, MArgument *Args, MArgument Res){
        MathematicaTryCatch(io.SetWolframLibraryData(libData);
                            const std::string key = MArgument_getUTF8String(Args[0]);
                            MArgument_setReal(Res, io.Get<double>(key));
                            )
    }
    
	DLLEXPORT int SetScalar(WolframLibraryData libData, mint Argc, MArgument *Args, MArgument Res){
        MathematicaTryCatch(io.SetWolframLibraryData(libData);
                            const std::string key = MArgument_getUTF8String(Args[0]);
                            const double val = MArgument_getReal(Args[1]);
                            io.Set<double>(key, val);
                            )
    }

	// GetString/SetString
	DLLEXPORT int GetString(WolframLibraryData libData, mint Argc, MArgument *Args, MArgument Res)
	{
        MathematicaTryCatch(io.SetWolframLibraryData(libData);
                            const std::string key = MArgument_getUTF8String(Args[0]);
                            const std::string stringSTD = io.GetString(key);
                            MArgument_setUTF8String(Res, (char*)stringSTD.c_str());
        )
//			char * string = new char[stringSTD.size() + 1];
//			std::copy(stringSTD.begin(), stringSTD.end(), string);
	}
	DLLEXPORT int SetString(WolframLibraryData libData, mint Argc, MArgument *Args, MArgument Res){
        MathematicaTryCatch(io.SetWolframLibraryData(libData);
                            const std::string key = MArgument_getUTF8String(Args[0]);
                            const std::string val = MArgument_getUTF8String(Args[1]);
                            io.SetString(key, val);
                            )
    }

	// GetVector/SetVector
	DLLEXPORT int GetVector(WolframLibraryData libData, mint Argc, MArgument *Args, MArgument Res){
        MathematicaTryCatch(io.SetWolframLibraryData(libData);
                            const std::string key = MArgument_getUTF8String(Args[0]);
                            MArgument_setMTensor(Res, io.VectorToMTensor<double>(io.GetVector<double>(key)));
                            )
	}
    
	DLLEXPORT int SetVector(WolframLibraryData libData, mint Argc, MArgument *Args, MArgument Res)
	{
        MathematicaTryCatch(io.SetWolframLibraryData(libData);
                            const std::string key = MArgument_getUTF8String(Args[0]);
                            MTensor val = MArgument_getMTensor(Args[1]);
                            io.SetVector<double>(key, io.MTensorToVector<double>(val));
                            )
	}

	// GetArray/SetArray
	DLLEXPORT int GetArray(WolframLibraryData libData, mint Argc, MArgument *Args, MArgument Res){
        MathematicaTryCatch(io.SetWolframLibraryData(libData);
                            const std::string key = MArgument_getUTF8String(Args[0]);
                            const int d = MArgument_getInteger(Args[1]); // Array dimension
                            switch (d) {
                                case 1:MArgument_setMTensor(Res, (io.GetMTensor<double, 1>(key)));break;
                                case 2:MArgument_setMTensor(Res, (io.GetMTensor<double, 2>(key)));break;
                                case 3:MArgument_setMTensor(Res, (io.GetMTensor<double, 3>(key)));break;
                                case 4:MArgument_setMTensor(Res, (io.GetMTensor<double, 4>(key)));break;
                                case 5:MArgument_setMTensor(Res, (io.GetMTensor<double, 5>(key)));break;
                                case 6:MArgument_setMTensor(Res, (io.GetMTensor<double, 6>(key)));break;
                                default:
                                    ExceptionMacro("GetArray error : unsupported dimension");
                            }
                            )
	}
    DLLEXPORT int SetArray(WolframLibraryData libData, mint Argc, MArgument *Args, MArgument Res){
        /*MathematicaTryCatch( // Macro fails to compile for unknown reason.
         io.SetWolframLibraryData(libData);
         const std::string key = MArgument_getUTF8String(Args[0]);
         const MTensor val = MArgument_getMTensor(Args[1]);
         int d = libData->MTensor_getRank(val); // Array dimension
         if(libData->MTensor_getFlattenedLength(val)==0){
         ExceptionMacro("SetArray Error: Empty tensor for key " << key << "\n");}
         switch (d) {
         case 1: io.SetArray<double, 1>(key, io.MTensorToArray<double, 1>(val)); break;
         case 2: io.SetArray<double, 2>(key, io.MTensorToArray<double, 2>(val)); break;
         case 3: io.SetArray<double, 3>(key, io.MTensorToArray<double, 3>(val)); break;
         case 4: io.SetArray<double, 4>(key, io.MTensorToArray<double, 4>(val)); break;
         case 5: io.SetArray<double, 5>(key, io.MTensorToArray<double, 5>(val)); break;
         case 6: io.SetArray<double, 6>(key, io.MTensorToArray<double, 6>(val)); break;
         default:
         ExceptionMacro("GetArray error : unsupported dimension");
         }
         )*/
        
        try{
            io.SetWolframLibraryData(libData);
            const std::string key = MArgument_getUTF8String(Args[0]);
            const MTensor val = MArgument_getMTensor(Args[1]);
            int d = libData->MTensor_getRank(val); // Array dimension
            if(libData->MTensor_getFlattenedLength(val)==0){
                ExceptionMacro("SetArray Error: Empty tensor for key " << key << "\n");}
            switch (d) {
                case 1: io.SetArray<double, 1>(key, io.MTensorToArray<double, 1>(val)); break;
                case 2: io.SetArray<double, 2>(key, io.MTensorToArray<double, 2>(val)); break;
                case 3: io.SetArray<double, 3>(key, io.MTensorToArray<double, 3>(val)); break;
                case 4: io.SetArray<double, 4>(key, io.MTensorToArray<double, 4>(val)); break;
                case 5: io.SetArray<double, 5>(key, io.MTensorToArray<double, 5>(val)); break;
                case 6: io.SetArray<double, 6>(key, io.MTensorToArray<double, 6>(val)); break;
                default:
                    ExceptionMacro("GetArray error : unsupported dimension");
            }
            return LIBRARY_NO_ERROR;
        } catch (std::logic_error & e){
            IO::WarnMsg() << "Exception caught. " << e.what(); \
            return LIBRARY_FUNCTION_ERROR;
        }
    }
}//End of Extern "C"


#endif /* MathematicaIOExternC_h */
