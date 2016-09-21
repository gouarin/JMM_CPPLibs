//
//  Macros.h
//  OptimalTransport
//
//  Created by Jean-Marie Mirebeau on 12/03/2014.
//  Copyright (c) 2014 Jean-Marie Mirebeau. All rights reserved.
//

#ifndef OptimalTransport_Macros_h
#define OptimalTransport_Macros_h

#include "PrintContainers.h"
#include "PrintAsMathematicaTable.h"
#include "EnumToString.h"

#include "sys/time.h"
#include "time.h"

#define ExportVarArrow(name) \
<< '"' << #name << '"' << " -> "  <<  name  << ","

#define ExportEnumArrow(name) \
<< '"' << #name << '"' << " -> " << '"' << enumToString(name) << '"' << ","

#define ExportArrayArrow(name) \
<< '"' << #name << '"' << " -> "; print_container(os, name); os << ","

#define ExportArrayOfArraysArrow(indices,name) \
<< '"' << #name << '"' << " -> "; PrintArrayOfArrays(os, indices, name); os << ","

// Note that CPU time adds up the contribution of all threads, unlike not physical time.
#define ExecuteAndPushCPUTime(instructions, array) \
{ \
clock_t top = -clock(); \
instructions; \
top+= clock(); \
array.push_back(top/double(CLOCKS_PER_SEC)); \
}

#define ExecuteAndPushPhysicalTime(instructions,array) \
{ \
timeval start, end; \
gettimeofday(&start, nullptr); \
instructions; \
gettimeofday(&end, nullptr); \
array.push_back((end.tv_sec-start.tv_sec)+(end.tv_usec-start.tv_usec)/1000000.); \
}

#endif
