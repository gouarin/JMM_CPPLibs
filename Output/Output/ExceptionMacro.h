// Copyright 2017 Jean-Marie Mirebeau, University Paris-Sud, CNRS, University Paris-Saclay
// Distributed WITHOUT ANY WARRANTY. Licensed under the Apache License, Version 2.0, see http://www.apache.org/licenses/LICENSE-2.0

#ifndef ExceptionMacro_h
#define ExceptionMacro_h

#include <sstream>

#define ExceptionMacro(msg) \
{ \
std::ostringstream message; \
message << msg << "\n" << \
"File: " << __FILE__ << ", Line: " << __LINE__ << "\n"; \
throw std::logic_error(message.str()); \
}

#endif /* ExceptionMacro_h */
