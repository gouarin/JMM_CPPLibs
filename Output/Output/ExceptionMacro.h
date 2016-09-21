//
//  ExceptionMacro.h
//  LiftedFastMarching_MathematicaInterface
//
//  Created by Jean-Marie Mirebeau on 14/09/2016.
//
//

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
