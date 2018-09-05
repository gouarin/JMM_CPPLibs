// Copyright 2017 Jean-Marie Mirebeau, University Paris-Sud, CNRS, University Paris-Saclay
// Distributed WITHOUT ANY WARRANTY. Licensed under the Apache License, Version 2.0, see http://www.apache.org/licenses/LICENSE-2.0

#ifndef RedeclareTypesMacro_h
#define RedeclareTypesMacro_h

/* Usage : 
 #define FromA(x) A:: x

 struct A {
    typedef int Integer;
    typedef double Scalar;
    cont int d=2;
 };
 
 struct B {
    Redeclare2Types(FromA,Integer,Scalar)
    Redeclare1Constant(FromA,d)
 };
 
 // For constants, in combination with templates, the linker may unfortunately require one additional line such as
 constexp typeof(A::d) B::d
 */

#define FromContainingClass(x) ContainingClass:: x
#define FromSuperclass(x) Superclass:: x
#define FromTraits(x) Traits:: x


#define Redeclare1Type(From,__a) typedef typename From(__a) __a;
#define Redeclare2Types(From,__a,...) Redeclare1Type(From,__a) Redeclare1Type(From,__VA_ARGS__)
#define Redeclare3Types(From,__a,...) Redeclare1Type(From,__a) Redeclare2Types(From,__VA_ARGS__)
#define Redeclare4Types(From,__a,...) Redeclare1Type(From,__a) Redeclare3Types(From,__VA_ARGS__)
#define Redeclare5Types(From,__a,...) Redeclare1Type(From,__a) Redeclare4Types(From,__VA_ARGS__)
#define Redeclare6Types(From,__a,...) Redeclare1Type(From,__a) Redeclare5Types(From,__VA_ARGS__)
#define Redeclare7Types(From,__a,...) Redeclare1Type(From,__a) Redeclare6Types(From,__VA_ARGS__)
#define Redeclare8Types(From,__a,...) Redeclare1Type(From,__a) Redeclare7Types(From,__VA_ARGS__)


#define Redeclare1Constant(From,__a) static constexpr decltype(From(__a)) __a = From(__a);
#define Redeclare2Constants(From,__a,...) Redeclare1Constant(From,__a) Redeclare1Constant(From,__VA_ARGS__)
#define Redeclare3Constants(From,__a,...) Redeclare1Constant(From,__a) Redeclare2Constants(From,__VA_ARGS__)
#define Redeclare4Constants(From,__a,...) Redeclare1Constant(From,__a) Redeclare3Constants(From,__VA_ARGS__)
#define Redeclare5Constants(From,__a,...) Redeclare1Constant(From,__a) Redeclare4Constants(From,__VA_ARGS__)


#endif /* RedeclareTypesMacro_h */
