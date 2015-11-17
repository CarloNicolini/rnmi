cmake_minimum_required(VERSION 2.8 FATAL_ERROR)
project(rnmi)

add_executable(rnmi rnmi.cpp zrg.cpp nmi.h)
add_executable(nmi nmi.cpp zrg.cpp nmi.h)

if(MATLAB_SUPPORT)
    add_mex(rnmi_mx rnmi_mx.cpp zrg.cpp nmi.h)
endif(MATLAB_SUPPORT)

if(OCTAVE_SUPPORT)
    octave_add_oct(rnmi_oct SOURCES rnmi_mx.cpp zrg.cpp nmi.h "mex")
endif(OCTAVE_SUPPORT)
