/* File: epg.i */
%module epg
%feature("autodoc", "3");

%{
#define SWIG_FILE_WITH_INIT
#include "epg_cpmg.h"
%}

%include "numpy.i"
%init %{
import_array();
%}

//%apply ( double* ARGOUT_ARRAY1, int DIM1 ){(double* signal, int Nechos)};
%apply ( double* INPLACE_ARRAY1, int DIM1 ){(double* signal, int Nechos)};
void cpmg_epg(double* signal, int Nechos, double rf_90, double rf_180, double T1, double T2, double Techo);

//%apply ( double* ARGOUT_ARRAY1, int DIM1 ){(double* signal, int Nechos)};
%apply ( double* INPLACE_ARRAY1, int DIM1 ){(double* signal, int Nechos)};
void cpmg_epg_b1(double* signal, int Nechos, double rf_90, double rf_180, double T1, double T2, double Techo, double B1scale);




