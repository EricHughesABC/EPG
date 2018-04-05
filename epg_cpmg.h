#pragma once

#include <iostream>     // std::cout
#include <complex>      // std::complex
#include <array>
#include <new>

void print_complex_matrix(std::complex<double> *mmm, int nrows, int ncols);
void print_double_matrix(double *mmm, int nrows, int ncols);
std::complex<double> *rf_pulse_matrix(double alpha, double phi);
int mult_pulse_FpFmZ(std::complex<double> *rfpulse, std::complex<double> *FpFmZ, int num_cols);
int mult_ee_FpFmZ(double *ee, std::complex<double> *FpFmZ, int num_cols);
void epg_grad(std::complex<double> *FpFmZ, int num_cols);
void cpmg_epg(double *signal, int Nechos, double rf_90, double rf_180, double T1, double T2, double T);
void cpmg_epg_b1(double *signal, int Nechos, double rf_90, double rf_180, double T1, double T2, double T, double B1scale);



