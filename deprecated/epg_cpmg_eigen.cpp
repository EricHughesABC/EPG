// EigenConsoleApplication1.cpp : Defines the entry point for the console application.
//


//#include "stdafx.h"
#include <ctime>
#include <Eigen/Dense>
#include <iostream>
#include <complex>      // std::complex

//using Eigen::MatrixXd;
//using Eigen::Matrix3cd;

void rf_pulse_matrix(Eigen::Matrix3cd *rfpulse, double alpha, double phi)
{
	const double PI = 3.1415926535897932384626433832795;

	double alpha_half = alpha / 2.0;
	double c_a_h = cos(alpha_half);
	double c_a_h2 = c_a_h * c_a_h;
	double s_a_h = sin(alpha_half);
	double s_a_h2 = s_a_h * s_a_h;
	double c_a = cos(alpha);
	double s_a = sin(alpha);

	std::complex<double> i_phi = std::complex<double>(0.0, phi);
	std::complex<double> i_half = std::complex<double>(0.0, 0.5);

	(*rfpulse)(0, 0) = c_a_h2;
	(*rfpulse)(0, 1) = exp(i_phi*2.0)*s_a_h2;
	(*rfpulse)(0, 2) = -1.0*std::complex<double>(0.0, 1.0)*exp(i_phi)*s_a;

	(*rfpulse)(1, 0) = exp(-2.0 * i_phi)*s_a_h2;
	(*rfpulse)(1, 1) = c_a_h2;
	(*rfpulse)(1, 2) = std::complex<double>(0.0, 1.0) * exp(-1.0*i_phi)*s_a;

	(*rfpulse)(2, 0) = -1.0*i_half*exp(-1.0*i_phi)*s_a;
	(*rfpulse)(2, 1) = 1.0*i_half*exp(1.0*i_phi)*s_a;
	(*rfpulse)(2, 2) = c_a;

	//std::cout << *rfpulse << std::endl;
}









void epg_grad(Eigen::MatrixXcd *FpFmZ, int num_cols)
{

	for (int i = 0; i < num_cols - 1; i++)
	{
		//std::cout << i << "\t" << num_cols - 1 - i << "\t" << num_cols - 1 - i - 1 << "\t" << (*FpFmZ)(0,i) << std::endl;
		(*FpFmZ)(0, num_cols - 1 - i) = (*FpFmZ)(0, num_cols - 1 - i - 1);
	}


	(*FpFmZ)(0, 0) = std::complex<double>(0.0, 0.0);

	for (int i = 0; i < num_cols - 1; i++)
		(*FpFmZ)(1, i) = (*FpFmZ)(1, i + 1);

	(*FpFmZ)(1, num_cols - 1) = std::complex<double>(0.0, 0.0);

	(*FpFmZ)(0, 0) = conj((*FpFmZ)(1, 0));
}


void  cpmg_epg(double *signal, int Nechos, double rf_180, double T1, double T2, double Techo)
{
	const double PI = 3.1415926535897932384626433832795;
	int    num_cols = 2 * Nechos;

	Eigen::MatrixXcd FpFmZ(3, 34);
	Eigen::Matrix3cd rfpulse90y = Eigen::Matrix3cd::Zero();
	Eigen::Matrix3cd rfpulse180x = Eigen::Matrix3cd::Zero();
	Eigen::Matrix3d ee = Eigen::Matrix3d::Zero();
	//Eigen::MatrixXd signal(1, 17);




	FpFmZ = FpFmZ.Zero(3, 34);



	// ****************************************
	// Initialize diagonal of relaxation matrix
	// ****************************************

	ee(0, 0) = exp(-Techo / 2.0 / T2);
	ee(1, 1) = exp(-Techo / 2.0 / T2);
	ee(2, 2) = exp(-Techo / 2.0 / T1);


	// ********************************
	// Set Z magnetization to 1.0 +0.0j
	// ********************************

	FpFmZ(0, 0) = 0.0;
	FpFmZ(1, 0) = 0.0;
	FpFmZ(2, 0) = 1.0;

	// *********************************
	// Create 90(y) rotation matrix
	// *********************************

	rf_pulse_matrix(&rfpulse90y, PI / 2.0, PI / 2.0);

	// *********************************
	// Create 180(x) rotation matrix
	// *********************************

	rf_pulse_matrix(&rfpulse180x, rf_180*PI / 180.0, 0.0);

	// *********************************
	// Apply 90 degree pulse
	// *********************************

	FpFmZ = rfpulse90y * FpFmZ;

	// ****************************************
	// Apply 180 pulses of CPMG sequence
	// ****************************************

	for (int i = 0; i < Nechos; i++)
	{

		FpFmZ = ee * FpFmZ;
		FpFmZ(2, 0) = FpFmZ(2, 0) + 1.0 - ee(2, 2);
		epg_grad(&FpFmZ, 34);

		FpFmZ = rfpulse180x * FpFmZ;

		FpFmZ = ee * FpFmZ;
		FpFmZ(2, 0) = FpFmZ(2, 0) + 1.0 - ee(2, 2);
		epg_grad(&FpFmZ, 34);

		signal[i] = real(FpFmZ(0, 0));
	}


}



//int main()
//{
//	const double PI = 3.1415926535897932384626433832795;
//	std::clock_t    start;
//
//	double T1 = 3000.0;
//	double T2 = 50.0;
//	double Techo = 10.0;
//	int    Nechos = 17;
//	int    num_cols = 2 * Nechos;
//
//	Eigen::MatrixXcd FpFmZ(3, num_cols);
//	Eigen::Matrix3cd rfpulse90y = Eigen::Matrix3cd::Zero();
//	Eigen::Matrix3cd rfpulse180x = Eigen::Matrix3cd::Zero();
//	Eigen::Matrix3d ee = Eigen::Matrix3d::Zero();
//
//	double *signal = new  double [Nechos] { 0.0 };
//
//	start = std::clock();
//
//	for (int j = 0; j < 10000; j++)
//		cpmg_epg(signal, Nechos, PI, T1, T2, Techo);
//
//	std::cout << "Time: " << (std::clock() - start) / (double)(CLOCKS_PER_SEC / 1000) << " ms" << std::endl;
//
//	for (int i = 0; i<Nechos; i++)
//		std::cout << signal[i] << std::endl;
//
//
//
//}
