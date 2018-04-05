// epg_cpmg_1.cpp
//C++ module for calculating EPG (Extended Phase Graph) signals
//Based on Matlab code by Brian Hargreaves
//Diffusion effects not implemented
//
// All 2-D matrices coded as a 1-D array.  Elements accessed by calculating an offset:
// ie  for a 3x3 rotation matrix rrr
// rrr[1][2] = val is given by the following
//
// rrr[1*num_cols+2] = val
//
// where num_cols = 3, the number of columns in the matix.
//
// The reason for using this approach is for speed, we have looked at the library eigen and the code was half the speed
// as this approach.

#include "epg_cpmg.h"



void print_complex_matrix(std::complex<double> *mmm, int nrows, int ncols)
{
	// mmm is a 2-D complex matrix of size nrows by ncols coded as a 1-D array of length nrows x ncols

	for (int r = 0; r < nrows; r++) {
		for (int c = 0; c < ncols; c++) {
			std::cout << r << " " << c << " " << mmm[r*ncols + c] << " ";
		}
		std::cout << "\n";
	}
	std::cout << "\n";
}


void print_double_matrix(double *mmm, int nrows, int ncols)
{
	// mmm is a 2-D double matrix of size nrows by ncols coded as a 1-D array of length nrows x ncols

	for (int r = 0; r < nrows; r++) {
		for (int c = 0; c < ncols; c++) {
			//FpFmZ[r*ncols + c] = std::complex<double>(r, c);
			std::cout << r << " " << c << " (" << mmm[r*ncols + c] << ") ";
		}
		std::cout << "\n";
	}
	std::cout << "\n";
}


 std::complex<double> *rf_pulse_matrix(double alpha, double phi)
{
	// Calculates the elements of a 3x3 rotation matrix based on the paper
	// by Mathias Weigel, Extended phase graphs: Dephasing, RF pulses, and echoes - pure and simple
	// J. Magn. Reson. Imaging 2015;41:266–295
	// alpha is angle of rotation in radians
	// phi is the phase angle in radians
	// The 2-D 3x3 matrix is represented as a 1-D array of 9 elements

	std::complex<double> *rfpulse = new  std::complex<double>[9];

	double alpha_half = alpha / 2.0;
	double c_a_h = cos(alpha_half);
	double c_a_h2 = c_a_h * c_a_h;
	double s_a_h = sin(alpha_half);
	double s_a_h2 = s_a_h * s_a_h;
	double c_a = cos(alpha);
	double s_a = sin(alpha);
	std::complex<double> i_phi = std::complex<double>(0.0, phi);
	std::complex<double> i_half = std::complex<double>(0.0, 0.5);

	rfpulse[0] = c_a_h2;
	rfpulse[1] = exp(2.0*i_phi)*s_a_h2;
	rfpulse[2] = -1.0*std::complex<double>(0.0, 1.0)*exp(i_phi)*s_a;

	rfpulse[3 + 0] = exp(-2.0 * i_phi)*s_a_h2;
	rfpulse[3 + 1] = c_a_h2;
	rfpulse[3 + 2] = std::complex<double>(0.0, 1.0) * exp(-1.0*i_phi)*s_a;

	rfpulse[6 + 0] = -1.0*i_half*exp(-1.0*i_phi)*s_a;
	rfpulse[6 + 1] = 1.0*i_half*exp(1.0*i_phi)*s_a;
	rfpulse[6 + 2] = c_a;

	return rfpulse;
}


 int mult_pulse_FpFmZ(std::complex<double> *rfpulse, std::complex<double> *FpFmZ, int num_cols_FpFmZ)
{

	// This function does the matrix multiplication of rfpulse x FpFmZ
	// The result is returned in FpFmZ
	// rfpulse is a 2-D matrix of 3 rows and 3 columns represented as a 1-D array
	// FpFmZ is a 2-Dmatrix of 3 rows and num_cols_FpFmZ columns represented as a 1-D array

	const int num_rows=3;
	const int num_cols_rfpulse = 3;

	std::complex<double> *result = new  std::complex<double>[num_rows * num_cols_FpFmZ];
	std::complex<double> sum = std::complex<double>(0.0, 0.0);

	for (int r = 0; r < num_rows; r++) // row of rfpulse and result
	{
		for (int c = 0; c < num_cols_FpFmZ; c++)   //columns of result and FpFmZ
		{
			for (int k = 0; k < num_rows; k++) // row index of FpFmZ
			{
				sum = sum + rfpulse[r * num_cols_rfpulse + k] * FpFmZ[k*num_cols_FpFmZ + c];
			}
			result[r*num_cols_FpFmZ + c] = sum;
			sum = std::complex<double>(0.0, 0.0);  // reset sum to zero
		}
	}
	for (int i = 0; i < num_rows * num_cols_FpFmZ; i++)
		FpFmZ[i] = result[i];

	delete[] result;
	return 0;
}


 int mult_ee_FpFmZ(double *ee, std::complex<double> *FpFmZ, int num_cols_FpFmZ)
{
	// This function does the matrix multiplication of ee x FpFmZ
	// The result is returned in FpFmZ
	// ee is a 2-D matrix of 3 rows and 3 columns represented as a 1-D array
	// FpFmZ is a 2-Dmatrix of 3 rows and num_cols_FpFmZ columns represented as a 1-D array

	const int num_rows=3;
	const int num_cols_rfpulse = 3;

	std::complex<double> *result = new  std::complex<double>[3 * num_cols_FpFmZ];
	std::complex<double> sum = std::complex<double>(0.0, 0.0);

	for (int r = 0; r < num_rows; r++)
	{
		for (int c = 0; c < num_cols_FpFmZ; c++)
		{
			for (int k = 0;  k < num_rows;  k++)
			{
				sum = sum + ee[r * num_cols_rfpulse + k] * FpFmZ[k*num_cols_FpFmZ + c];
			}
			result[r*num_cols_FpFmZ + c] = sum;
			sum = std::complex<double>(0.0, 0.0);
		}
	}
	for (int i = 0; i < num_rows * num_cols_FpFmZ; i++)
		FpFmZ[i] = result[i];

	delete[] result;
	return 0;
}


 void epg_grad(std::complex<double> *FpFmZ, int num_cols)
{

	for (int i = 0; i<num_cols - 1; i++)
		FpFmZ[num_cols - 1 - i] = FpFmZ[num_cols - 1 - i - 1];

	FpFmZ[0] = std::complex<double>(0.0, 0.0);

	for (int i = 0; i<num_cols - 1; i++)
		FpFmZ[num_cols + i] = FpFmZ[num_cols + i + 1];

	FpFmZ[num_cols + num_cols - 1] = std::complex<double>(0.0, 0.0);

	FpFmZ[0] = conj(FpFmZ[num_cols]);
}


 void cpmg_epg(double *signal, int Nechos, double rf_90, double rf_180, double T1, double T2, double T)
{
	// signal : final calculated EPG signal, Nechos in length
	// Nechos : number of echos to calculate
	// rf_180 : rf pulse value for cpmg given in degrees
	// T1 : value for T1 relaxation time given in ms
	// T2 : value for T2 relaxation time given in ms
	// T  : value for CPMG echo spacing given in ms

	const double PI = 3.1415926535897932384626433832795;

	int nrows = 3;  // Fp, Fm, Z
	int ncols = 2 * Nechos;

	std::complex<double> *FpFmZ = new  std::complex<double>[nrows * ncols];
	double *ee = new double[9]{ 0.0 };
	std::complex<double> *r90;
	std::complex<double> *r180;

	// ****************************************
	// Initialize diagonal of relaxation matrix
	// ****************************************

	ee[0] = exp(-T /2.0/ T2);
	ee[4] = exp(-T /2.0/ T2);
	ee[8] = exp(-T /2.0/ T1);

	// ********************************
	// Set Z magnetization to 1.0 +0.0j
	// ********************************

	FpFmZ[2 * ncols + 0] = std::complex<double>(1.0, 0.0);

	// *********************************
	// Create 90(y) rotation matrix
	// *********************************

	r90 = rf_pulse_matrix(rf_90*PI/180.0, PI / 2.0);

	// *********************************
	// Create 180(x) rotation matrix
	// *********************************

	r180 = rf_pulse_matrix(rf_180*PI/180.0, 0.0);

	// *********************************
	// Apply 90 degree pulse
	// *********************************

	mult_pulse_FpFmZ(r90, FpFmZ, ncols);

	// ****************************************
	// Apply 180 pulses of CPMG sequence
	// ****************************************

	for (int i = 0; i<Nechos; i++)
	{
		//Apply relaxation for  first half of half echo time
		mult_ee_FpFmZ(ee, FpFmZ, ncols);
		FpFmZ[2 * ncols ] = FpFmZ[2 * ncols ] + 1.0 - ee[8];


		//Apply gradient
		epg_grad(FpFmZ, ncols);

		// Apply effective 180 degree pulse
		mult_pulse_FpFmZ(r180, FpFmZ, ncols);

		//Apply relaxation for  second half of half echo time
		mult_ee_FpFmZ(ee, FpFmZ, ncols);
		FpFmZ[2 * ncols ] = FpFmZ[2 * ncols ] + 1.0 - ee[8];

		//Apply gradient
		epg_grad(FpFmZ, ncols);

		//Save signal
		signal[i] = real(FpFmZ[0]);
	}




	delete[] FpFmZ;

	delete[] r90;
	delete[] r180;
	delete[] ee;

	//return signal;
}




 void cpmg_epg_b1(double *signal, int Nechos, double rf_90, double rf_180, double T1, double T2, double T, double B1scale)
{
	// signal : final calculated EPG signal, Nechos in length
	// Nechos : number of echos to calculate
	// rf_180 : rf pulse value for cpmg given in degrees
	// T1 : value for T1 relaxation time given in ms
	// T2 : value for T2 relaxation time given in ms
	// T  : value for CPMG echo spacing given in ms

	const double PI = 3.1415926535897932384626433832795;

	int nrows = 3;  // Fp, Fm, Z
	int ncols = 2 * Nechos;

	std::complex<double> *FpFmZ = new  std::complex<double>[nrows * ncols];
	double *ee = new double[9]{ 0.0 };
	std::complex<double> *r90;
	std::complex<double> *r180;

	// ****************************************
	// Initialize diagonal of relaxation matrix
	// ****************************************

	ee[0] = exp(-T /2.0/ T2);
	ee[4] = exp(-T /2.0/ T2);
	ee[8] = exp(-T /2.0/ T1);

	// ********************************
	// Set Z magnetization to 1.0 +0.0j
	// ********************************

	FpFmZ[2 * ncols + 0] = std::complex<double>(1.0, 0.0);

	// *********************************
	// Create 90(y) rotation matrix
	// *********************************

	r90 = rf_pulse_matrix(B1scale*rf_90*PI/180.0, PI / 2.0);

	// *********************************
	// Create 180(x) rotation matrix
	// *********************************

	r180 = rf_pulse_matrix(B1scale*rf_180*PI/180.0, 0.0);

	// *********************************
	// Apply 90 degree pulse
	// *********************************

	mult_pulse_FpFmZ(r90, FpFmZ, ncols);

	// ****************************************
	// Apply 180 pulses of CPMG sequence
	// ****************************************

	for (int i = 0; i<Nechos; i++)
	{
		//Apply relaxation for  first half of half echo time
		mult_ee_FpFmZ(ee, FpFmZ, ncols);
		FpFmZ[2 * ncols ] = FpFmZ[2 * ncols ] + 1.0 - ee[8];


		//Apply gradient
		epg_grad(FpFmZ, ncols);

		// Apply effective 180 degree pulse
		mult_pulse_FpFmZ(r180, FpFmZ, ncols);

		//Apply relaxation for  second half of half echo time
		mult_ee_FpFmZ(ee, FpFmZ, ncols);
		FpFmZ[2 * ncols ] = FpFmZ[2 * ncols ] + 1.0 - ee[8];

		//Apply gradient
		epg_grad(FpFmZ, ncols);

		//Save signal
		signal[i] = real(FpFmZ[0]);
	}




	delete[] FpFmZ;

	delete[] r90;
	delete[] r180;
	delete[] ee;

	//return signal;
}


//int main()
//{
//
//
//	const double PI = 3.1415926535897932384626433832795;
//
//	double* result;
//
//	int Nechos = 17;
//
//	result = cpmg_epg(Nechos, 120.*PI/180.0, 3000.0, 50.0, 10.0/2.0);
//
//	std::cout << "\nresult\n";
//
//	for (int i = 0; i < Nechos; i++)
//		std::cout << result[i] << ",\n";
//
//	delete[] result;
//
//	return 0;
//}

