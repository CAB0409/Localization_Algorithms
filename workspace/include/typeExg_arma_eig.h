#ifndef TYPEEXG_ARMA_EIG_H
#define TYPEEXG_ARMA_EIG_H

// Copyright (C) 2017 Kyaw Kyaw Htike @ Ali Abdul Ghafur. All rights reserved.

/*
By default, Eigen currently supports standard floating-point types 
(float, double, std::complex<float>, std::complex<double>, long double), 
as well as all native integer types (e.g., int, unsigned int, short, etc.), and bool.
*/

#include "armadillo"
#include "Eigen/Dense"
#include <Eigen/Core>
#include <Eigen/Geometry>
#include <unsupported/Eigen/MatrixFunctions>//pow matrix
#include <vector>

#define EigenMatrix Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>
#define EIGEN_NO_DEBUG // to remove bound checking to get faster element access

using namespace Eigen;

template <typename Derived>
void print_size(const EigenBase<Derived>& b)
{
  std::cout << "size (rows, cols): " << b.size() << " (" << b.rows()
            << ", " << b.cols() << ")" << std::endl;
}

template <typename Derived>
void print_contents(const MatrixBase<Derived>& a)
{
    std::cout << "\n ------------------------------- \n";
    IOFormat CommaInitFmt(StreamPrecision, DontAlignCols, ", ", ", ", "", "", " Contents: ", ";");
    std::cout << a.format(CommaInitFmt) << "\n ------------------------------- \n";
}
/*
template <typename Derived>
void print_contents(const VectorBase<Derived>& a)
{
    std::cout << "\n ------------------------------- \n";
    IOFormat CommaInitFmt(StreamPrecision, DontAlignCols, ", ", ", ", "", "", " Contents: ", ";");
    std::cout << a.format(CommaInitFmt) << "\n ------------------------------- \n";
}*/

/////////// *************** ///////////////////
/////////// *************** ///////////////////
/////////// CORRECTNESS TEST PASSED ///////////
/////////// *************** ///////////////////
/////////// *************** ///////////////////

// works for 2D matrices (real numbers, not complex)
// assumes that underlying matrix data is stored in column major order.
template <typename T>
void eigen2arma(const EigenMatrix& matIn, arma::Mat<T> &matOut)
{	
	int nrows = matIn.rows();
	int ncols = matIn.cols();
	
	matOut.set_size(nrows, ncols);
	
	T * ptr = matOut.memptr();
	unsigned long count = 0;
	
	for (int j = 0; j < ncols; j++)
		for (int i = 0; i < nrows; i++)
			ptr[count++] = matIn(i,j);		
}

// works for 3D matrices (real numbers, not complex)
// a 3D matrix can be thought of as a 2D matrix with any number of channels
// assumes that underlying matrix data is stored in column major order.
template <typename T>
void eigen2arma(const std::vector<EigenMatrix> &matIn, arma::Cube<T> &matOut)
{
	int nrows = matIn[0].rows();
	int ncols = matIn[0].cols();
	int nchannels = matIn.size();
	
	matOut.set_size(nrows, ncols, nchannels);
	
	T * ptr = matOut.memptr();
	unsigned long count = 0;
	
	EigenMatrix temp;
	
	for (int k = 0; k < nchannels; k++)
	{
		temp = matIn[k];
		for (int j = 0; j < ncols; j++)
		{
			for (int i = 0; i < nrows; i++)
			{
				ptr[count++] = temp(i,j);				
			}				
		}	
	}
}


// works for 2D matrices (real numbers, not complex)
// assumes that underlying matrix data is stored in column major order.
template <typename T>
void arma2eigen(const arma::Mat<T> &matIn, EigenMatrix &matOut)
{		
	int nrows = matIn.n_rows;
	int ncols = matIn.n_cols;
	
	matOut.setZero(nrows, ncols);
	
	T *dst_pointer = (T*)matOut.data();
	unsigned long count = 0;

	for (int j = 0; j < ncols; j++)
		for (int i = 0; i < nrows; i++)
			dst_pointer[count++] = matIn.at(i,j);
}

// works for 3D matrices (real numbers, not complex)
// a 3D matrix can be thought of as a 2D matrix with any number of channels
// assumes that underlying matrix data is stored in column major order.
template <typename T>
void arma2eigen(const arma::Cube<T> matIn, std::vector<EigenMatrix> &matOut)
{
	
	int nrows = matIn.n_rows;
	int ncols = matIn.n_cols;
	int nchannels = matIn.n_slices;

	matOut.resize(nchannels);	
	
	for (int k = 0; k < nchannels; k++)
	{
		matOut[k].setZero(nrows, ncols);
		T *dst_pointer = (T*)matOut[k].data();
		unsigned long count = 0;
		for (int j = 0; j < ncols; j++)
			for (int i = 0; i < nrows; i++)
				dst_pointer[count++] = matIn.at(i,j,k);		
	}		

}


#undef EigenMatrix
#undef EIGEN_NO_DEBUG

#endif