#include "matrix.h"

//Functions of the matrix structure


//constructor
matrix2D::matrix2D(int M, int N)
{
	data = new double*[M];
	for(int i=0; i < M; ++i)
		data[i] = new double[N];
	nrow = M;
	ncol = N;
}

//destructor
matrix2D::~matrix2D()
{
	for(int i=0; i < nrow; ++i)
		delete [] data[i];
	delete [] data;
}

//Functions for setting common values for every element in the matrix
void matrix2D::setInitvalue(double initval)
{
	for(int i=0; i < nrow; ++i)
		for(int j=0; j < ncol; ++j)
			data[i][j]=  initval;
}

//Functions for accessing matrix elements
double & matrix2D::operator()(int i, int j)
{
	return data[i][j];
}
const double & matrix2D::operator()(int i, int j) const
{
	return data[i][j];	
}
