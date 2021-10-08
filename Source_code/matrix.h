#pragma once

//Structure to support simple matrix format
struct matrix2D
{
	//data
	double **data;
	int nrow;
	int ncol;
	
	//fuctions
	matrix2D(int M, int N);
	~matrix2D();

	void setInitvalue(double initval);
	double &operator()(int i, int j);
	const double &operator()(int i, int j) const;
};