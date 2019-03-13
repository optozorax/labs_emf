#pragma once

#include <string>
#include <vector>
#include "matrix.h"
#include "common.h"

//-----------------------------------------------------------------------------
class Vector
{
public:
	Vector();
	Vector(int size, real fill = 0);
	Vector(const Matrix& a);

	void loadFromFile(std::string fileName);
	void saveToFile(std::string fileName) const;

	void save(std::ostream& out) const;
	void load(std::istream& in);

	void generate(int n, int min, int max);
	void generate(int n);

	//-------------------------------------------------------------------------
	void resize(int n, real fill = 0);

	void negate(void);

	void zero(void);

	void toDenseMatrix(Matrix& dense, bool isVertical) const;

	//-------------------------------------------------------------------------
	int size(void) const;
	
	real& operator()(int i);
	const real& operator()(int i) const;

private:
	std::vector<real> mas;
};

//-----------------------------------------------------------------------------
sumreal sumAllElementsAbs(const Vector& a);
bool sum(const Vector& a, const Vector& b, Vector& result);
bool mul(const Matrix& a, const Vector& b, Vector& result);
real calcNorm(const Vector& a);