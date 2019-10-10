// File: Matrix08.h (version: Week 08.1)
//   What: Header file for matrix class, including overloading = and *
// Author: Niall Madden
//   Date: Week 08 2019
// Note: Implementation for a square matrix

#ifndef _MATRIX_H_INCLUDED
#define _MATRIX_H_INCLUDED

#include "Vector08.h"

class Matrix {
private:
   double *entries;
   unsigned int N;
public:
   Matrix (unsigned int Size=2);
   Matrix (const Matrix &m); // New copy constructor
   ~Matrix(void);

   Matrix &operator=(const Matrix &B); // overload assignment operator

   unsigned int size(void) {return (N);};
   double getij (unsigned int i, unsigned int j);
   void setij (unsigned int i, unsigned int j, double x);

   Vector operator*(Vector u);
   void print(void);
};

#endif
