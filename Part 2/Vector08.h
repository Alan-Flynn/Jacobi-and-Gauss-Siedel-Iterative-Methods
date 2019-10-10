// File:     Vector08.h (Version W08.1)
// Author:   Niall Madden (NUI Galway) Niall.Madden@NUIGalway.ie
// Date:     Week 8 of 1819-CS319
// What:     Header file for vector class
// See also: Vector08.cpp

// In case this header files included multiple times.
#ifndef _VECTOR_H_INCLUDED
#define _VECTOR_H_INCLUDED

class Vector {
private:
  double *entries;
  unsigned int N;
public:
  Vector (unsigned int Size=2);
  Vector (const Vector &v); // Arg must be passed by value. Why?
  ~Vector(void);

  Vector &operator=(const Vector &b); // overload assignment operator
  
  unsigned int size(void) {return N;};
  double geti (unsigned int i);
  void seti (unsigned int i, double x);
  
  Vector operator+(Vector b);
  Vector operator-(Vector b); // binary
  Vector operator-(void); // unary
   
  void print(void);
  double norm(void); // Compute the 2-norm of a vector
  void zero(void); // Set entries of vector to zero.
};

#endif
