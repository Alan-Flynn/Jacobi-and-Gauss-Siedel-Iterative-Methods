//Name: Alan Flynn
//ID: 16331333

#include <iostream>
#include <stdlib.h>

#include "Vector.h"
#include "Matrix.h"
#include <cmath>


//  c=alpha*a + beta*b where a,b are Vectors; alpha, beta are scalars
void VecAdd (Vector &c, Vector &a, Vector &b,
             double alpha=1.0, double beta=1.0);

//  Compute a matrix-Vector product (v=A*u)
void MatVec(Matrix &A, Vector &u, Vector &v);

// Use Gauss' method to solve a least-squares problem
void Gauss(Matrix &A, Vector &b, Vector &v,
	    unsigned int &count, double tol, int max_iteration, bool &max_reachedG);

void Jacobi(Matrix &A, Vector &b, Vector &v,
	    unsigned int &count, double tol, int max_iteration, bool &max_reachedJ);

int main(void )
{
  unsigned int N=3;

  Matrix A(N);
  Vector x(N), b(N);

  // Set up the linear system to solve, which is
  // 9*x_1 + 3*x_2 + 3*x_3 = 15
  // 3*x_1 + 9*x_2 + 3*x_3 = 15
  // 3*x_1 + 3*x_2 + 9*x_3 = 15
  for (unsigned int i=0; i<N; i++)
  {
    for (unsigned int j=0; j<N; j++)
      if (i==j)
	A.setij(i, j, (double)(3*N));
      else A.setij(i,j, (double)(1+i+j));
    b.seti(i, 15.0); // right-hand side
    x.seti(i, 1.0); // true solution
  }
  std::cout << "A=" << std::endl;
  A.print();
  std::cout << "b=" << std::endl;
  b.print();
  std::cout << "x=" << std::endl;
  x.print();
  std::cout << "\n" << std::endl;

  double tol_input;
  std::cout << "Please enter a tolerance value: ";
  std::cin >> tol_input;
  std::cout << "\n";

  int max_iteration_input;
  std::cout << "Please enter the maximum number of iterations: ";
  std::cin >> max_iteration_input;
  std::cout << "\n";

  Vector est(N);
  unsigned int iterations;

  //boolean values check if maximum number of iterations have been reached
  bool max_reachedG = false;
  bool max_reachedJ = false;
  est.zero();

  Gauss(A, b, est, iterations, tol_input, max_iteration_input, max_reachedG);

  //If maximum number of iterations hasn't been reached print results
  if(max_reachedG == false){
    std::cout << "Gaussien took " << iterations << " iterations using overloading." << std::endl;
    std::cout << "Estimate is : " << std::endl;
    est.print();
    std::cout << "\n";
  }



  est.zero();
  Jacobi(A, b, est, iterations, tol_input, max_iteration_input, max_reachedJ);

  //If maximum number of iterations hasn't been reached print results
  if(max_reachedJ == false){
    std::cout << "Jacobi took " << iterations << " iterations." <<std::endl;
    std::cout << "Estimate is : " << std::endl;
    est.print();
  }

  return (0);
}

//////////////////
//  Set v=A*u
void MatVec(Matrix &A, Vector &u, Vector &v)
{
  unsigned int N;
  N = A.size();

  if ( (N != u.size()) || ( N != v.size() ) )
    std::cerr << "dimension mismatch in MatVec " << std::endl;
  else
  {
    for (unsigned int i=0; i<N; i++)
    {
      double x=0;
      for (unsigned int j=0; j<N; j++)
	x += A.getij(i,j)*u.geti(j);
      v.seti(i,x);
    }
  }
}

//  alpha*a + beta*b where a,b are Vectors; alpha, beta are scalars
void VecAdd (Vector &c, Vector &a, Vector &b, double alpha, double beta)
{
  unsigned int N;
  N = a.size();

  if ( (N != b.size()) )
    std::cerr << "dimension mismatch in VecAdd " << std::endl;
  else
  {
    for (unsigned int i=0; i<N; i++)
      c.seti(i, alpha*a.geti(i)+beta*b.geti(i) );
  }
}

// Use Gaus' method to solve Ax=b,
// On entry : x is the initial guess
// On exit  : x is the estimate for the solution
void Gauss(Matrix &A, Vector &b, Vector &x,
	    unsigned int &count, double tol, int max_iteration, bool &max_reachedG)
{
  unsigned int N=A.size();
  count=0;
  if ( (N != b.size()) || (N != x.size() ) )
    std::cout << "Jacobi: error - A must be the same size as b,x"
	      << std::endl;

  Matrix Dinv(N), T(N);   // The diagonal and off-diagonal matrices
  for (unsigned int i=0; i<N; i++)
    for (unsigned int j=0; j<N; j++)
      if (j != i)
      {
	T.setij(i,j, -A.getij(i,j));
	Dinv.setij(i,j,0.0);
      }
      else
      {
	T.setij(i,j, 0.0);
	Dinv.setij(i,j, 1.0/A.getij(i,j));
      }

  Matrix L(N), U(N);
  //Gets value of L
    for (unsigned int i=0; i<N; i++){
        for (unsigned int j=0; j<N; j++){
            if ( i >= j){
                L.setij(i,j, A.getij(i,j));
                U.setij(i,j, 0);
            }
            else{
                U.setij(i,j, -A.getij(i,j));
                L.setij(i,j, 0);
            }
        }
    }
    Vector r(N);

    x.zero();

    do
    {

      count++;

      //Overloading being used to calculate values
      x = (b + U*x)/L;
      r = b-A*x; // set r=b-A*r

      if (count > max_iteration){
        std::cout << "Maximum number of iterations reached before Gaussian could finish." << std::endl;
        max_reachedG = true;
        break;
    }
  }   while ( r.norm() > tol);

}

// Use Jacobi's method to solve Ax=b,
// On entry : x is the initial guess
// On exit  : x is the estimate for the solution
void Jacobi(Matrix &A, Vector &b, Vector &x,
	    unsigned int &count, double tol, int max_iteration, bool &max_reachedJ)
{
  unsigned int N=A.size();
  count=0;
  if ( (N != b.size()) || (N != x.size() ) )
    std::cout << "Jacobi: error - A must be the same size as b,x"
	      << std::endl;

  Matrix Dinv(N), T(N);   // The diagonal and off-diagonal matrices
  for (unsigned int i=0; i<N; i++)
    for (unsigned int j=0; j<N; j++)
      if (j != i)
      {
	T.setij(i,j, -A.getij(i,j));
	Dinv.setij(i,j,0.0);
      }
      else
      {
	T.setij(i,j, 0.0);
	Dinv.setij(i,j, 1.0/A.getij(i,j));
      }

  // Now implement the algorithm:
  Vector d(N), r(N);
  do
  {
    count++;
    MatVec(T,x,d);      // Set d=T*x
    VecAdd(d, b, d);    // set d=b+d (so d=b+T*x)
    MatVec(Dinv, d, x); // set x = inverse(D)*(b+T*x)

    MatVec(A, x, r);    // set r=A*x
    VecAdd(r, b, r, 1.0, -1.0); // set r=b-A*r

    if (count > max_iteration){
        std::cout << "Maximum number of iterations reached before Jacobi could finish." << std::endl;
        max_reachedJ = true;
        break;

    }

  }   while ( r.norm() > tol);

}
