/*
 * Cosamp algorithm
 *  Input
 *      K : sparsity of Sest
 *      Phi : measurement matrix
 *      u: measured vector
 *       tol : tolerance for approximation between successive solutions. 
 *  Output
 *      Sest: Solution found by the algorithm
 *
 * Algorithm as described in "CoSaMP: Iterative signal recovery from 
 * incomplete and inaccurate samples" by Deanna Needell and Joel Tropp.
 */ 

#include "f2c.h"
#include "stdio.h"
#include "clapack.h"

#define SIZE 4

void MAIN_(){}
void MAIN__(){}
void _MAIN_(){}

main( )
{
  char JOBU;
  char JOBVT;

  int i;
    
  integer M = SIZE;
  integer N = SIZE;
    
  integer LDA = M;
  integer LDU = M;
  integer LDVT = N;
  integer LWORK;
  integer INFO;
   
  integer mn = min( M, N );
    
  integer MN = max( M, N );
     
  double a[SIZE*SIZE] = { 16.0, 5.0, 9.0 , 4.0, 2.0, 11.0, 7.0 , 14.0, 3.0, 10.0, 6.0, 15.0, 13.0, 8.0, 12.0, 1.0};
  double s[SIZE];
  double wk[201];
  double uu[SIZE*SIZE];
  double vt[SIZE*SIZE];
     
  JOBU = 'A';
     
  JOBVT = 'A';
     
  LWORK =  201;
    
  /* Subroutine int dgesvd_(char *jobu, char *jobvt, integer *m, integer *n, 
     doublereal *a, integer *lda, doublereal *s, doublereal *u, integer *
     ldu, doublereal *vt, integer *ldvt, doublereal *work, integer *lwork, 
     integer *info)
  */

  dgesvd_( &JOBU, &JOBVT, &M, &N, a, &LDA, s, uu, 
	   &LDU, vt, &LDVT, wk, &LWORK, &INFO);
          
  printf("\n INFO=%d", INFO );          

  for ( i= 0; i< SIZE; i++ ) {
    printf("\n s[ %d ] = %f", i, s[ i ] );
  }

  return 0;
}     
