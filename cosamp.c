#include <stdlib.h>
#include "blaswrap.h"
#include "f2c.h"
#include <stdio.h>
#include "clapack.h"

#define MAX_ITER 10
int i,j;

void print_matrix(doublereal *p, int m, int n){
  for(i=0;i<m;i++){
    for(j=0;j<n;j++){
      printf("%f ",(double)p[j+i*j]);
    }
    printf("\n");
  }
  printf("\n");
}


int cmp (const void *a, const void *b){
  doublereal pa = *(const doublereal*) a;
  doublereal pb = *(const doublereal*) b;
  pa = pa - pb;
  if(pa == 0)
    return 0;
  if(pa > 0)
    return 1;
  return 0;
}

main(int argc, char **argv){
  /* given */
  int k = 2;//sparsity
  integer m = 6;
  integer n = 16;
  doublereal Phi[m*n]; //measurement matrix
  doublereal u[m]; //measured vector
  doublereal tol = 0.01; //tolerance for approx between successive solns. 
  /* end given */

  // populate Phi and u
  int phi_sample[96] = {1,1,1,0,1,1,
		  0,0,1,0,0,1,
		  1,1,1,1,0,0,
		  0,0,1,0,0,0,
		  1,1,1,1,1,1,
		  0,1,1,0,0,1,
		  1,0,1,1,0,0,
		  0,0,0,0,0,1,
		  0,0,0,0,0,0,
		  1,0,1,0,0,1,
		  0,0,1,0,1,1,
		  1,1,1,1,1,1,
		  0,1,1,0,1,0,
		  1,1,1,1,1,0,
		  1,0,1,1,1,0,
		  1,0,1,0,1,1};

  for(i=0;i<64;i++)
    Phi[i] = phi_sample[i];

  int y_sample[6] = {0,0,28,0,5,5};

  for(i=0;i<6;i++)
    u[i] = y_sample[i];

  doublereal *T = malloc(3*n*sizeof(doublereal));// this should probably be 3k

  //copy u to v
  doublereal *v = malloc(m*sizeof(doublereal));
  for(i=0;i<m;i++){
    v[i] = u[i];
  }
  
  int t = 0;
  integer incx = 1; // increment (usually 1)

  char trans = 'C';
  doublereal alpha = 1;
  doublereal beta = 0;
  integer lda = m;//1;
  // does this actually need to be malloced? look at abs
  doublereal *y = malloc(n*sizeof(doublereal));

  integer info;
  integer nrhs = 1; // number of columns of matrices B and X in B = AX
  integer lwork = max(1,min(m,n) + max(min(m,n),nrhs));
  doublereal *work = malloc( lwork*sizeof(doublereal));
  integer la = m;
  integer lb = max(m,n);
  doublereal *x = malloc(k*sizeof(doublereal));
 
  // COSAMP Starts Here
  while((t < MAX_ITER) && 
	(dnrm2_(&m, v, &incx)/
	 dnrm2_(&m, u, &incx) > tol)){
    printf("%d\n",t);

    // Phi* *v
    trans = 'C';
    dgemv_(&trans,&m,&n,&alpha,Phi,&lda,v,&incx,&beta,y,&incx);
    
    print_matrix(y,n,1);

    // y = abs(Phi* *v)
    // may have to loop with dcabs1_(doublecomplex z)
    for (i=0;i<n;i++){
      if(y[i] < 0)
	y[i] = -1 * y[i];
    }



    // sort y
    qsort(y, n, sizeof(doublereal), cmp);
    
    // merge top 2k with T (this can be a lot better)
    // may need to change to indicies
    for(i=k;i<3*k;i++){
      T[i] = y[i];
    }
    qsort(T, 3*k, sizeof(doublereal), cmp);

    
    // reduce system size and solve for least square
    // ? = Phi y
    // answer stored in T
    trans = 'N';
    dgels_(&trans,&m,&n,&nrhs,Phi,&la,T,&lb,work,&lwork,&info);
 
    // abo(T) then sort T
    for (i=0;i<3*k;i++){
      if(T[i] < 0)
	T[i] = -1 * T[i];
    }
    qsort(T, 3*k, sizeof(doublereal), cmp);

    //A*x_1
    trans = 'N';
    dgemv_(&trans,&m,&n,&alpha,v,&lda,Phi,&incx,&beta,x,&incx);

    //populate v / compute residual
    for (i=0;i<k;i++){
      v[i] = u[i] - x[i];
    }

    //print_matrix(x,16,1);

    t++;
  }
}     
