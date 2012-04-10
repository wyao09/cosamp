#include <stdlib.h>
#include "blaswrap.h"
#include "f2c.h"
#include <stdio.h>
#include "clapack.h"

#define MAX_ITER 1
int i,j;

typedef struct{
  doublereal value;
  int index;
} tuple;

// converts from column major order
void print_matrix(doublereal *p, int m, int n){
  for(i=0;i<m;i++){
    for(j=0;j<n;j++){
      printf("%f ",(double)p[i+j*m]);
    }
    printf("\n");
  }
  printf("\n");
}

//need to optimize this
int cmp_tuple (const void *a, const void *b){
  tuple pa = *(const tuple*) a;
  tuple pb = *(const tuple*) b;
  pa.value = pa.value - pb.value;
  if(pa.value == 0)
    return 0;
  if(pa.value > 0)
    return 1;
  return 0;
}

//make bit array later/use better set implementation
void reset(int *set, int n){
  for(i=0;i<n;i++)
    set[i] = 0;
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
  double phi_sample[96] = {1,1,1,0,1,1,
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

  for(i=0;i<m*n;i++)
    Phi[i] = phi_sample[i];

  double y_sample[6] = {0,0,28,0,5,5};

  for(i=0;i<m;i++)
    u[i] = y_sample[i];

  int T[(int)n]; //stores indicies
  reset(T,n);

  //copy u to v
  doublereal *v = (doublereal*) malloc(m*sizeof(doublereal));
  for(i=0;i<m;i++){
    v[i] = u[i];
  }

  doublereal b[m*n];

  
  int t = 0;
  integer incx = 1; // increment (usually 1)

  char trans = 'C';
  doublereal alpha = 1;
  doublereal beta = 0;
  integer lda = m;//1;
  // need to keep original indicies
  tuple *y = (tuple*) malloc(n*sizeof(tuple));
  doublereal y_i[n];

  integer info;
  integer nrhs = 1; // number of columns of matrices B and X in B = AX
  integer lwork = m + n;// max(1,min(m,n) + max(min(m,n),nrhs));
  doublereal *work = (doublereal*)  malloc( lwork*sizeof(doublereal));
  integer la = m;
  integer lb = m;//max(m,n);
  doublereal *x = (doublereal*) malloc(k*sizeof(doublereal));
 
  // COSAMP Starts Here
  while((t < MAX_ITER) && 
	(dnrm2_(&m, v, &incx)/
	 dnrm2_(&m, u, &incx) > tol)){
    printf("%d\n",t);

    // Phi* *v
    trans = 'C';
    dgemv_(&trans,&m,&n,&alpha,Phi,&lda,v,&incx,&beta,y_i,&incx);

    // y = abs(Phi* *v) and add index
    // may have to loop with dcabs1_(doublecomplex z)
    for (i=0;i<n;i++){
      if(y_i[i] < 0)
	y[i].value = -1 * y_i[i];
      else
	y[i].value = y_i[i];
      y[i].index = i;
    }

    for(i=0;i<n;i++){
      printf("%f %d\n",y[i].value, y[i].index);
    }

    // sort y
    qsort(y, n, sizeof(doublereal), cmp_tuple);



    // merge top 2k with T (this can be a lot better)
    // may need to change to indicies
    for(i=0;i<2*k;i++){
      T[y[i].index] = 1;
      printf("%d ",y[i].index);
    }
    printf("\n");
    
    int k = 0;
    for(i=0;i<n;i++){
      for(j=0;j<m;j++){
	if(T[j]){
	  b[k] = Phi[i*n+j];
	  k++;
	}
      }
    }

    //print_matrix(b,m,4);

    /*

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

    */
    //print_matrix(x,16,1);

    t++;
  }
}     
