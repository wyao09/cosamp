#include <stdlib.h>
#include "blaswrap.h"
#include "f2c.h"
#include <stdio.h>
#include "clapack.h"

#define MAX_ITER 1
int i,j,l;

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

void print_tuple(tuple *t, int n){
  printf("\n");
  for(i=0;i<n;i++){
    printf("%f %d\n", t[i].value, t[i].index);
  }
}

//need to optimize this
int cmp (const void *a, const void *b){
  tuple pa = *(const tuple*) a;
  tuple pb = *(const tuple*) b;
  pa.value = pa.value - pb.value;
  if(pa.value == 0)
    return 0;
  if(pa.value < 0)
    return 1;
  return 0;
}

/* SET Starts */

//make bit array later/use better set implementation
void reset(int *set, int n){
  for(i=0;i<n;i++)
    set[i] = 0;
}

int size(int *set, int n){
  l = 0;
  for(i=0;i<n;i++)
    if (set[i])
      l++;
  return l;
}

/* SET Ends */

main(int argc, char **argv){
  /* given */
  int k = 2;//sparsity
  integer m = 6;
  integer n = 16;
  int mn = max(m,n);//max of m,n
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
		  1,1,1,1,1,0,
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

  //copy u to v; size of v is max(m,n) since it will store b later on
  doublereal *v = (doublereal*) malloc(mn*sizeof(doublereal));
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
  tuple b_tuple[mn];

  integer info;
  integer nrhs = 1; // number of columns of matrices B and X in B = AX
  integer la = m;
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

    // sort y
    qsort(y, n, sizeof(tuple), cmp);

    // merge top 2k with T (this can be a lot better)
    // may need to change to indicies
    for(i=0;i<2*k;i++){
      T[y[i].index] = 1;
      printf("%d ",y[i].index);
    }
    printf("\n");
    
    //reduce Phi
    int l = 0;
    for(i=0;i<n;i++){
      for(j=0;j<m;j++){
	if(T[i]){
	  b[l] = Phi[i*m+j];
	  l++;
	}
      }
    }

    print_matrix(v,m,1);

    // reduce system size and solve for least square
    // ? = Phi y
    // answer stored in v, which is a problem since we want to multiply by u each time
    trans = 'N';
    integer ni = (integer) size(T,n);
    integer lb = (integer) max((int)ni,(int)m);
    integer lwork = min(m,ni) + max(min(m,ni),nrhs);
    doublereal *work = (doublereal*)  malloc( lwork*sizeof(doublereal));
    printf("nrhs:%d,m:%d,n:%dlda:%d,ldb:%d,lwork:%d \n\n",
	   (int)nrhs,(int)m,(int)ni,(int)la,(int)lb,(int)lwork);

    //this is broken
    dgels_(&trans,&m,&ni,&nrhs,b,&la,v,&lb,work,&lwork,&info);

    free(work);

    printf("\n info: %d\n\n", (int)info);

    print_matrix(v,ni,1);

    //move abs(v) to b_tuple
    for(i=0;i<mn;i++){
      if(v[i] < 0)
	b_tuple[i].value = -1*v[i];
      else
	b_tuple[i].value = v[i];
      b_tuple[i].index = i;
    }

    // include numbericalprecision?
 
    // abo(b) then sort T, swith between 2k and 3k or mn??
    if(t==0)
      qsort(b_tuple, 2*k, sizeof(doublereal), cmp);
    else
      qsort(b_tuple, 3*k, sizeof(doublereal), cmp);

    //print_tuple(b_tuple, 2*k);
    reset(T,n);

    for(i=0;i<k;i++){
      printf("%d\n",i);
      T[b_tuple[i].index] = 1;
    }

    /*

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
