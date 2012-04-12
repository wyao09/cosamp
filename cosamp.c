#include <stdlib.h>
#include <string.h>
#include "blaswrap.h"
#include "f2c.h"
#include <stdio.h>
#include "clapack.h"

#define MAX_ITER 10

/* Todo: - get complex version to work
 *       - fix load
 *       - organize arguments
 */

/* Globals */
int i,j,l;

typedef struct{
  doublereal value;
  int index;
} tuple;

/* Functions */
void loadComplexMatrixFromFile(double *buffer, int m, int n, char *filename);
void print_matrix(doublereal *p, int m, int n);
void print_t(int *T, int n);
void print_tuple(tuple *t, int n);
int cmp (const void *a, const void *b);
void reset(int *set, int n);
int size(int *set, int n);

/*
  y = Phi * x
  given y and Phi, we want to solve for x
 */

main(int argc, char **argv){
  if(argc !=4 && argc !=6){
    printf("%d\n",argc);
    printf("usage: cosamp [sparsity] [m] [n] [Phi filename] [y filename]\n");
    return 1;
  }
  char Phi_file[64];
  char y_file[64];
  if(argc == 6){
    //file io
    strcpy(Phi_file, argv[4]);
    strcpy(y_file,argv[5]);
  }

  /* given */
  int k = atoi(argv[1]);//sparsity
  integer m = atoi(argv[2]);
  integer n = atoi(argv[3]);
  doublereal Phi[m*n]; //measurement matrix
  doublereal y[m]; //measured vector

  /* auxiliary */
  int mn = max(m,n);
  int T[(int)n]; //stores indicies
  int Ti[(int)n]; //stores indicies of indicies
 
  /* constants */
  doublereal tol = 0.01; //tolerance for approx between successive solns. 
  
  // populate Phi and y from file and reset indicies
  loadComplexMatrixFromFile(Phi, m, n, Phi_file);
  loadComplexMatrixFromFile(y, m, 1, y_file);
  reset(T,n);
  reset(Ti,n);

  //copy y to v and w; size of v is max(m,n) since it will store b later on
  doublereal *v = (doublereal*) malloc(mn*sizeof(doublereal));
  doublereal *w = (doublereal*) malloc(mn*sizeof(doublereal));
  for(i=0;i<m;i++){
    v[i] = y[i];
  }

  doublereal b[m*n];
  doublereal Phi_reduced[m*n];
  doublereal b_reduced[m*n];

  int t = 0;
  integer incx = 1; // increment (usually 1)

  char trans = 'C';
  doublereal alpha = 1;
  doublereal beta = 0;
  integer lda = m;//1;
  // need to keep original indicies
  tuple *y_t = (tuple*) malloc(n*sizeof(tuple));
  doublereal y_i[n];
  tuple b_tuple[mn];

  integer info;
  integer nrhs = 1; // number of columns of matrices B and X in B = AX
  integer la = m;
  doublereal *x = (doublereal*) malloc(k*sizeof(doublereal));
 
  // COSAMP Starts Here
  while((t < MAX_ITER) && 
	(dnrm2_(&m, v, &incx)/
	 dnrm2_(&m, y, &incx) > tol)){
    printf("%d\n",t);

    // Phi* *v
    trans = 'C';
    dgemv_(&trans,&m,&n,&alpha,Phi,&lda,v,&incx,&beta,y_i,&incx);

    // y = abs(Phi* *v) and add index
    // may have to loop with dcabs1_(doublecomplex z)
    for (i=0;i<n;i++){
      if(y_i[i] < 0)
	y_t[i].value = -1 * y_i[i];
      else
	y_t[i].value = y_i[i];
      y_t[i].index = i;
    }

    // sort y_t
    qsort(y_t, n, sizeof(tuple), cmp);

    //
    double val = y_t[2*k-1].value;

    // merge top 2k with T (this can be a lot better)
    // may need to change to indicies
    for(i=0;i<n;i++){
      if(y_t[i].value >= val){
	T[y_t[i].index] = 1;
      }
      else
	break;
    }
    
    //reduce Phi
    int l = 0;
    for(i=0;i<n;i++){
      for(j=0;j<m;j++){
	if(T[i]){
	  b[l] = Phi[i*m+j];
	  Phi_reduced[l] = Phi[i*m+j];
	  l++;
	}
      }
    }

    // reduce system size and solve for least square
    // ? = Phi y
    // answer stored in v, which is a problem since we want to multiply by u each time
    trans = 'N';
    integer ni = (integer) size(T,n);
    integer lb = (integer) max((int)ni,(int)m);
    integer lwork = min(m,ni) + max(min(m,ni),nrhs);
    doublereal *work = (doublereal*)  malloc( lwork*sizeof(doublereal));

    //make copy of y
    for(i=0;i<m;i++){
      w[i] = y[i];
    }

    dgels_(&trans,&m,&ni,&nrhs,b,&la,w,&lb,work,&lwork,&info);

    free(work);

    if((int)info!=0){
      printf("matrix has illegal value or does not have full rank\n");
      return 1;
    }

    //move abs(w) to b_tuple
    for(i=0;i<mn;i++){
      if(w[i] < 0)
	b_tuple[i].value = -1*w[i];
      else
	b_tuple[i].value = w[i];
      b_tuple[i].index = i;
    }

    // include numbericalprecision?
 
    // abo(b) then sort T, swith between 2k and 3k? probably best to find size every time
    qsort(b_tuple, size(T,n), sizeof(tuple), cmp);
    
    //instead of resetting, need to create a new Ti for this and use original Phi
    reset(Ti,n);

    val = b_tuple[k-1].value;

    // this needs to be fixed
    for(i=0;i<mn;i++){
      if(b_tuple[i].value >= val){
	Ti[b_tuple[i].index] = 1;
      }
      else
	break;
    }
    l = 0;
    for(i=0;i<n;i++){
      if(T[i]){
	if(!Ti[l])
	  T[i] = 0;
	l++;
      }
    }
    
    //reduce b (aka w)
    l=0;
    for(i=0;i<n;i++){
      if(Ti[i]){
	b_reduced[l] = w[i];
	l++;
      }
    }

    doublereal Phi_reduced2[sizeof(doublereal)*m*size(T,n)];

    //reduce Phi (could be done using T)
    l= 0;
    for(i=0;i<ni;i++){ //or is it n?
      if(Ti[i]){
	for(j=0;j<m;j++){
	  Phi_reduced2[l] = Phi_reduced[i*m+j];
	  l++;
	}
      }
    }

    trans = 'N';
    ni = size(T,n);
    dgemv_(&trans,&m,&ni,&alpha,Phi_reduced2,&lda,b_reduced,&incx,&beta,y_i,&incx);

    //populate v / compute residual
    for (i=0;i<m;i++){
      v[i] = y[i] - y_i[i];
    }

    t++;
  }

  printf("recovered x:\n");
  l=0;
  for(i=0;i<n;i++){
    if(T[i]){
      while(!Ti[l])
	l++;
      printf(" %f\n",w[l]);
      l++;
    }
    else
      printf(" 0\n");      
  }
}

/* Function Implementations */

void loadComplexMatrixFromFile(double *buffer, int m, int n, char *filename){
  int filesize = m*n;
  double tmp[filesize*2];
  FILE *f;
  int r;

  f = fopen(filename, "rb");
  if (f){
    r = fread(tmp, filesize*sizeof(double)*2, 1, f);
  }
  // Error opening file
  else{
    printf("ERROR: read failed\n");
    exit(2);
  }
  // Copy real portion to buffer
  for(r=0;r<filesize;r++){
    buffer[r] = tmp[r*2];
  }
}

// Prints matrix in row major order (stored in column major)
void print_matrix(doublereal *p, int m, int n){
  for(i=0;i<m;i++){
    for(j=0;j<n;j++){
      printf("%f ",(double)p[i+j*m]);
    }
    printf("\n");
  }
  printf("\n");
}

void print_t(int *T, int n){
  for(i=0;i<n;i++){
    if(T[i])
      printf("T:%d\n",i+1);
  }
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

/* SET Functions */

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
