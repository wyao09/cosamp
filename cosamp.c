#include <stdlib.h>
#include <string.h>
#include "blaswrap.h"
#include "f2c.h"
#include <stdio.h>
#include "clapack.h"
#include <sys/time.h>
#include <sys/resource.h>

#define MAX_ITER 10

/* Todo: - get complex version to work
 *       - fix load
 *       - organize arguments
 */

/* Globals */
int i,j,l;

//timer
double timer_start1[10];
double timer_end1[10];
double timer_start2[10];
double timer_end2[10];
double timer_start3[10];
double timer_end3[10];
double timer_start4[10];
double timer_end4[10];
double timer_start5[10];
double timer_end5[10];
double timer_start6[10];
double timer_end6[10];
double timer_start7[10];
double timer_end7[10];
double timer_start8[10];
double timer_end8[10];
double timer_start9[10];
double timer_end9[10];


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
double get_time();


typedef struct forest_node_simple_t {
    void* value;
    struct forest_node_simple_t* parent;
} forest_node_simple;

forest_node_simple* FindSimple(forest_node_simple* node) {
    while (node->parent != NULL) {
        node = node->parent;
    }
    return node;
}

int SizeSimple(forest_node_simple* node){
  int n = 1;
   while (node->parent != NULL) {
     node = node->parent;
     n++;
   }
   return n;
}

void UnionSimple(forest_node_simple* node1, forest_node_simple* node2) {
    node2->parent = node1; /* or node1->parent = node2; */
    // need to loop through and check for duplicate
}

forest_node_simple* MakeSetSimple(void* value) {
    forest_node_simple* node = malloc(sizeof(forest_node_simple));
    node->value = value;
    node->parent = NULL;
    return node;
}


/*
  y = Phi * x
  given y and Phi, we want to solve for x
  3 15 40
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
  int iter = 0;
  int k = atoi(argv[1]);//sparsity
  integer m = atoi(argv[2]);
  integer n = atoi(argv[3]);
  doublereal Phi[m*n]; //measurement matrix
  doublereal y[m]; //measured vector

  /* auxiliary */
  int mn = max(m,n);
  int T[n]; //stores indicies
  int Ti[n]; //stores indicies of indicies
  doublereal v[n]; //working copy of y
  doublereal w[mn]; //working copy of y;  replaced during least square
  doublereal Phi_reduced1[m*n];
  doublereal Phi_reduced2[m*n];
  doublereal b_reduced[m*n];//get rid of this one?
  doublereal guess[n];
  doublereal x[k];
  tuple guess_t[n];
  tuple b_tuple[mn];
  doublereal work[m+n];

  integer one = 1;
  char trans = 'C';
  doublereal alpha = 1;
  doublereal beta = 0;
  integer info;
  double val;
  integer ni;
  integer lba;
  integer lwork;

  /* constants */
  doublereal tol = 0.01; //tolerance for approx between successive solns. 

  // populate Phi and y from file and reset indicies
  loadComplexMatrixFromFile(Phi, m, n, Phi_file);
  printf("hit here\n");
  loadComplexMatrixFromFile(y, m, 1, y_file);
  reset(T,n);
  reset(Ti,n);

  // copy y to v
  for(i=0;i<m;i++){
    v[i] = y[i];
  }
 
  /* start timer */
  double start = get_time();

  // COSAMP Starts Here
  while((iter < MAX_ITER) && (dnrm2_(&m,v,&one)/dnrm2_(&m,y,&one) > tol)){
    timer_start1[iter] = get_time();
    trans = 'C';
    dgemv_(&trans,&m,&n,&alpha,Phi,&m,v,&one,&beta,guess,&one);

    // y = abs(Phi* *v) and add index
    // may have to loop with dcabs1_(doublecomplex z)
    for (i=0;i<n;i++){
      if(guess[i] < 0)
	guess_t[i].value = -1 * guess[i];
      else
	guess_t[i].value = guess[i];
      guess_t[i].index = i;
    }
    timer_end1[iter] = get_time();
    timer_start2[iter] = get_time();
    // sort guess_t
    qsort(guess_t, n, sizeof(tuple), cmp);
    timer_end2[iter] = get_time();
    timer_start3[iter]=get_time();
    // merge top 2k with T (this can be a lot better)
    // may need to change to indicies
    val = guess_t[2*k-1].value;
    for(i=0;i<n;i++){
      if(guess_t[i].value >= val){
	T[guess_t[i].index] = 1;
      }
      else
	break;
    }
    timer_end3[iter]=get_time();
    timer_start4[iter]=get_time();
    //reduce Phi
    l = 0;
    for(i=0;i<n;i++){
      for(j=0;j<m;j++){
	if(T[i]){
	  Phi_reduced1[l] = Phi[i*m+j];
	  Phi_reduced2[l] = Phi[i*m+j];
	  l++;
	}
      }
    }
    timer_end4[iter]=get_time();
    timer_start5[iter]=get_time();
    // reduce system size and solve for least square, store answer in w
    trans = 'N';
    ni = size(T,n);
    lba = max(ni,m);
    lwork = min(m,ni) + max(min(m,ni),one);

    //make copy of y since w will be replaced
    for(i=0;i<m;i++){
      w[i] = y[i];
    }
    timer_end5[iter]=get_time();
    timer_start6[iter]=get_time();
    dgels_(&trans,&m,&ni,&one,Phi_reduced1,&m,w,&lba,work,&lwork,&info);

    if(info!=0){
      printf("matrix has illegal value or does not have full rank\n");
      return 1;
    }

    timer_end6[iter]=get_time();
    timer_start7[iter]=get_time();
    //move abs(w) to b_tuple
    for(i=0;i<mn;i++){
      if(w[i] < 0)
	b_tuple[i].value = -1*w[i];
      else
	b_tuple[i].value = w[i];
      b_tuple[i].index = i;
    }
    
    // include numbericalprecision?
 
    // abo(b) then sort T
    qsort(b_tuple, ni, sizeof(tuple), cmp);
    timer_end7[iter]=get_time();
    timer_start8[iter]=get_time();
    //create a new Ti for storing indicies of indicies
    reset(Ti,n);

    val = b_tuple[k-1].value;

    for(i=0;i<mn;i++){
      if(b_tuple[i].value >= val)
	Ti[b_tuple[i].index] = 1;
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

    //reduce Phi (could be done using T)
    l= 0;
    for(i=0;i<ni;i++){
      if(Ti[i]){
	for(j=0;j<m;j++){
	  Phi_reduced1[l] = Phi_reduced2[i*m+j];
	  l++;
	}
      }
    }
    timer_end8[iter]=get_time();
    timer_start9[iter]=get_time();
    trans = 'N';
    ni = size(T,n);
    dgemv_(&trans,&m,&ni,&alpha,Phi_reduced1,&m,
	   b_reduced,&one,&beta,guess,&one);

    //populate v / compute residual
    for (i=0;i<m;i++){
      v[i] = y[i] - guess[i];
    }
    timer_end9[iter]=get_time();
    iter++;
  }
  /* end timer */
  double end = get_time();

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
  printf("\nBenchmark: %f sec, %d iterations\n",end-start,iter);
  for(i=0;i<10;i++){
    printf("1: %f\n",timer_end1[i]-timer_start1[i]);
  }
  for(i=0;i<10;i++){
    printf("2: %f\n",timer_end2[i]-timer_start2[i]);
  }
  for(i=0;i<10;i++){
    printf("3: %f\n",timer_end3[i]-timer_start3[i]);
  }
  for(i=0;i<10;i++){
    printf("4: %f\n",timer_end4[i]-timer_start4[i]);
  }
  for(i=0;i<10;i++){
    printf("5: %f\n",timer_end5[i]-timer_start5[i]);
  }
  for(i=0;i<10;i++){
    printf("6: %f\n",timer_end6[i]-timer_start6[i]);
  }
  for(i=0;i<10;i++){
    printf("7: %f\n",timer_end7[i]-timer_start7[i]);
  }
  for(i=0;i<10;i++){
    printf("8: %f\n",timer_end8[i]-timer_start8[i]);
  }
  for(i=0;i<10;i++){
    printf("9: %f\n",timer_end9[i]-timer_start9[i]);
  }
}

/* Function Implementations */

void loadComplexMatrixFromFile(double *buffer, int m, int n, char *filename){
  int filesize = m*n;
  //segfaults here if filesize is too big
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
  //return pb.value - pa.value;
  pa.value = pa.value - pb.value;
  if(pa.value == 0)
    return 0;
  if(pa.value < 0)
    return 1;
  return 0;
}

double get_time(){
  struct timeval t;
  struct timezone tzp;
  gettimeofday(&t, &tzp);
  return t.tv_sec + t.tv_usec*1e-6;
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
