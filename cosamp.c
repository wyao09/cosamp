/* Matlab Implementation here:

   % Cosamp algorithm
   %   Input
   %       K : sparsity of Sest
   %       Phi : measurement matrix
   %       u: measured vector
   %       tol : tolerance for approximation between successive solutions. 
   %   Output
   %       Sest: Solution found by the algorithm
   %
   % Algorithm as described in "CoSaMP: Iterative signal recovery from 
   % incomplete and inaccurate samples" by Deanna Needell and Joel Tropp.
   % 

   % This script creates a random signal S, and samples the signal with Phi.
   % Later it uses Cosamp to recover the random signal.

   S = zeros(1000,1);
   K = 10;
   S(1:K) = rand(K,1);
   Phi = rand(100,1000);
   Phi(Phi >= 0.5) = 1;
   Phi(Phi < 0.5) = -1;
   Phi = Phi/sqrt(100);
   u = Phi*S;

   tol = 0.01;
   maxiterations = 10;

   Sest = zeros(size(Phi,2),1);
   v = u;
   t = 1; 
   numericalprecision = 1e-12;
   T = [];

   while (t <= maxiterations) && (norm(v)/norm(u) > tol)
   %% make a guess on nonzero locations
   y = abs(Phi'*v);
  
   %% pick top 2K components
   [vals,z] = sort(y,'descend');
   Omega = find(y >= vals(2*K) & y > numericalprecision);
  
   %% merge the selected indices with those from the previous iteration
   T = union(Omega,T);

   %% reduce the system sizem and solve least squares
   b = pinv(Phi(:,T))*u;
   
   %% pick top K components
   [vals,z] = sort(abs(b),'descend');
   Kgoodindices = (abs(b) >= vals(K) & abs(b) > numericalprecision);  
   T = T(Kgoodindices);
   Sest = zeros(size(Phi,2),1);
   b = b(Kgoodindices);
   Sest(T) = b;
  
   %% compute residuals
   v = u - Phi(:,T)*b;
  
   t = t+1;
   end

*/ 

#include "blas.h"
#include "f2c.h"
#include <stdio.h>
#include "clapack.h"

#define MAX_ITER 10
int i;

main(int argc, char **argv){
  /* given */
  doublereal *Phi; //measurement matrix
  integer vector_size;
  doublereal *u; //measured vector
  doublereal tol = 0.01; //tolerance for approx between successive solutions. 
  /* end given */

  //copy u to v
  doublereal *v = malloc(vector_size*sizeof(u));
  for(i=0;i<vector_size;i++){
    v[i] = u[i];
  }
  
  int t = 0;
  integer incx = 1; // increment (usually 1)

  while((t < MAX_ITER) && 
	(dnrm2_(&vector_size, v, &incx)/dnrm2_(&vector_size, u, &incx) > tol)){

    // make a guess on nonzero locations
    y = abs(Phi* *v);
    //use complex? cgemv.f
    
      t++;
    }
  
}     
