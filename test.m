% This script creates a random signal S, and samples the signal with Phi.
% Later it uses Cosamp to recover the random signal.
% 
% K = 10;
% S = zeros(1000,1);
% S(1:K) = rand(K,1);
% Phi = rand(100,1000);
% Phi(Phi >= 0.5) = 1;
% Phi(Phi < 0.5) = -1;
% Phi = Phi/sqrt(100);

K=4;

Phi = read_complex_binary('A.dat');
Phi = reshape(Phi,15,50);

S = read_complex_binary('x.dat');

u = Phi*S;

tol = 0.01;
maxiterations = 10;

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

Sest = zeros(size(Phi,2),1);
v = u;
t = 1; 
numericalprecision = 1e-12;
T = [];
Ti = [];

 while (t <= maxiterations) && (norm(v)/norm(u) > tol)
  %% make a guess on nonzero locations
  y = abs(Phi'*v);
  
  %% pick top 2K components
  [vals,z] = sort(y,'descend');
  Omega = find(y >= vals(2*K) & y > numericalprecision);
  
  %% merge the selected indices with those from the previous iteration
  T = union(Omega,Ti);

  %% reduce the system sizem and solve least squares
  b = pinv(Phi(:,T))*u;
   
   %% pick top K components
   [vals,z] = sort(abs(b),'descend');
   Kgoodindices = (abs(b) >= vals(K) & abs(b) > numericalprecision);  
   Ti = T(Kgoodindices);
   Sest = zeros(size(Phi,2),1);
   b = b(Kgoodindices);
   Sest(Ti) = b;
   
   %% compute residuals
   v = u - Phi(:,Ti)*b;
  
   t = t+1;
 end
