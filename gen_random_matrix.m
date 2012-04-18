%%% THIS FILE GENERATES RANDOM MATRIX AND STORE EACH COLUMN TO A DIFFERENT
%%% FILE, THE FILES ARE LATER PASSED TO A C PROGRAM TO SYNTHESIZE THE OFDM
%%% SYMBOLS

clear;
m = 150;
n = 600;
k = 15;

%% y = A x
A = rand(m,n);
A(A<0.5) = -1;
A(A>=0.5) = 1;

x = zeros(n,1);
ind = randperm(n);
x(ind(1:k)) = rand(k,1)*10;

y = A*x;

filename = './A.dat';
A_vec = reshape(A,[],1);
write_complex_binary(filename, A_vec.');

v = read_complex_binary (filename);
A_read = reshape(v,m,n);

filename = './y.dat';
write_complex_binary(filename, y.');

y_read = read_complex_binary (filename);

filename = './x.dat';
write_complex_binary(filename, x.');

x_read = read_complex_binary (filename);

