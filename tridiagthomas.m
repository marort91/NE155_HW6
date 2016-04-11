function x = tridiagthomas(A,f)

% tridiagthomas - Tridiagonal Matrix Solver using Thomas' Algorithm
%
% Solves Ax = f for x where A is a tridiagonal matrix.
%
% Inputs
% A - N by N matrix
% f - N by 1 vector
%
% Outputs
% x - N by 1 vector

a = diag(A);
b = diag(A,-1);
c = diag(A,1);

N = length(f);
v = zeros(N,1);
y = v;

w = a(1);
y(1) = f(1)/w;

for i = 2:N
    
    v(i-1) = c(i-1)/w;
    w = a(i) - b(i-1)*v(i-1);
    y(i) = (f(i) - b(i-1)*y(i-1))/w;
    
end

for j = N-1:-1:1
    
    y(j) = y(j) - v(j)*y(j+1);
    
end

x = y;

return