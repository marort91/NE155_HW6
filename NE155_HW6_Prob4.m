%NE 155 - Introduction to Numerical Methods in Radiation Transport
%Homework 6
%Problem 4

clc, clf, clear

%Problem Parameters
a = 4;
D = 1;
siga = 0.7;
nusigf = 0.6;

h = 0.1;

% Calculation of calculation points
N = 2*a/h;
x = -a+h:h:a-h;

% Generation of Finite Difference Matrix to be Solved
M = (-D/h^2).*full(gallery('tridiag',N-1,1,-2,1)) + siga*eye(N-1,N-1);

% Initial Flux Guess
phi0 = ones(N-1,1);

% Tolerance Criteria
ktol = 1e-4;
phitol = 1e-10;

% Intial Eigenvalue Guess
k = 1.0;

for iter = 1:1000
    
    kprev = k;
    
    [psi, iteration] = GaussSeidelSolve(M,nusigf.*phi0,100000,phitol,'absolute');
    fprintf('Inner iterations on flux converged in: %i \n',iteration);
    
    k = sum(nusigf.*psi)/sum(nusigf*phi0);
    k_residual = abs(k-kprev);
    
    phi0 = (1/k).*psi;
    
    fprintf('k-effective residual: %f\n',k_residual);
    fprintf('Iteration: %i \n',iter);
    
    if ( k_residual <= ktol )
        
        break
        
    end
    
end

plot(x,phi0,'ko')
xlabel('x (cm)');
ylabel('Flux')
title('Criticality Eigenvalue Calculation for Slab of Width 8, h = 0.1');
