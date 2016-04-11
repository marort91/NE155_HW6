%NE 155 - Introduction to Numerical Methods in Radiation Transport
%Homework 6
%Problem 3

clc, clf, clear

%Problem Parameters
a = 4;
D = 1;
siga = 0.2;
S0 = 8;
L = sqrt(D/siga);

%Analytical Expression for Comparison (Using Constant Source)

h = [1 0.5 0.25 0.1 0.05 0.01];

for i = 1:length(h)
    
    [x1,phi,N(i)] = DiffusionEqnSolver(a,D,siga,@(x) S0,h(i));
    phi_constsrc_analytic = ((-S0/siga)/(exp(-a/L)+exp(a/L)))...
        .*(exp(-x1./L)+exp(x1./L)) + S0/siga;
    
    max_rel_err(i) = max((phi_constsrc_analytic - phi')./phi');
    
end

loglog(N,max_rel_err)
xlabel('Number of Grid Points')
ylabel('Maximum Relative Error between Calculated and True Solution')
title('Verification of Error Decrease')
grid on

m = (log(max_rel_err(end)) - log(max_rel_err(1)))/(log(h(end))-log(h(1)));

s = sprintf('Slope = %f \n',m);
legend(s)
