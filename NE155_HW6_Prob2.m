%NE 155 - Introduction to Numerical Methods in Radiation Transport
%Homework 6
%Problem 2

clc, clf, clear

%Problem Parameters
a = 4;
D = 1;
siga = 0.2;
h = 0.1;
S0 = 8;

L = sqrt(D/siga);

%Numerical Solutions of Neutron Diffusion Equation
[x1,phi_zerosrc] = DiffusionEqnSolver(a,D,siga,@(x) 0,h);
[x2,phi_constsrc] = DiffusionEqnSolver(a,D,siga,@(x) S0,h);
[x3,phi_cossrc] = DiffusionEqnSolver(a,D,siga,@(x) cos(x),h);

%Analytical Expressions for Comparison

phi_zerosrc_analytic = zeros(1,length(x1));
phi_constsrc_analytic = ((-S0/siga)/(exp(-a/L)+exp(a/L)))...
    .*(exp(-x2./L)+exp(x2./L)) + S0/siga;

C1 = (-1/(exp(-a/L)+exp(a/L)))*(cos(a)/(D*(1+1/(L^2))));

phi_cossrc_analytic = C1.*(exp(-x3./L)+exp(x3./L)) + cos(x3)./(D*(1+1/L^2));

%Plot Comparison of Numerical and Analytic Solution to NDE for Constant
%Source

figure(1)
hold on
plot(x1,phi_zerosrc,'k-');
plot(x1,phi_zerosrc_analytic,'bo');
grid on
legend('Numerical Solution','Analytical Solution')
xlabel('x (cm)')
ylabel('Flux')
title('No Source with h = 0.1');
hold off

figure(2)
hold on
plot(x1,phi_constsrc,'k-');
plot(x1,phi_constsrc_analytic,'bo');
grid on
legend('Numerical Solution','Analytical Solution')
xlabel('x (cm)')
ylabel('Flux')
title('Constant (S = 8) Source with h = 0.1');
hold off

figure(3)
hold on
plot(x3,phi_cossrc,'k-')
plot(x3,phi_cossrc_analytic,'bo')
grid on
legend('Numerical Solution','Analytical Solution')
xlabel('x (cm)')
ylabel('Flux')
title('Cosine Source with h = 0.1');
hold off
