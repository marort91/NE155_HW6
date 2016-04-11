function [x, phi, N] = DiffusionEqnSolver(a,D,siga,src,h)

% DiffusionEqnSolver - Neutron Diffusion Equation Solver
%
% Solves one-dimensional slab neutron diffusion equation with source using 
%the finite difference method.
%
% Inputs
% a - width of slab
% D - diffusion constant
% siga - macroscopic absorption cross section
% src - source function (function handle)
% h - discretization width
%
% Outputs
% x - location of calculated flux values
% phi - neutron flux

% Calculation of calculation points
N = 2*a/h;
x = -a+h:h:a-h;

% Generation of Finite Difference Matrix to be Solved
M = (-D/h^2).*full(gallery('tridiag',N-1,1,-2,1)) + siga*eye(N-1,N-1);

% Generation of Source Term

S = zeros(N-1,1);

for i = 1:N-1

    S(i) = src(x(i));
    
end

phi = tridiagthomas(M,S);

return