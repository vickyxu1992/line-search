% function to compute alpha that min. the function f(x+alpha*d)  
%  
% The algorithm is as given in Practical Optimization by R. Fletcher,
% page 33 through 39.
% To use this code, one needs to provide the subroutines fun() and grad()
% for evaluating function value and gradient of f.
% 
% Input:
%
% x          -  current point     
% f          -  function value of f at x   
% g          -  gradient of f at x   
% d          -  line search direction
% alphamax   -  upper bound on alpha
%
%
% Output:   
%
% retcode    -  number indicating whether step length is found or not 
% ax         -  new point 
% af         -  function value of f at ax  
% ag         -  gradient of f at ax
% alpha      -  the step length found by this code
%

function [retcode, ax, af, ag, alpha] = linesearch(x, f, g, d, alpha0);

% define parameters
  
sigma   = 0.1;               %% Wolfe Powell parameter
rho     = 0.00001;            %% Goldstein parameter
tau1    = 9.0;               %% preset factor for jump of alpha
tau2    = 0.1;               %% preset factor for alpha in sectioning
tau3    = 0.5;               %% preset factor for alpha in sectioning
lowerbd = -1.0e20;           %% lower bound of objective function

[retcode, ax, af, ag, alpha] = bracketing(x, f, g, d, alpha0, sigma, rho, tau1, tau2, tau3, lowerbd);
