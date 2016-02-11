function [out] = mylsqlin_gradFunction(x,Prob)
% Gradient of the objective function for mylsqlin
%   Solves the derviative of the function 
%   0.5*(NORM(C*X-D)).^2
%   = 0.5*(ABS(C*X-D)).^2

C = Prob.user.C;
D = Prob.user.D;

out = 0.5*(abs(C*x-D)).^2;

