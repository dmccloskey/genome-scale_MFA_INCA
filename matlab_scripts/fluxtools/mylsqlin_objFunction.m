function [out] = mylsqlin_objFunction(x,Prob)
% Objective function for mylsqlin
%   Solves the function 
%   0.5*(NORM(C*X-D)).^2

C = Prob.user.C;
D = Prob.user.D;

out = 0.5*(norm(C*x-D)).^2;

