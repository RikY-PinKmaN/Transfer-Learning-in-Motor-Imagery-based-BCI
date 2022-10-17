% RCSP Nonlinear Constraint
%   by Mahnaz Arvaneh(stuma@i2r.a-star.edu.sg)
function [c, ceq] = confun(x)
global sigmac
global sigmat
c=[];
ceq=[x*sigmat*x'-1];