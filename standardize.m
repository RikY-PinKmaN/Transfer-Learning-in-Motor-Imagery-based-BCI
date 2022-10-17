function [D,sm,sd]=standardize(d)
% STANDARDIZE Standardizes data
%   by Ang Kai Keng (kkang@pmail.ntu.edu.sg)
%
%   Syntax: 
%     [D,sm,sd]=standardize(d)
%   where
%     D  : standardized data
%     sm : standard mean
%     sd : standard deviation
%     d  : data set matrix with size [sN,sD] where
%          sN          = number of data tuples,
%          sD-1        = number of input dimensions,
%          lastcol     = class id [0..nN-1].
%
%   See also STANDARDIZEW.

% Determine the size of original data
[sN sD]=size(d);

% Compute mean
sm=mean(d);

% Compute standard deviation
sd=std(d);

% Subtract mean
td=d-repmat(sm,sN,1);

% Divide by standard deviation
D=td./repmat(sd,sN,1);