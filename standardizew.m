function D=standardizew(d,sm,sd)
%STANDARDIZEW Standardizes data with supplied mean and deviation
%   by Ang Kai Keng (kkang@pmail.ntu.edu.sg)
%
%   Syntax: 
%     D=standardizew(d,sm,sd)
%   where
%     D  : standardized data
%     d  : data set matrix with size [sN,sD] where
%          sN          = number of data tuples,
%          sD-1        = number of input dimensions,
%          lastcol     = class id [0..nN-1].
%     sm : standard mean
%     sd : standard deviation
%
%   See also STANDARDIZE.

% Determine the size of original data
[sN sD]=size(d);

% Subtract mean
tD=d-repmat(sm,sN,1);

% Divide by standard deviation
D=tD./repmat(sd,sN,1);
