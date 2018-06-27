function C = covar2(X,t)
% fast but crude covariance estimate of n-dim matrix X
%
% iteratively smooths and resamples until small enough to quickly take the
% inner product
%
% if using for spm meeg dataset, try something like
% cv = covar(D(D.indchantype('MEGPLANAR'),:,:),1);
% 
%
% AS

while prod(size(X)) > 1e7
    % limit to n million points
    X = HighResMeanFilt(X,.25,4);
end

% optional X or X'
try t; catch t = 0; end


% cat to 2D
X     = VecRetainDim(X);



if t;
    X = X';
end



[n,m] = size(X);

Xm    = mean(X,2);
Xc    = X - Xm*ones(1,m);

C     = (1/m)*Xc'*Xc;