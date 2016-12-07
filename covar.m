function C = covar(X,t)

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