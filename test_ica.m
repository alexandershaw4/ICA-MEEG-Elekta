
% make a fake signal
global noisemodel
noisemodel.x = randn(1,2000);

o =     make_oscillations(50,2,4);
o = o + make_oscillations(20,2,4);

% ica
[C,A,W]=fastica(o'*o);

% ~ reconstruction
P = ( C'*pinv(W)' )';

P = sqrt(diag(P));

% correct neg values
i = find(o<0);
P(i) = P(i)*-1;

% plot
plot(o,'b');hold on; plot(P,'r--');