
% make a fake signal
global noisemodel
noisemodel.x = randn(1,2000);

o =     make_oscillations(50,1,4);
o = [o;  make_oscillations(20,1,4)];
o = [o;  make_oscillations(80,1,4)];
o = [o;  make_oscillations(40,1,4)];



if any(size(o)==1)
    
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

else
    % ica
    [C,A,W]=fastica(o);
    
    P  = ( C'*pinv(W)' )';
    ns = size(o,1);
    
    for i = 1:ns
        subplot(ns,1,i),...
            plot(o(i,:),'b'); hold on;
            dW = pinv(W);
            plot(mean( (C(i,:)'*dW(:,i)')'),'r');

    end
end

