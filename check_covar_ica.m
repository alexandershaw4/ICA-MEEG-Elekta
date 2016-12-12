function check_covar_ica(D,C)

while prod(size(D)) > 1e6
    D = HighResMeanFilt(D,.2,1);
    C = HighResMeanFilt(C,.2,1);
end

Cv1 = covar(D);
Cv2 = covar(C);



title('Covariance among channels')
subplot(131),imagesc(Cv1); title('Pre ICA')
subplot(132),imagesc(Cv2); title('Post ICA');
subplot(133),imagesc(Cv1-Cv2); title('Pre - Post');

