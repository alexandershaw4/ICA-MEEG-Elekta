function check_covar_ica(D,C)

while prod(size(D)) > 1e4
    D = HighResMeanFilt(D,.2,1);
    C = HighResMeanFilt(C,.2,1);
    fprintf('rescaling...\n');
end

Cv1 = covar(full(D));
Cv2 = covar(full(C));


title('Covariance among channels [may be rescaled]')
subplot(131),imagesc(Cv1); title('Pre ICA')
subplot(132),imagesc(Cv2); title('Post ICA');
subplot(133),imagesc(Cv1-Cv2); title('Pre - Post');

