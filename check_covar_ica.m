function check_covar_ica(D,C)


Cv1 = covar(D);
Cv2 = covar(C);


title('Covariance among channels')
subplot(131),imagesc(Cv1); title('Pre ICA')
subplot(132),imagesc(Cv2); title('Post ICA');
subplot(133),imagesc(Cv1-Cv2); title('Pre - Post');