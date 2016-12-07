

D = spm_eeg_load('Mmafespmeeg_Rov123_C1_S1_sss'); DD = D;
I = spm_eeg_load('ica_Mmafespmeeg_Rov123_C1_S1_sss'); II = I;

% find MEEG channels
MEG = strfind(D.chanlabels,'MEG');
MEG = find(~cellfun(@isempty,MEG));

EEG = strfind(D.chanlabels,'EEG');
EEG = find(~cellfun(@isempty,EEG));

% Neurophys channels only
D = D(sort(unique([EEG MEG])),:,:);
I = I(sort(unique([EEG MEG])),:,:);

xD = covar(D);
xI = covar(I);

cD = covar(D,1);
cI = covar(I,1);

figure; fs = 20;

[m,M] = minmax(spm_vec([xD xI]));

subplot(221),imagesc(xD);caxis([m M]); title('per sample covariance, no ica','fontsize',fs);
subplot(222),imagesc(xI);caxis([m M]); title('per sample covariance, with ica','fontsize',fs);

[m,M] = minmax(spm_vec([cD cI]));

subplot(223),imagesc(cD);caxis([m M]); title('channel covariance, no ica','fontsize',fs);
subplot(224),imagesc(cI);caxis([m M]); title('channel covariance, with ica','fontsize',fs);


% decompose covariance to channel by first 3 eigen components
figure,subplot(121),imagesc(PEig(cD,1:3));subplot(122),imagesc(PEig(cI,1:3));

% decompose covariance to samples by first 3 eigen components
figure,subplot(121),imagesc(PEig(xD,1:3));subplot(122),imagesc(PEig(xI,1:3));

% MAP for t{1}
DJ=DD.inv{1}.inverse.J{1};
IJ=II.inv{1}.inverse.J{1};

figure,
%subplot(121),imagesc(DJ); subplot(122),imagesc(IJ);
subplot(121),plot(DJ); subplot(122),plot(IJ);

cDJ = covar(full(DJ),1);
cIJ = covar(full(IJ),1);