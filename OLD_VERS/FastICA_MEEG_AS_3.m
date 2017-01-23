function S = FastICA_MEEG_AS_3(ID,NC,UL,fname,time,bonf)
% Window based ICA for SPM MEEG
% - Temporal & spatial ica
% - ID is an spm meeg file
% - NC is maximum number of components
% - UL is upper limit on number to remove
% - fname is the prepend to the filename, eg. 'ica_'
% - this version longer epochs [spec time in s]:
%   concatenates n trials / epochs to make time window, runs ica & correls,
%   adjusts data, then places back into epoch / trial positions in data
% AS2016

addpath(genpath('/home/as08/old_spm12/'));
addpath('/home/as08/Downloads/FastICA_25/'); review = 0; plotcor=1;
if ischar(ID); ID = spm_eeg_load(ID); end

EOG(1) = find(strcmp(ID.chanlabels,'EOG061'));
EOG(2) = find(strcmp(ID.chanlabels,'EOG062'));
EOG    = ID(EOG,:,:);


% Find MEG & EEG channels
MEG = strfind(ID.chanlabels,'MEG');
MEG = find(~cellfun(@isempty,MEG));

EEG = strfind(ID.chanlabels,'EEG');
EEG = find(~cellfun(@isempty,EEG));

% Neurophys channels only
D = ID(sort(unique([EEG MEG])),:,:);


% Channel locations
S.EEG = ID.sensors('EEG').chanpos;
S.MEG = ID.sensors('MEG').chanpos;
S.REF = [S.EEG; S.MEG];

% Channel labels & topographies
LAB = ID.chantype(MEG);
TOP = load('MEGArtifactTemplateTopographies');
PLN = strfind(LAB,'MEGPLANAR');
PLN = find(~cellfun(@isempty,PLN));
MAG = find(~ismember(MEG,PLN));


% % Max no. components 
% if nargin < 2 
%     NC = size(D,1);
% end

% Max number rejected components
if nargin < 3
    UL = [];
end

% new file name [clone]
if nargin < 4 || isempty(fname)
    fname = 'nica_';
end

% time search
if nargin < 5
    time = 10 ;
end

% bonferroni or not
if nargin < 6
     thrp = .05;
else thrp = .05 / NC;
end

% sort num trials to concat
t     = ID.time;
tocat = round(time/t(end))+1;

nc    = round(size(ID,3)/tocat)-1;
win   = 1;
for t = 1:nc
    clear p_grad p_temp p_mag
    
    fprintf('Finding components in time window %d\n',t);
    
    cD  = squeeze(D(:,:,[win:win+tocat-1])); % Chan x Samps x [ntrials in time]
    cD  = reshape(cD,[size(cD,1) size(cD,2)*size(cD,3)]);
    e   = EOG(:,:,[win:win+tocat-1]);
    e   = reshape(e,[size(e,1) size(e,2)*size(e,3)]);
    
    fprintf('including trials %d to %d in this window\n',win,win+tocat);
    win = win + tocat;
    
    if nargin < 2 || isempty(NC)
         [C, A, W] = fastica(cD);
    else [C, A, W] = fastica(cD,'numOfIC',NC);
    end
    
    % A = Chans x copmonents
    % C = Components x samples
    
    
    % Temporally correlated components
    %----------------------------------
    for i = 1:size(C,1)
        
        %e = EOG(:,:,t);
        c = C  (i,:);
        
        [Q,p] = corr(c', e');
        
        if any(p < thrp)
            fprintf('found EOG correlated component: %d ... \n',i);
            p_temp(i,:) = p;
        end
    end
    

    
    % Topography correlated
    %----------------------------------
    iW  = pinv(W);
    iWP = iW(PLN,:);
    iWM = iW(MAG,:);
    
    PL = [TOP.VEOG{2} TOP.HEOG{2}];
    MG = [TOP.VEOG{1} TOP.HEOG{1}];
    
    for i = 1:size(iW,2);       
        
        % Planar topographies
        [Q,p] = corr(iWP(:,i),PL);
        
        if any(p < thrp)
            fprintf('topography correlated GRAD component: %d ... \n',i);
            p_grad(i,:) = p;
        end
        
        
        % Magnetometer topogrpahies
        [Q,p] = corr(iWM(:,i),MG);
        
        if any(p < thrp)
            fprintf('topography correlated MAG component: %d ... \n',i);
            p_mag(i,:) = p;
        end        
    end
    
    
    % sort out common spatial and temporal component correlates
    %-----------------------------------------------------------
    try p_grad ; catch p_grad = []; end
    try p_mag  ; catch p_mag  = []; end
    try p_temp ; catch p_temp = []; end
    
    [i1 i2]    = find(p_grad);
    [i3 i4]    = find(p_mag);
    [tmpc,pos] = find(p_temp);
    topos      = unique([i1(:);i3(:)]);
        
    forkill = tmpc(ismember(tmpc,topos));
    forkill = unique(forkill);
    if ~isempty(UL)
        try forkill = forkill(1:UL); end
    end
    
    fprintf('removing %d components\n',length(forkill));
    
    try AllGone(t,:) = length(forkill); end
    
    if plotcor && ~isempty(forkill)
        subplot(411), plot(e(1,:));hold on;plot(e(2,:));title('EOGs');
        subplot(412), plot(cD'),title('Original');
        subplot(413), plot( (C(forkill,:)'*iW(:,forkill)' ) ); title('Correlated components');
        C(forkill,:)  = 0;
        iW(:,forkill) = 0;        
        subplot(414), plot( (C'*iW') ); title('With Removal');
        drawnow
    else
        C(forkill,:)  = 0;
        iW(:,forkill) = 0;
    end
    
    
    
    % review?
    if review && ~isempty(forkill)
        check_covar_ica(cD,(C'*iW')');
        drawnow;
    end
    

    
    % store
    if ~isempty(forkill)
        %D(:,:,t) = ( C'*iW')';
        new =reshape( ( C'*iW')' ,[size(cD,1) size(ID,2) (tocat)]);
        D(:,:,[win:win+tocat-1]) = new;
    end
    
end

fprintf('Removed an average of %d components per trial\n',mean(AllGone));

% return
S = clone(ID,[fname ID.fname]);     % clone input
S(sort(unique([EEG MEG])),:,:) = D;  % update selected channel subspace
S.save;

end