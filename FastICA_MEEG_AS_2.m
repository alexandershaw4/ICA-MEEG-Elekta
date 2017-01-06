function S = FastICA_MEEG_AS_2(ID,NC,UL,fname,bonf)
% ICA for SPM MEEG
% - Temporal & spatial ica
% - ID is an spm meeg file
% - NC is maximum number of components
% - UL is upper limit on number to remove
% - fname is the prepend to the filename, eg. 'ica_'
% AS2016

addpath(genpath('/home/as08/old_spm12/'));
addpath('/home/as08/Downloads/FastICA_25/'); review = 1; ;
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
if nargin < 4
    fname = 'nica_';
end

% bonferroni or not
if nargin < 5
     thrp = .05;
else thrp = .05 / NC;
end


nc    = size(D,3);
for t = 1:nc
    
    fprintf('Finding components in trial %d\n',t);
    cD = squeeze(D(:,:,t)); % Chan x Samps
    
    if nargin < 2 || isempty(NC)
         [C, A, W] = fastica(cD);
    else [C, A, W] = fastica(cD,'numOfIC',NC);
    end
    
    % A = Chans x copmonents
    % C = Components x samples
    
    
    % Temporally correlated components
    %----------------------------------
    for i = 1:size(C,1)
        
        e = EOG(:,:,t);
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
    
    C(forkill,:)  = 0;
    iW(:,forkill) = 0;
    
    
    % review?
    if review && ~isempty(forkill)
        check_covar_ica(cD,(C'*iW')');
        drawnow;
    end
    
    % store
    if ~isempty(forkill)
        D(:,:,t) = ( C'*iW')';
    end
    
end

fprintf('Removed an average of %d components per trial',mean(AllGone));

% return
S = clone(ID,[fname ID.fname]);     % clone input
S(sort(unique([EEG MEG])),:,:) = D;  % update selected channel subspace
S.save;

end