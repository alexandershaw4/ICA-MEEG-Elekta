function S = FastICA_MEEG_AS(ID,NC)

addpath('~/Downloads/FastICA_25/'); review = 1;

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


% Max no. components 
if nargin < 2;
    NC = size(D,1);
end


nc    = size(D,3);
thrp  = .5 / NC;
for t = 1:nc
    
    fprintf('Finding components in trial %d\n',t);
    cD = squeeze(D(:,:,t)); % Chan x Samps
    
    
    [C, A, W] = fastica(cD,'numOfIC',NC);
    
    % A = Chans x copmonents
    % C = Components x samples
    
    
    % Temporally correlated components
    %----------------------------------
    for i = 1:size(C,1)
        
        e = EOG(:,:,t);
        c = C  (i,:);
        
        [~,p]=corr(c', e');
        
        if any(p < thrp)
            fprintf('found EOG correlated component: %d ... ',i);
            
            while any(p < thrp)
                fprintf('robust removing\n');
                C(i,:) = Orthog(c,e,p,thrp);
                c      = C(i,:);
                [~,p]  = corr(c', e');
            end
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
        [~,p] = corr(iWP(:,i),PL);
        
        if any(p < thrp)
            fprintf('topography correlated GRAD component: %d ... ',i);
            
            while any(p < thrp)
                fprintf('robust removing\n');
                iWP(:,i) = Orthog(iWP(:,i)',PL',p,thrp);
                [~,p]    = corr(iWP(:,i),PL);
            end
        end
        
        
        % Magnetometer topogrpahies
        [~,p] = corr(iWM(:,i),MG);
        
        if any(p < thrp)
            fprintf('topography correlated MAG component: %d ... ',i);
            
            while any(p < thrp)
                fprintf('robust removing\n');
                iWM(:,i) = Orthog(iWM(:,i)',MG',p,thrp);
                [~,p]    = corr(iWM(:,i),MG);
            end
        end        
    end
    
    
    
    % review?
    if review
        check_covar_ica(cD,(C'*iW')');
        drawnow;
    end
    
    % store
    D(:,:,t) = ( C'*iW')';
    
end


% return
S = clone(ID,['ica_' ID.fname]);     % clone input
S(sort(unique([EEG MEG])),:,:) = D;  % update selected channel subspace
S.save;

end

function O = Orthog(c,e,p,thrp)

v = find(p<thrp); % vert or horz

for j = 1:length(v)
    % Regression apparatus
    
    try   O = O + (.5*c)+(.5*(c-polyval(polyfit(e(v(j),:),c,1),e(v(j),:))));
    catch O =     (.5*c)+(.5*(c-polyval(polyfit(e(v(j),:),c,1),e(v(j),:))));
    end
end

O = O/j;

end

