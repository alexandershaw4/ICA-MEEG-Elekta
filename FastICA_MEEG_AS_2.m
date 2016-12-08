function S = FastICA_MEEG_AS_2(ID,NC)

addpath('~/Downloads/FastICA_25/'); review = 1; w = 2;

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
        
        [~,p] = corr(c', e');
        
        if any(p < thrp)
            fprintf('found EOG correlated component: %d ... \n',i);
            C(i,:) = c*0;

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
            fprintf('topography correlated GRAD component: %d ... \n',i);
            iWP(:,i) = iWP(:,i)*0;
        end
        
        
        % Magnetometer topogrpahies
        [~,p] = corr(iWM(:,i),MG);
        
        if any(p < thrp)
            fprintf('topography correlated MAG component: %d ... \n',i);
            iWM(:,i) = iWM(:,i)*0;
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
S = clone(ID,['nica_' ID.fname]);     % clone input
S(sort(unique([EEG MEG])),:,:) = D;  % update selected channel subspace
S.save;

end

% function O = Orthog(c,e,p,thrp,W)
% 
% v = find(p<thrp); % vert or horz
% 
% 
% if ~isempty(v)
%     for j = 1:length(v)
%         % Regression apparatus
%         
%         s       = HighResMeanFilt(polyval(polyfit(e(v(j),:),c,1),e(v(j),:)),1,8);  
%         try   O = O + s;
%         catch O =     s;
%         end
%     end
%     
%     O = O/j;
% else
%     O = c;
% end
% 
% end


