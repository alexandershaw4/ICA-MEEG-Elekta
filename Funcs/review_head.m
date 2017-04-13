function review_head(S)

%scatter3(S.REF(:,1),S.REF(:,2),S.REF(:,3)); hold on;
scatter3(S.EEG(:,1),S.EEG(:,2),S.EEG(:,3),'b','filled'); hold on
scatter3(S.MEG(:,1),S.MEG(:,2),S.MEG(:,3),'r','filled'); 

set(gca,'visible','off');