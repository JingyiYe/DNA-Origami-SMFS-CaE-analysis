function LA_hires_passive_split()
%% This function splits up data files out of hi res data that has been processed (*_bounds.mat)

% saves it as a "pdata" structure so can use with scripts for minis or
% hi-res

% Lisa, 160218

[dfile,dpath] = uigetfile('/Users/Lisa/Desktop/Analysis/062615/passive/*_bounds.mat', 'MultiSelect','on');
if ~iscell(dfile)
    dfile={dfile};
end


for nFile = 1:numel(dfile)

    disp(['Splitting ' dfile{nFile}]);
    load([dpath dfile{nFile}]);
    
    t=[];
    f=[];
    
    
    for i = 1:length(data.indices)
        t = data.time(data.indices(i, 1):data.indices(i, 2));
        f = data.force(data.indices(i, 1):data.indices(i, 2));
        
        [labels] = strcat(dfile{nFile}(1:end-11), '_', num2str(i), '_split.mat');
        
        pdata.Y_force = f;
        pdata.Time = t;
        pdata.labels = labels;
        
        save([dpath dfile{nFile}(1:end-11) '_' num2str(i) '_split.mat'], 'pdata');
    end
end