function LA_hires_passive_processing()
%% This function is a first attempt at processing hi-res passive data
% It will (hopefully): 
% 1. concatenate the data of a hires file and plot force vs time.
% 2. Allow selection of passive regions

% a processed file is saved of concatenated time and force and the
% boundaries of the regions of passive sections. It has a _bounds.mat
% extension.

% Lisa, 150630


[dfile,dpath] = uigetfile('/Users/Lisa/Desktop/Analysis/062515/passive/*.mat', 'MultiSelect','on');
if ~iscell(dfile)
    dfile={dfile};
end

for nFile = 1:numel(dfile)
    disp(['Processing ' dfile{nFile}]);
    load([dpath dfile{nFile}]);
    
    
    ncell = length(stepdata.time);
    time_concat =[];
    force_concat=[];
    ext_concat=[];
    
    for i = 1:ncell % hi res traces are divided into cells, need to concatenate.
        time_concat = [time_concat stepdata.time{i}];
        force_concat = [force_concat stepdata.force{i}];
        ext_concat = [ext_concat stepdata.extension{i}];
    end
    
    filt_t=FilterAndDecimate(time_concat, 10);
    filt_f=FilterAndDecimate(force_concat, 10);
    
    scrsz = get(0,'ScreenSize');
    g = figure('Position',[1 scrsz(4) scrsz(3) scrsz(4)]); %make it full screen!
    plot(time_concat, force_concat, 'Color',[0.6 0.6 0.6])
    hold on
    plot(filt_t, filt_f, 'b')
    set(gca, 'FontSize', 16, 'FontWeight', 'bold')
    xlabel('Time (s)')
    ylabel('Force (pN)')
    title('Passive data at 250 Hz')
    
    regions = inputdlg('How many force regions to be analyzed?');
    regions = str2num(regions{1});
    
    bounds = zeros(regions, 2); %make a variable of the boundaries of the regions of interest;
    % each row is time1 and time2; want region in
    % between those 2 times.
    indices = zeros(regions, 2); %make a variable of the indices of the regions of interest;
    % each row is index1 and index2; will avoid
    % having to run loop again in future.
    temp_index=[]; % temporary index for each region
    
    for j = 1: regions
        figure(g)
        [x, ~] = ginput(2);
        bounds(j, 1) = x(1);
        bounds(j, 2) = x(2);
        temp_index = find(time_concat > bounds(j, 1) & time_concat< bounds(j, 2));
        indices(j, 1) = temp_index(1);
        indices(j, 2)=temp_index(end);
        
    end
    
    %only saving fields we want. Save file with new extension so have original file.
    data=struct;
    data.bounds = bounds;
    data.time = time_concat;
    data.force = force_concat;
    data.extension = ext_concat;
    data.indices=indices;
    
    save([dpath dfile{nFile}(1:end-4) '_bounds.mat'],'data');

end

  
% if multiselect is on, place all code within this loop and have dfile reference nFile:

% for nFile = 1:numel(dfile)
%     disp(['Processing ' dfile{nFile}]);
%     load([dpath dfile{nFile}]); 
% end