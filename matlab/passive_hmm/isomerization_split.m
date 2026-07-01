function isomerization_split()
%% this function takes in a "_split.mat" file, queries for an isomerization
% event, and then split the data into trans and cis files

% input box will ask for number of isomerization events- 0 will go to next
% file. It does allow you to select two or more regions as isomerized.

% remember that trans state is normal hopping, cis is impaired.


[dfile,dpath] = uigetfile('/Users/Lisa/Desktop/Analysis/*_split.mat', 'MultiSelect','on');
if ~iscell(dfile)
    dfile={dfile};
end

for nFile = 1:numel(dfile)

    disp(['Splitting ' dfile{nFile}]);
    load([dpath dfile{nFile}]);
    
    temp_pdata=pdata; % for each file, use the pdata as a temp_pdata 
    %variable because want to save pdata in final file.
    
    t = FilterAndDecimate(temp_pdata.Time, 10);
    f = FilterAndDecimate(temp_pdata.Y_force, 10);
    
    scrsz = get(0,'ScreenSize');
    g = figure('Position',[1 scrsz(4) scrsz(3) scrsz(4)]); %make it full screen!
    plot(t, f, 'b')
    set(gca, 'FontSize', 16, 'FontWeight', 'bold')
    xlabel('Time (s)')
    ylabel('Force (pN)')
    title('Passive data at 250 Hz')
    
    regions = inputdlg('How many isomerization events are there? Click beginning and end of each.');
    regions = str2num(regions{1});
    
    previous= 1; % start at index 1, will update each time
    
    if regions > 0
        time_index = zeros(regions, 2);
        index = zeros(regions, 2);
        temp_index=[]; % temporary index for each region        
        for j = 1: regions
            figure(g)
            [x, ~] = ginput(2);
            time_index(j, 1) = x(1);
            time_index(j, 2) = x(2);
            temp_index = find(temp_pdata.Time > time_index(j, 1) & temp_pdata.Time< time_index(j, 2));
            index(j, 1) = temp_index(1); %beginning of isomerization
            index(j, 2) = temp_index(end); % end of isomerization
            
            % each cycle save one cis and one trans event
            pdata.Time = temp_pdata.Time(previous:index(j, 1)-1);
            pdata.Y_force = temp_pdata.Y_force(previous:index(j, 1)-1);
            save([dpath dfile{nFile}(1:end-4) '_trans', num2str(j), '.mat'],'pdata');
            
            pdata.Time = temp_pdata.Time(index(j, 1):index(j, 2));
            pdata.Y_force = temp_pdata.Y_force(index(j, 1):index(j, 2));
            save([dpath dfile{nFile}(1:end-4) '_cis', num2str(j), '.mat'],'pdata');
            
            previous = index(j, 2)+1;
        end
        pdata.Time = temp_pdata.Time(previous:end);
        pdata.Y_force = temp_pdata.Y_force(previous:end);
        save([dpath dfile{nFile}(1:end-4) '_trans', num2str(regions+1), '.mat'],'pdata');
    end
end
end

        