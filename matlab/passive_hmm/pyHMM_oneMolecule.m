
function LA_pyHMM_oneMolecule()
%% this function uses the python HMM library to analyze multiple data files.

% the function needs _split.mat files which have been cropped to one trap
% position. It runs the python script through a matlab shell.

% the dependencies are: LA_pyHMM2.m, LA_hmm_lifetimes_n_v2.m

% the function tests different numbers of states and different initial Tij
% values. It picks the best one (based on BIC) and keeps it.

% This makes a Bell plot, which is nice for checking visually but not
% necessary since the oneConstruct script makes an aggregate weighted Bell
% plot. So you could comment out this part for speed, but be sure to still
% leave the loop before the figure on since it checks that the mu's are in
% order and reorders if necessary.

% Lisa, April 19, 2016

%pyversion('/Users/Lisa/miniconda3/bin/python'); only need this line the
%first time run the code per session
py.importlib.import_module('hmmlearn');

minN = 3; %will test between these numbers of states.
maxN = 3; %3;

bias = [0.1 0.05 0.01 0.005 0.001 0.0005 0.0001]; % larger number -> faster transitions favored. 
% p=0.2 corresponds to average lifetime of ~5 data points, so is barely detectable

[dfile,dpath] = uigetfile('/Users/Lisa/Desktop/Analysis/*.mat', 'MultiSelect','on');
if ~iscell(dfile)
    dfile={dfile};
end

name = dfile{1}(6:14); % this should be the date of the file. Format is phage121715N02_1_split.mat 

x = numel(dfile);
mu = cell(x, 1); % will have all mu's for each file
sigma = cell(x,1); 
trajectory =cell(x,1); % will have all trajectories
k = cell(x,1); %initial rates output from function; only has k1, k2, etc (1 column)
rates = cell(x,1); %final "rate matrix" (0s on diagonal, not true rate matrix);
Nij = cell(x,1);
CI = cell(x,1); %confidence interval. 
data = cell(x,1);
error1 = zeros(x, 1); %for plotting error bars.
error2 = zeros(x, 1); %for plotting error bars

for nFile = 1:numel(dfile)
    disp(['Analyzing ' dfile{nFile}]);
    load([dpath dfile{nFile}]);
    
    data{nFile}=pdata;
    
    tempBIC = zeros(maxN-minN+1, length(bias));
    tempMu = cell(maxN-minN+1, length(bias));
    tempTraj = cell(maxN-minN+1, length(bias));
    tempSigma = cell(maxN-minN+1, length(bias));
    tempTij = cell(maxN-minN+1, length(bias));
    
    for n =1:maxN-minN+1
        for i = 1:length(bias)
            [tempMu{n, i}, tempSigma{n, i}, tempTraj{n, i}, tempBIC(n, i), tempTij{n,i}] = LA_pyHMM2(pdata.Y_force, minN+n-1, bias(i));
        end
    end
        
%     disp(tempBIC)
    [~,ind]=min(tempBIC(:)); %find the minimum BIC value
    [bestN, bestBias] = ind2sub(size(tempBIC), ind); % get the subscripts of that value
%     
%     tempTij{bestN, bestBias}
    
    %remember that "bestN" starts at 1, not minN, so need to add minN-1 to
    %value
    
    mu{nFile} = tempMu{bestN, bestBias};
    sigma{nFile} = tempSigma{bestN, bestBias};
    trajectory{nFile} = tempTraj{bestN, bestBias};
    [k{nFile}, Nij{nFile}, CI{nFile}] = LA_hmm_lifetimes_n_v2(trajectory{nFile}, bestN+minN-1);
    
    %plots the data with the fit overlaid. Again, not necessary but nice to
    %check.
    figure
    hold on
    plot(pdata.Time, pdata.Y_force, pdata.Time, mu{nFile}(round(trajectory{nFile}+1)))
    %scatter(pdata.Time, pdata.Y_force)
    plot(pdata.Time, mu{nFile}(round(trajectory{nFile}+1)), 'r')
    hold on

    % this part makes all the rates square matrices. For 2 states does not
    % change values, just adds zero diagonals. for >2 it separates k12
    % and k13 from k1obs  
    for i = 1:bestN+minN-1
        rates{nFile}(i, :) = k{nFile}(i).*Nij{nFile}(i,:)./sum(Nij{nFile}(i,:));
    end

end

%% construct a Bell plot

force1 = zeros(x,1); %will extract from the mu's. Probably not necessary per se but makes plotting script more legible.
force2 = zeros(x,1);
k12 = zeros(x,1);
k21 = zeros(x,1);
% make aggregate plot of data (just 1<->2 transitions)
for i = 1:numel(dfile)  
     % want to make sure every state 1 is lowest state so need to check order.
    if any(mu{i} ~= sort(mu{i}))
       disp('Forces are not in order, reordering now.')
       [mu{i}, ind] = sort(mu{i});
       sigma{i} = sigma{i}(ind);
       trajectory{i} = ind(trajectory{i}+1)-1; % trajectory starts at 0, not 1.
       rates{i} = rates{i}(ind, ind);
       CI{i} = CI{i}(ind, :);
       Nij{i} = Nij{i}(ind, ind);
    end
    force1(i) = mu{i}(1);
    k12(i) = rates{i}(1,2);   
    error1(i) = CI{i}(1,2)-rates{i}(1,2); % this is 2*sigma
     % the interval from the fit is symmetric, just find the magnitude of the error.

    force2(i) =  mu{i}(2);
    k21(i) =  rates{i}(2,1);
    error2(i) = CI{i}(2,2)-rates{i}(2,1);   
end 

nodata = find(k12==0); %remove points with no rate so doesn't mess up when
 % divide by k later (sometimes zero if <5 observations)
k12(nodata)=[];
force1(nodata) =[];
error1(nodata) =[];

nodata2 = find(k21==0); %remove points with no rate so doesn't mess up when
 % divide by k later (sometimes zero if <5 observations)
k21(nodata2)=[];
force2(nodata2) =[];
error2(nodata2) =[];

figure; hold on % (2) S.Y. 11.13.16
plot(force1, log(k12), '.r', 'MarkerSize',18)
hold on
plot(force2, log(k21), 'ob', 'MarkerSize',6)
 % note that the error on k is found from the fit, but we are plotting ln(k)
 % the error of ln(k) can be approximated as error of k divided by k
 % this is still 2*sigma (95% CI)
errorbar(force1, log(k12), error1./k12, '.r')
errorbar(force2, log(k21), error2./k21, '.b')
set(gca, 'FontSize', 16, 'FontWeight', 'bold')

% linear fit 1
fitvars = polyfit(force1, log(k12), 1);
m = fitvars(1);
c = fitvars(2);
linearx=0:0.1:15;
linearf=m*linearx+c;
plot(linearx,linearf, '-r')
hold on
% linear fit 2
fitvars = polyfit(force2, log(k21), 1);
m = fitvars(1);
c = fitvars(2);
linearx=0:0.1:15;
linearf=m*linearx+c;
plot(linearx,linearf, '--b')
xlabel('Force (pN)  ')
ylabel('ln(k)  ')
title(['Bell plot for molecule ', name, '  '])

results=struct;
results.data = data;
results.mu= mu;
results.sigma=sigma;
results.trajectory=trajectory;
results.rates = rates;
results.CI = CI;
results.Nij = Nij; %added this in, early files will not have this field.

save([dpath name '_pyfit.mat'], 'results');