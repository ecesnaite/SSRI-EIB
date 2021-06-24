% (C) Elena Cesnaite, email: e.cesnaite@gmail.com, page: https://www.researchgate.net/profile/Elena-Cesnaite

% This code was created to analyze data described in a paper:'One-week escitalopram intake shifts excitation-inhibition balance 
% in healthy female brain: implications for individual cortical responsivity to SSRIs' by Zsido & Molloy et al. 
% The code was adapted from a version developed by Dr Keyvan Mahjoory.
% applies cluster analysis according to
% http://www.sciencedirect.com/science/article/pii/S0165027007001707

%OUTPUTS:
%   CL   Output of permutation
%   CLNoPerm    unpermuted data

clear all, close all

load('') % load EEG derivatives
parameter = slope; % slope, broadpow, iaf UPDATE
parameter_name = 'slope';

day1 = 'baseline'
day2 = 'steady state' % 'steady state' baseline' 'single dose'

group = 'Placebo' % 'SSRI' 'Placebo'

saveDir = '';

load(''); % channel labels
load('')  % 2D channel locations
chanNeis = load('');%channel neighbourhood information
nRnd = 1000; % number of permutations

% match the order of channel neighbour and channel label
[indx pl] = ismember(chanLabels,{chanNeis.neighbours.label});

chanNeis.neighbours = chanNeis.neighbours(pl);

% group
groups = readtable('');% group information: drug or placebo

% group into SSRI and Placebo
if strcmp(group, 'SSRI')
    ind_ssri = ismember(groups.Group,'SSRI')
    single_group = groups(ind_ssri,:);
    
    data_single_group = parameter(ind_ssri,:);%get parameters for one group (e.g. ssri)
    
elseif strcmp(group, 'Placebo')
    ind_plac = ismember(groups.Group,'Placebo')
    single_group = groups(ind_plac,:)
    
    data_single_group = parameter(ind_plac,:);%get parameters for one group (e.g. placebo)
else
    error('wrong group input')
end

ind_baseline = find(single_group.Day==1)
ind_single = find(single_group.Day==2)
ind_steady = find(single_group.Day==3)

% take day info
if strcmp(day1, 'baseline')
    data_day1 = data_single_group(ind_baseline,:); % get parameters for baseline alone
elseif strcmp(day1, 'single dose')
    data_day1 = data_single_group(ind_single,:); % get parameters
elseif strcmp(day1, 'steady state')
    data_day1 = data_single_group(ind_steady,:); % get parameters
end

% take day2 info
if strcmp(day2, 'baseline')
    data_day2 = data_single_group(ind_baseline,:); % get parameters for baseline alone
    
elseif strcmp(day2, 'single dose')
    data_day2 = data_single_group(ind_single,:); % get parameters
    
elseif strcmp(day2, 'steady state')
    data_day2 = data_single_group(ind_steady,:); % get parameters
end

nsub = size(ind_baseline,1);
nchan = 27;

%% Generating random indices:

indxRnd = [];  % Dimenstion of the matrix: nsubjects * nRnd
for i = 1:nRnd
    indxRnd = [indxRnd randperm(nsub)'];
end

%% Unpermuted clusters
fprintf(' Detecting clusters on non-permuted data \n')

for ichan = 1:nchan
    [p,h,stats] = signrank(data_day2(:,ichan), data_day1(:,ichan),'method', 'approximate');%Wilcoxon signed rank test
    hvals0(ichan,1) = h;
    pvals0(ichan,1) = p;
    zvals0(ichan,1) = stats.zval;
end

data = [];
data.chanLabels = chanLabels;
data.chanNeighbours = chanNeis;
data.pvals = pvals0;
data.tvals = zvals0;
cfg.sigThresh= 0.05;
cfg.chanLocs2D = locs_2D_new; %chanlocs
CL0 = k1_find_clusters_sensor_SSRI(data,cfg);  %% detecting clusters on non-permuted data

if isempty(CL0.neg) && ~isempty(CL0.pos)
    CLNoPerm.tmax= CL0.pos.tvalMax;
elseif ~isempty(CL0.neg) && isempty(CL0.pos)
    CLNoPerm.tmax= CL0.neg.tvalMax;
elseif ~isempty(CL0.neg) && ~isempty(CL0.pos)
    if abs(CL0.neg.tvalMax) > CL0.pos.tvalMax
        CLNoPerm.tmax = CL0.neg.tvalMax;
    elseif abs(CL0.neg.tvalMax) <= CL0.pos.tvalMax
        CLNoPerm.tmax = CL0.pos.tvalMax;
    end
else
    error('tmax computation failed!')
end

CLNoPerm.pvals = pvals0;
CLNoPerm.tvals = zvals0;
CLNoPerm.clusters = CL0;
CLNoPerm.nsigChan = CL0.nsigChan;

if CLNoPerm.nsigChan ==0 % No need to run permutation
    CL = [];
    warning('no significant cluster has been found')
else  % Run permutation
    
    %% Permutation analysis
    fprintf(' Permutation over questionnaire data')
    fprintf('\n Iteration:     ')
    
    %% Randomization of indices
    for iRnd = 1:nRnd
        
        tvals = zeros(nchan,1);
        pvals = zeros(nchan,1);
        
        ind_perm = randperm(nsub);
        
        data_day1_permute = data_day1(ind_perm,:);
        
        for ichan=1:nchan
            
            [p,h,stats] = signrank(data_day2(:,ichan), data_day1_permute(:,ichan),'method', 'approximate');% Wilcoxon signed-rank test
            hvals(ichan,1) = h;
            pvals(ichan,1) = p;
            zvals(ichan,1) = stats.zval;
            
        end
        
        data = [];
        data.chanLabels = chanLabels;
        data.chanNeighbours = chanNeis;
        data.pvals = pvals;
        data.tvals = zvals;
        CL(iRnd).perm= k1_find_clusters_sensor_SSRI(data,cfg); 
        
        %%% COMPUTING max t-value
        if isempty(CL(iRnd).perm.neg) &&  ~isempty(CL(iRnd).perm.pos)
            CL(iRnd).tmax = CL(iRnd).perm.pos.tvalMax;
        elseif ~isempty(CL(iRnd).perm.neg) &&  isempty(CL(iRnd).perm.pos)
            CL(iRnd).tmax = CL(iRnd).perm.neg.tvalMax;
        elseif ~isempty(CL(iRnd).perm.neg) &&  ~isempty(CL(iRnd).perm.pos)
            if abs(CL(iRnd).perm.neg.tvalMax) > CL(iRnd).perm.pos.tvalMax
                CL(iRnd).tmax = CL(iRnd).perm.neg.tvalMax;
            elseif abs(CL(iRnd).perm.neg.tvalMax) <= CL(iRnd).perm.pos.tvalMax
                CL(iRnd).tmax = CL(iRnd).perm.pos.tvalMax;
            end
        else
            error('t-max computation failed!');
        end
        
        
        fprintf([repmat('\b',1,length(num2str(iRnd))), '%s'],num2str(iRnd))
        
        
    end  %%  permutaion
    
    fprintf('\n');
end

%gives out: CLNoPerm CL indxRnd
clNoPerm = CLNoPerm;
indxPerm = indxRnd;

plotTitle = [parameter_name,' ', day1, ' vs ', day2, ' ', group];

if ~isempty(CL) % No significant cluster
    CLCell = struct2cell(CL');
    pvalClust = length(find( abs(cell2mat(CLCell(2,:))) > abs(clNoPerm.tmax) ))/nRnd;
    plotTitle = [plotTitle,'  P= ' num2str(pvalClust)]; % significance of a cluster
end

