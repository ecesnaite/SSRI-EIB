% (C) Elena Cesnaite, email: e.cesnaite@gmail.com, page: https://www.researchgate.net/profile/Elena-Cesnaite

% This code was created to analyze data described in a paper:'One-week escitalopram intake shifts excitation-inhibition balance 
% in healthy female brain: implications for individual cortical responsivity to SSRIs' by Zsido & Molloy et al. 
% The code was adapted from a version developed by Dr Keyvan Mahjoory.
% The code applies cluster analysis according to
% http://www.sciencedirect.com/science/article/pii/S0165027007001707


%OUTPUTS:
%   CL   Output of permutation
%   CLNoPerm    unpermuted data

clear all, close all
saveDir = '';

load(''); % channel labels
load('')  % 2D channel locations
chanNeis = load('');%channel neighbourhood information
nRnd = 1000; % number of permutations

% match the order of channel neighbour and channel label
[indx pl] = ismember(chanLabels,{chanNeis.neighbours.label});

chanNeis.neighbours = chanNeis.neighbours(pl);

load('') % load EEG derivatives
parameter_name = 'slope';% slope, broadpow, iaf UPDATE
parameter = slope; % slope, broadpow, iaf UPDATE

% group 
groups = readtable(''); % group information: drug or placebo

day = 1; %1 - baseline, 2 - D1 or 3- D7 
ind_day = find(groups.Day==day)%
ind_group = groups(ind_day,:)

nsub = size(ind_group,1);
nchan = 27;

%% Generating random indices:

indxRnd = [];  % Dimension of the matrix: nsubjects * nRnd
for i = 1:nRnd
    indxRnd = [indxRnd randperm(nsub)'];
end

%% Unpermuted clusters
fprintf(' Detecting clusters on non-permuted data \n')
data_day = parameter(ind_day,:);%get parameters of that day alone

for ichan = 1:nchan
    [p,h,stats] = ranksum(data_day(strcmp(ind_group.Group,'SSRI'),ichan), data_day(strcmp(ind_group.Group,'Placebo'),ichan)); %Wilcoxon rank sum test
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
cfg.chanLocs2D = locs_2D_new; 
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
        
        pvals = zeros(nchan,1); %just to make computation faster
        zvals = zeros(nchan,1); %same
        
        perm_group = randperm(nsub);
        ind_rand = [ind_group(:,1) ind_group(perm_group,2)];
        
        for ichan=1:nchan
            [p,h,stats] = ranksum(data_day(strcmp(ind_rand.Group,'SSRI'),ichan), data_day(strcmp(ind_rand.Group,'Placebo'),ichan)); %Wilcoxon rank sum test
            zvals(ichan,1) = stats.zval;
            pvals(ichan,1) = p;
        end
        
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

plotTitle = ['Day ', num2str(day), ' SSRI vs Placebo'];

if ~isempty(CL) % No significant cluster
    CLCell = struct2cell(CL');
    pvalClust = length(find( abs(cell2mat(CLCell(2,:))) > abs(clNoPerm.tmax) ))/nRnd;
    plotTitle = [plotTitle,'  P= ' num2str(pvalClust)]; % significance of a cluster
end
