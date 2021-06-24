% (C) Elena Cesnaite, email: e.cesnaite@gmail.com, page: https://www.researchgate.net/profile/Elena-Cesnaite

% This code was created to analyze data described in a paper:'One-week escitalopram intake shifts excitation-inhibition balance 
% in healthy female brain: implications for individual cortical responsivity to SSRIs' by Zsido & Molloy et al. 
% The code allows for manual inspection of noisy datasets. Later, bad data segments are rejected to calculate ICA weights. 
% Last updated 22.06.2021

clc; clear all; close all
eeglab

% Load data %
saveDir = ''
dataDir = ''
files=dir([dataDir, '*.set']);

%Get the right dataset from the Manual_Inspection folder based on the
%iteration number in the main script
p = []
while isempty(p)
    iterationNb  = input('Which iteration would you like to load? (e.g. 1 for BDI001_bad_segments_1.set) ', 's'); 
        day  = input('Which day would you like to load? (BL, D1 or D7) ', 's'); 

    idx = endsWith({files.name},['_',day, '_manual_inspection_bad_segments_', iterationNb, '.set']);
    p = find(idx);
end

%Load Data
EEG = pop_loadset(files(p).name, dataDir);

% Take only bad segments that were automatically detected
V_Rejected_Sample_Range = [];
rejDataNew = [];
rejData = [];

rejData =  [rejData ([EEG.event(:,(ismember({EEG.event.type}, {'auto_start'}))).latency])';];
rejData =  [rejData ([EEG.event(:,(ismember({EEG.event.type}, {'auto_end'}))).latency])';];
rowNb = find(ismember({EEG.event.type}, {'auto_start'})==1);

%Mark manually bad segments but do not reject
pop_eegplot_adjust(EEG, 1, 1, 0)

figure; topoplot([],EEG.chanlocs, 'style', 'blank',  'electrodes', 'numpoint', 'chaninfo', EEG.chaninfo);

fig = figure;
pop_spectopo(EEG, 1, [0 length(EEG.data)], 'EEG' , 'freqrange',[1 45], 'title', [char(files(p).name),  '_' ,char(day)])

rows = ones(length(rowNb),1);

% display intervals on the command window

auto_seg = EEG.event(:,(ismember({EEG.event.type}, {'auto_start'})),:);

disp('Seconds of marked bad segments: ')

for k = 1:length(auto_seg) 
    test_seg = [auto_seg(k).latency : (auto_seg(k).duration + auto_seg(k).latency)];
    disp([num2str(round(auto_seg(k).latency./500)), ' sec ', '; ', 'row number: ', num2str(k), ]);
end

keep_all = input('Would you like to keep all suggested cuts to the data? Please answer 1 - Yes or 0 - No: ') ;% 1 or 0

if keep_all ==1
    rejDataNew = rejData;
else
    numrows = ['The total amount of detected bad segments: ',num2str(length(rowNb))];
    disp(numrows)
    delete_all_rows = input('Would you like to delete all suggested cuts to the data? Please answer 1 - Yes or 0 - No: ');% input which rows you want to keep as marked bad segments, 1st row should also be kept!
    if delete_all_rows == 1
        EEG.event(:,(ismember({EEG.event.type}, {'auto_start'})),:)=[];
        EEG.event(:,(ismember({EEG.event.type}, {'auto_end'})),:)=[];
        
    else
        num_rows_remove = str2num(input('Where (rows in numbers - above) is good data falsely detected as bad? Please input their row number from above!  ', 's')); % input which rows you want to keep as marked bad segments, 1st row should also be kept!
        
        for isb = 1:length(num_rows_remove)
            rows(num_rows_remove(isb)) = 0;
        end
        rows = logical(rows);
        rejDataNew = rejData(rows,:); %TO DO
        auto_seg = auto_seg(rows);
        EEG.event(:,(ismember({EEG.event.type}, {'auto_start'})),:)=[];
        EEG.event(:,(ismember({EEG.event.type}, {'auto_end'})),:)=[];
        if isempty(auto_seg)
        else
            for u=1:length(auto_seg)
                l=length(EEG.event);
                EEG.event(l+1) = auto_seg(u);
            end
        end
%         EEG.event(l+1,:) = auto_seg(:,:);
    end
end


% If any segments were manually rejected, update the rejDataNew and
% EEG.event
Rejected_Sample=any(V_Rejected_Sample_Range);
if any(Rejected_Sample)
    for r = 1:size(V_Rejected_Sample_Range,1)
        l = size(rejDataNew,1);
        n = length(EEG.event);
        rejDataNew(l+1,1) = V_Rejected_Sample_Range(r,1) %beginning of marked manually marked bad segment
        rejDataNew(l+1,2) = V_Rejected_Sample_Range(r,2) %end of manually maked bad segment
        
        EEG.event(n+1).latency = V_Rejected_Sample_Range(r,1);
        EEG.event(n+1).duration = V_Rejected_Sample_Range(r,2)-V_Rejected_Sample_Range(r,1); % + half a second to compensate for the start extension and + 1 second to compensate for the end extension
        EEG.event(n+1).type = 'manual_start';
        
    end
end
if any(rejDataNew)
    [~,idx] = sort(rejDataNew(:,1)); % sort just the first column
    rejDataNew = rejDataNew(idx,:);
    if rejDataNew(1,1) < 1
        rejDataNew(1,1) = 1;
    end
      
        isb = 1;
        while size(rejDataNew,1) > 1 && isb < size(rejDataNew,1)
            prev_minus_next = rejDataNew(isb+1,1) - rejDataNew(isb,2);
            if (prev_minus_next) < 3 * EEG.srate
                rejDataNew(isb,2) = rejDataNew(isb+1,2);
                rejDataNew(isb+1,:) = [];
                isb = isb - 1;
            end
            isb = isb + 1 ;
        end
        
        no_seg = ['The total amount of detected bad segments after the cleaning: ',num2str(size(rejDataNew,1))];
        warning(no_seg)
        % add manual marks to EEG.event
        EEG_rej = pop_select(EEG, 'nopoint', [rejDataNew(:,1) rejDataNew(:,2)]); %reject data for ICA. TO DO: for manual artifact detection
        
        fig = figure;
        pop_spectopo(EEG, 1, [0 length(EEG.data)], 'EEG' , 'freqrange',[1 45], 'title', [char(files(p).name),  '_' ,char(day)]), set(fig, 'Visible', 'off');
        
        % This part changes the color of frontal channels to red
        colors = hot(64);
        axhandle = gca;
        lines = findobj(axhandle,'type','line');
        
        set(lines(length(lines)-1:length(lines)),'color',colors(37,:));
        %savefigure
        export_fig([dataDir, 'spectograms'], '-pdf', '-append', fig); % h is
        
        
        close Figure 2
        %Check if set is present in manual inspection folder and delete old
        
        % calculate ICA weights only on good segments (remove bad ones)
        EEG_ica = pop_runica(EEG_rej, 'extended',1,'interupt','on');
        EEG_28ch = pop_editset(EEG, 'icachansind', 'EEG_ica.icachansind', 'chanlocs', '{EEG_ica.chanlocs EEG_ica.chaninfo EEG_ica.urchanlocs }', 'icaweights', 'EEG_ica.icaweights', 'icasphere', 'EEG_ica.icasphere');
    
else
        fig = figure;
        pop_spectopo(EEG, 1, [0 length(EEG.data)], 'EEG' , 'freqrange',[1 45], 'title', [char(files(p).name),  '_' ,char(day)]), set(fig, 'Visible', 'off');
        
        % This part changes the color of frontal channels to red
        colors = hot(64);
        axhandle = gca;
        lines = findobj(axhandle,'type','line');
        
        set(lines(length(lines)-1:length(lines)),'color',colors(37,:));
        %savefigure
        export_fig([dataDir, 'spectograms'], '-pdf', '-append', fig); % h is
        close Figure 2
        
        warning('You are rejecting 0 data!')
        EEG_ica = pop_runica(EEG, 'extended',1,'interupt','on');
        EEG_28ch = pop_editset(EEG, 'icachansind', 'EEG_ica.icachansind', 'chanlocs', '{EEG_ica.chanlocs EEG_ica.chaninfo EEG_ica.urchanlocs }', 'icaweights', 'EEG_ica.icaweights', 'icasphere', 'EEG_ica.icasphere');
end
%Save set in pre_ICA folder and remove iteration number
saveName = [files(p).name(1:14), 'preICA']

EEG_28ch = pop_saveset(EEG_28ch, 'filename', saveName,'filepath', saveDir);


