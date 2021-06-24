% (C) Elena Cesnaite, email: e.cesnaite@gmail.com, page: https://www.researchgate.net/profile/Elena-Cesnaite

% This code was created to analyze data described in a paper:'One-week escitalopram intake shifts excitation-inhibition balance 
% in healthy female brain: implications for individual cortical responsivity to SSRIs' by Zsido & Molloy et al. 
% This code calculates ICA weights of a common average reference data. 
% Last updated 22.06.2021

clc; clear all; close all
eeglab

dataDir = '';
saveDir = ''
files = dir([dataDir, '*.set']);

for isa = 1:length(files)
    
    EEG = pop_loadset(files(isa).name, dataDir); 
    
    EEG = eeg_checkset(EEG);
    
    %re-reference to common average
    EEG.data = EEG.data-repmat(mean(EEG.data),size(EEG.data,1),1);
    EEG.icaact = []
    EEG.icawinv = []
    EEG.icasphere = []
    EEG.icaweights = []
    EEG.icachansind = []
    
    events = EEG.event(:,(ismember({EEG.event.type}, [{'auto_start'}, {'manual_start'}])),:);
    if isempty(events)
        EEG_rej= EEG;
    else
        rejData = [[events.latency]', [[events.latency] + [events.duration]]'];
        EEG_rej = pop_select(EEG, 'nopoint', [rejData(:,1) rejData(:,2)]); %reject data for ICA.
    end
     
    EEG_ica = pop_runica(EEG_rej, 'extended',1,'interupt','on');
    EEG.icachansind = EEG_ica.icachansind;
    EEG.icaweights = EEG_ica.icaweights;
    EEG.icasphere = EEG_ica.icasphere;
    EEG = pop_saveset(EEG, 'filename',files(isa).name(1:end-4),'filepath', saveDir); 
    
end