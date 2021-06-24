% (C) Elena Cesnaite, email: e.cesnaite@gmail.com, page: https://www.researchgate.net/profile/Elena-Cesnaite

% This code was created to analyze data described in a paper:'One-week escitalopram intake shifts excitation-inhibition balance 
% in healthy female brain: implications for individual cortical responsivity to SSRIs' by Zsido & Molloy et al. 
% This code removes the mastoid electrode
% Last updated 22.06.2021

clc; clear all; close all

eeglab

dataDir = '';
saveDir = '';

files = dir([dataDir, '*.set']);

for isa = 1:length(files)
    
    name = extractBetween(files(isa).name,"",".set")
    save_name = strcat(name, '_final');

    EEG = pop_loadset(files(isa).name, dataDir); 
    EEG = eeg_checkset(EEG);
    
     %remove the channel from the data
     EEG_2 = pop_select(EEG, 'nochannel', {'M2'});
     
     EEG_2 = pop_saveset(EEG_2, 'filename', char(save_name),'filepath', saveDir); %incl

end
