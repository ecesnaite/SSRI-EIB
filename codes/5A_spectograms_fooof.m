% (C) Elena Cesnaite, email: e.cesnaite@gmail.com, page: https://www.researchgate.net/profile/Elena-Cesnaite

% This code was created to analyze data described in a paper:'One-week escitalopram intake shifts excitation-inhibition balance 
% in healthy female brain: implications for individual cortical responsivity to SSRIs' by Zsido & Molloy et al. 
% This code calculates power spectral density for each channel which will
% be used as in input data for the FOOOF algorythm.
% Last updated 22.06.2021

clc; clear all; close all

eeglab

aftDir = '' % dataDir
saveDir = ''

files = dir([aftDir, '*.set']);

for isb = 1:length(files)
    
    EEG = pop_loadset(files(isb).name, aftDir);
    name = extractBefore(files(isb).name, '_prep')
        
   events = EEG.event(:,(ismember({EEG.event.type}, [{'auto_start'}, {'manual_start'}])),:);%Remove bad segments previously marked by us and 

    if isempty(events)
        EEG_rej= EEG;
    else
        rejData(:,1) = [events.latency];
        rejData(:,2) = [events.latency] + [events.duration];
        rejData;
        EEG_rej = pop_select(EEG, 'nopoint', [rejData(:,1) rejData(:,2)]);
    end
   
       [spect freq] = 6B_el_plot_spec_SSRI(EEG_rej.data',EEG_rej.srate,45)
       spect = spect';
       psd(isb).spect = spect
       psd(isb).freq = freq
       psd(isb).ID = name
       fig = figure('Visible', 'off');
       plot(freq, log10(spect)), title(name)
        export_fig([saveDir, 'spectograms_4'], '-pdf', '-append', fig); 
    clear rejData EEG EEG_rej spect freq name %UPDATE
end

save([saveDir, 'psd_and_freq_for_FOOOF_4'], 'psd')








