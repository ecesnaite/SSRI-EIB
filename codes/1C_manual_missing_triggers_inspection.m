% (C) Elena Cesnaite, email: e.cesnaite@gmail.com, page: https://www.researchgate.net/profile/Elena-Cesnaite

% This code was created to analyze data described in a paper:'One-week escitalopram intake shifts excitation-inhibition balance 
% in healthy female brain: implications for individual cortical responsivity to SSRIs' by Zsido & Molloy et al. 
% . 
% Last updated 22.06.2021

clc; clear all; close all

eeglab

% Load data %
saveDir = ''
dataDir = ''
textDir = ''

files=dir([dataDir, '*.set']);

%Get the right dataset from the Manual_Inspection folder based on the
%iteration number in the main script

for p = 1:length(files)
    
    %Load Data
    EEG = pop_loadset(files(p).name, dataDir);
    
    % Take only bad segments that were automatically detected
    V_Rejected_Sample_Range = [];
    fullname = files(p).name;
    subject = fullname(1:10);
    day = fullname(12:13);
    
    if length(EEG.data)/(EEG.srate*60) < 11 %inspect if resting state is less than 11min
        % personalized txt file
        fileID = fopen(fullfile(textDir,[char(subject),'.txt']),'a+'); %do not overwrite
        fprintf(fileID,'\n %s %s %s %12s %s', datestr(datetime), day,'; Resting state data is less than 11 min! It is: ', num2str(length(EEG.data)/(EEG.srate*60)));
        fclose(fileID);
        
        fileID = fopen('EEG_less_than_11min.txt','a+'); %do not overwrite
        fprintf(fileID,'\n %s %6s %s %12s %s',datestr(datetime), EEG.comments, day, '; Resting state data is less than 11 min! It is: ', num2str(length(EEG.data)/(EEG.srate*60)));
        fclose(fileID);
        
        error('the length of your data is too short!')
    else
        % personalized txt file
        fileID = fopen(fullfile(textDir,[char(subject),'.txt']),'a+'); %do not overwrite
        fprintf(fileID,'\n %s %s %s %12s %s', datestr(datetime), day,'; Length of recording in minutes: ', num2str(length(EEG.data)/(EEG.srate*60)));
        fclose(fileID);
        
    end
    
    
    % Load Channel locations
    EEG = pop_chanedit(EEG, 'lookup','/nobackup/kiribati2/EEG/Scripts/eeglab14_1_1b/plugins/dipfit2.3/standard_BESA/standard-10-5-cap385.elp');
    figure; topoplot([],EEG.chanlocs, 'style', 'blank',  'electrodes', 'numpoint', 'chaninfo', EEG.chaninfo);
    
    % Exclude EKG and eye movement channels before rejecting bad channels.
    indx_chan = ismember({EEG.chanlocs.labels}, {'VEOGu', 'VEOGo', 'HEOGr', 'HEOGl'});
    
    % Excluded channels should be saved
    if sum(indx_chan)~= 4
        fileID = fopen('EEG_channels.txt','a+'); %do not overwrite
        fprintf(fileID,'\n %s %6s %s %12s %s',datestr(datetime), EEG.comments, strcat(' iteration_', num2str(p)), '; No VEOGu, VEOGo, HEOGr, and/or HEOGl were found.');
        fclose(fileID);
        
        error('No peripheral channels were found or not all!')
        
    else
        EEG.EOG_channels = pop_select(EEG, 'channel', {'VEOGu', 'VEOGo', 'HEOGr', 'HEOGl'});
        EEG = pop_select(EEG, 'nochannel', {'VEOGu', 'VEOGo', 'HEOGr', 'HEOGl'});
    end
    
    
    fileID = fopen(fullfile(textDir,[char(subject),'.txt']),'a+');
    fprintf(fileID,'\n %s %s %12s', datestr(datetime), day,'; This subject was inspected manually because of the missing triggers. ');
    fclose(fileID);
    
    t = ['The total minute length of the recording is: ', num2str(length(EEG.data)/(EEG.srate*60))];
    disp(t)
    
    fig = figure;
    pop_spectopo(EEG, 1, [0 length(EEG.data)], 'EEG' , 'freqrange',[1 45], 'title', [char(files(p).name),  '_' ,char(day)]);%, set(fig, 'Visible', 'off');
    
    pop_eegplot_adjust(EEG, 1, 1, 0) % keep last imput as 0 for marking bad segments without creating a new dataset
    
    c = input('would you like to continue?', 's')
    
    if c == 'n'
        return
    end
       
    % If any segments were manually rejected, update the rejDataNew and
    Rejected_Sample=any(V_Rejected_Sample_Range);
    if any(Rejected_Sample)
        rejDataNew = []
        rejDataNew = V_Rejected_Sample_Range %beginning of marked manually marked bad segment
        EEG.event = []

        for r = 1:size(V_Rejected_Sample_Range,1)
         
            EEG.event(r).latency = V_Rejected_Sample_Range(r,1);
            EEG.event(r).duration = V_Rejected_Sample_Range(r,2)-V_Rejected_Sample_Range(r,1); % + half a second to compensate for the start extension and + 1 second to compensate for the end extension
            EEG.event(r).type = 'manual_start';
            
        end
       
        % add manual marks to EEG.event
        EEG_rej = pop_select(EEG, 'nopoint', [rejDataNew(:,1) rejDataNew(:,2)]); %reject data for ICA
        
        t2 = ['The total minute length of the recording after cleaning is: ', num2str(length(EEG_rej.data)/(EEG_rej.srate*60))];
        warning(t2)
              
        
        fig = figure;
        pop_spectopo(EEG_rej, 1, [0 length(EEG_rej.data)], 'EEG' , 'freqrange',[1 45], 'title', char(files(p).name));
        
        % This part changes the color of frontal channels to red
        colors = hot(64);
        axhandle = gca;
        lines = findobj(axhandle,'type','line');
        
        set(lines(length(lines)-1:length(lines)),'color',colors(37,:));
        %savefigure
        export_fig([dataDir, 'spectograms'], '-pdf', '-append', fig); % h is
        c = input('would you like to continue?', 's')
        
        if c == 'n'
            return
        end
        close Figure 2
        %Check if set is present in manual inspection folder and delete old
        
        % calculate ICA weights only on good segments (remove bad ones)
        EEG_ica = pop_runica(EEG_rej, 'extended',1,'interupt','on');
        EEG = pop_editset(EEG, 'icachansind', 'EEG_ica.icachansind', 'chanlocs', '{EEG_ica.chanlocs EEG_ica.chaninfo EEG_ica.urchanlocs }', 'icaweights', 'EEG_ica.icaweights', 'icasphere', 'EEG_ica.icasphere');
        
    else
        fig = figure;
        pop_spectopo(EEG, 1, [0 length(EEG.data)], 'EEG' , 'freqrange',[1 45], 'title', char(files(p).name)), set(fig, 'Visible', 'off');
        
        % This part changes the color of frontal channels to red
        colors = hot(64);
        axhandle = gca;
        lines = findobj(axhandle,'type','line');
        
        set(lines(length(lines)-1:length(lines)),'color',colors(37,:));
        %savefigure
        export_fig([dataDir, 'spectograms'], '-pdf', '-append', fig); % h is
        close Figure 2
        
        warning('You are rejecting 0 data!')
        pop_eegplot_adjust(EEG, 1, 1, 1)
        c = input('would you like to continue?', 's')
        
        if c == 'n'
            return
        end
        
        EEG_ica = pop_runica(EEG, 'extended',1,'interupt','on');
        EEG = pop_editset(EEG, 'icachansind', 'EEG_ica.icachansind', 'chanlocs', '{EEG_ica.chanlocs EEG_ica.chaninfo EEG_ica.urchanlocs }', 'icaweights', 'EEG_ica.icaweights', 'icasphere', 'EEG_ica.icasphere');
    end
    %Save set in pre_ICA folder and remove iteration number
    saveName = [files(p).name(1:14), 'preICA']
    
    EEG = pop_saveset(EEG, 'filename', saveName,'filepath', saveDir);
end



