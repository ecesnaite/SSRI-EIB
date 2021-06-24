% (C) Elena Cesnaite, email: e.cesnaite@gmail.com, page: https://www.researchgate.net/profile/Elena-Cesnaite

% This code was created to analyze data described in a paper:'One-week escitalopram intake shifts excitation-inhibition balance 
% in healthy female brain: implications for individual cortical responsivity to SSRIs' by Zsido & Molloy et al. 
% The code rejects bad data segments and calculates ICA weights for further pre-processing. 
% Last updated 22.06.2021

clear; clc; close all

eeglab

dataDir = ''
saveDir = ''
inspectDir = '' % datasets that need manual inspection will be saved here
textDir = ''

files=dir([dataDir]);
subjects={files(3:end).name}
days = {'BL', 'D1', 'D7'}

for p = 1:length(subjects)%
    
    for j = 1:length(days)
        dayDir = strcat(dataDir, subjects(p), '/', days(j), '/')
        dayDir = char(dayDir)
        vhdr_files = dir([dayDir, '*.vhdr'])
        
        if isempty(vhdr_files)
            clear EEG
            fileID = fopen('Missing_Data_List.txt','a+'); %do not overwrite
            fprintf(fileID,'\n %s %6s %s %12s',datestr(datetime),char(subjects(p)), strcat(' iteration_', num2str(p)),'; There are no files in the folder!');
            fclose(fileID);
            continue
        end
        
        EEG = pop_loadbv(dayDir, vhdr_files.name);
        
        % Inspect data for sampling rate and store the outcome if it's other than 1000 Hz%
        if EEG.srate~=1000
            fileID = fopen('EEG_problems_detected.txt','a+'); %do not overwrite
            fprintf(fileID,'\n %s %6s %s %12s %s',datestr(datetime),EEG.comments, strcat(' iteration_', num2str(p)),'; Sampling rate is not 1000 Hz! It is: ', num2str(EEG.srate));
            fclose(fileID);
            continue
        end
        % personalized txt file
        fileID = fopen(fullfile(textDir,[char(subjects(p)),'.txt']),'a+'); %do not overwrite
        fprintf(fileID,'\n %s %s %s %12s %s', datestr(datetime), strcat(' iteration_', num2str(p)),char(days(j)),'; Sampling rate: ', num2str(EEG.srate));
        fclose(fileID);
        
        % Band-pass filter combined with the notch
        [b1,a1] = butter(2,[1 45]/(EEG.srate/2));
        [b2,a2] = butter(2,[49 51]/(EEG.srate/2),'stop');
        
        EEG.data = filtfilt(conv(b1,b2),conv(a1,a2),double(EEG.data)')';
        
        % Downsample data to 500 Hz %
        EEG = pop_resample(EEG,EEG.srate/2);
        
        if length(EEG.event) < 5 
            fileID = fopen('EEG_no_markers.txt','a+'); %do not overwrite
            fprintf(fileID,'\n %s %6s %s %12s',datestr(datetime),EEG.comments, strcat(' iteration_', num2str(p)),'; There are no markers in the file!');
            fclose(fileID);
            EEG = pop_saveset(EEG, 'filename', [char(subjects(p)), '_' ,char(days(j)), '_manual_inspection_', num2str(p)],'filepath', inspectDir);
            continue
        end
        
        % set marker positions
        marker_position = EEG.event(:,(ismember({EEG.event.type}, {'S210'})),:);
        
        % cut data between start and end markers  + 2 sec     
        EEG = pop_select(EEG, 'point', [marker_position(1).latency, (marker_position(end).latency + EEG.srate * 2)]);
        
        % find S1
        
         s1_position = EEG.event(:,(ismember({EEG.event.type}, {'S  1'})),:);
         %break start S 210 +2sec before the S 1 
         break_start = (EEG.event(:,(ismember([EEG.event.bvmknum], [(s1_position(1).bvmknum) - 1])),:).latency) + EEG.srate * 2;
         break_end = s1_position(end).latency;
         
         %cut out the break from the data
         EEG = pop_select(EEG, 'nopoint', [break_start, break_end]);
         % add 'boundary_removed' into the EEG.event.code structure as a
         % comment
         boundary = EEG.event(:,ismember({EEG.event.type}, {'boundary'}),:);
         EEG.event(:,ismember([EEG.event.latency], [boundary(2).latency]),:).code = 'break_removed';
        
         if length(EEG.data)/(EEG.srate*60) < 11 %inspect if resting state is less than 11min
             % personalized txt file
             fileID = fopen(fullfile(textDir,[char(subjects(p)),'.txt']),'a+'); %do not overwrite
             fprintf(fileID,'\n %s %s %s %12s %s', datestr(datetime), strcat(' iteration_', num2str(p)),char(days(j)),'; Resting state data is less than 11 min! It is: ', num2str(length(EEG.data)/(EEG.srate*60)));
             fclose(fileID);
             
             fileID = fopen('EEG_less_than_11min.txt','a+'); %do not overwrite
             fprintf(fileID,'\n %s %6s %s %12s %s',datestr(datetime), EEG.comments, strcat(' iteration_', num2str(p)),'; Resting state data is less than 11 min! It is: ', num2str(length(EEG.data)/(EEG.srate*60)));
             fclose(fileID);
             EEG = pop_saveset(EEG, 'filename', [char(subjects(p)), '_' ,char(days(j)), '_manual_inspection_short_dataset_', num2str(p)],'filepath', inspectDir);
             continue
        else
            % personalized txt file
            fileID = fopen(fullfile(textDir,[char(subjects(p)),'.txt']),'a+'); %do not overwrite
            fprintf(fileID,'\n %s %s %s %s %12s %s', datestr(datetime), strcat(' iteration_', num2str(p)), char(days(j)),'; Length of recording in minutes: ', num2str(length(EEG.data)/(EEG.srate*60)));
            fclose(fileID);
            
        end
        
        % Load Channel locations
        EEG = pop_chanedit(EEG, 'lookup','');
        figure; topoplot([],EEG.chanlocs, 'style', 'blank',  'electrodes', 'numpoint', 'chaninfo', EEG.chaninfo);
        
        % Exclude EKG and eye movement channels before rejecting bad channels.
        indx_chan = ismember({EEG.chanlocs.labels}, {'VEOGu', 'VEOGo', 'HEOGr', 'HEOGl'});
        
        % Excluded channels should be saved
        if sum(indx_chan)~= 4
            fileID = fopen('EEG_channels.txt','a+'); %do not overwrite
            fprintf(fileID,'\n %s %6s %s %12s %s',datestr(datetime), EEG.comments, strcat(' iteration_', num2str(p)), '; No VEOGu, VEOGo, HEOGr, and/or HEOGl were found.');
            fclose(fileID);
            EEG = pop_saveset(EEG, 'filename', [char(subjects(p)), '_' ,char(days(j)), '_manual_inspection_missing_channels_', num2str(p)],'filepath', inspectDir);
            continue
        else
            EEG_28ch = pop_select(EEG, 'nochannel', {'VEOGu', 'VEOGo', 'HEOGr', 'HEOGl'});
            EEG_28ch.EOG_channels = pop_select(EEG, 'channel', {'VEOGu', 'VEOGo', 'HEOGr', 'HEOGl'});
        end

    % Mark artifacts that exceed threshold (3stdev from mean) for low
    % frequencies 1-15 Hz, 40 mV is the threshold for higher frequencies
    % 15-45 Hz
    
    [EEG_28ch rejData] = trimOutlier_adjust_Eoin(EEG_28ch, 80, 40, 500);   %trimOutlier_adjust also marks 2s at the beginning of recording
    
    auto_seg = EEG_28ch.event(:,(ismember({EEG_28ch.event.type}, {'auto_start'})),:);
    
    fileID = fopen(fullfile(textDir,[char(subjects(p)),'.txt']),'a+'); 
    fprintf(fileID,'\n %s %s %s %12s %s', datestr(datetime), strcat(' iteration_', num2str(p)),'; Number of artifacts detected by automatic algorythm: ', num2str(length(auto_seg)), '. Length in seconds: ', num2str(sum([auto_seg.duration])/500));
    fclose(fileID);
    
    
    if sum([auto_seg.duration])/EEG_28ch.srate <= 30 %trim_outlier - if length of auto detected bad segments is longer than 30 seconds - flag for manual review
       if rejData == 0
           EEG_rej = EEG_28ch;
       else
            EEG_rej = pop_select(EEG_28ch, 'nopoint', [rejData(:,1) rejData(:,2)]); %reject data for ICA. 
       end
    else
        fileID = fopen('EEG_bad_segments.txt','a+'); % if noisy segments (>30sec) are found, store it and continue
        fprintf(fileID,'\n %s %6s %s %12s %s',datestr(datetime), EEG.comments, strcat(' iteration_', num2str(p)), '; Bad segments were found that doesnt match Christians files. Load data again and inspect time series. Seconds detected:');
            for isb = 1:length(auto_seg)           
                fprintf(fileID,'\n %s %d %s %d',' row No.', isb, 'sec:', round(auto_seg(isb).latency./500));           
            end
        
        fclose(fileID);
        EEG = pop_saveset(EEG_28ch, 'filename', [char(subjects(p)), '_' ,char(days(j)), '_manual_inspection_bad_segments_', num2str(p)],'filepath', inspectDir);
        continue
         
    end
  % calculate ICA weights only on good segments (remove bad ones)
    EEG_ica = pop_runica(EEG_rej, 'extended',1,'interupt','on');
    EEG_28ch = pop_editset(EEG_28ch, 'icachansind', 'EEG_ica.icachansind', 'chanlocs', '{EEG_ica.chanlocs EEG_ica.chaninfo EEG_ica.urchanlocs }', 'icaweights', 'EEG_ica.icaweights', 'icasphere', 'EEG_ica.icasphere');
    
    EEG_28ch = pop_saveset(EEG_28ch, 'filename',[char(subjects(p)), '_' ,char(days(j)), '_preICA'],'filepath', saveDir); %include 3 channel data
     
   fig = figure; 
    pop_spectopo(EEG_rej, 1, [0 length(EEG_rej.data)], 'EEG' , 'freqrange',[1 45], 'title', [char(subjects(p)),  '_' ,char(days(j))]), set(fig, 'Visible', 'off');
    
    % This part changes the color of frontal channels to red
    colors = hot(64);
    axhandle = gca;
    lines = findobj(axhandle,'type','line');
    
    set(lines(length(lines)-1:length(lines)),'color',colors(37,:));
    %savefigure
    export_fig([inspectDir, 'spectograms'], '-pdf', '-append', fig); % h is
    
    
    close Figure 2
    
    end
end


%% Create another loop for the ICA weights and plotting


