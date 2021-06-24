% (C) Elena Cesnaite, email: e.cesnaite@gmail.com, page: https://www.researchgate.net/profile/Elena-Cesnaite

% This code was created to analyze data described in a paper:'One-week escitalopram intake shifts excitation-inhibition balance 
% in healthy female brain: implications for individual cortical responsivity to SSRIs' by Zsido & Molloy et al. 
% The code cleans data by rejecting ICA components that are associated with noise. 
% Last updated 22.06.2021

clc; clear all; close all

eeglab

dataDir = '';
saveDir = '';
probDir = '';
textDir = '';
cntu = [];

% adjust some eeglab options
files = dir([dataDir, '*.set']);
pop_editoptions( 'option_storedisk', 0, 'option_savetwofiles', 1, 'option_saveversion6', 1, 'option_single', 1, 'option_memmapdata', 0, 'option_eegobject', 0, 'option_computeica', 1, 'option_scaleicarms', 1, 'option_rememberfolder', 1, 'option_donotusetoolboxes', 0, 'option_checkversion', 1, 'option_chat', 0);

Fig1 = figure('Name','electrode positions');
Fig2 = figure('Name','spectrum of before any cleaning');
Fig3 = figure('Name',sprintf('fastICA Results - plot %d',1), 'units','normalized','outerposition',[0 0 1 1]);
Fig4 = figure('Name','spectrum of active/active&passive cleaned EEG'); set(gcf,'position',[10,10,1200,500]);
pause(0.5)
%%
for isa = 1:length(files)
    
    %Check if set is present in before ICA folder and continue if yes
    savedfiles=dir([saveDir, '*.set']);
    name = extractBetween(files(isa).name,"","_preICA.set")
    s_name = extractBetween({savedfiles.name},"","_prep.set");
% %     
     if sum(strcmp(s_name, char(name))) == 1
         continue
     else
        
        cntu = 'n';
        while cntu == 'n'
            clf(Fig1,'reset'); clf(Fig2,'reset'); clf(Fig3,'reset');clf(Fig4,'reset'); pause(0.5);
            EEG = pop_loadset(strcat(name,'_preICA.set'), dataDir); 

            EEG = eeg_checkset(EEG);
        % reject artifacts and run power spectrum on clean data
            events = EEG.event(:,(ismember({EEG.event.type}, [{'auto_start'}, {'manual_start'}])),:);
            if isempty(events)
                EEG_rej= EEG;
            else
                rejData(:,1) = [events.latency];
                rejData(:,2) = [events.latency] + [events.duration];
                rejData;
           
                EEG_rej = pop_select(EEG, 'nopoint', [rejData(:,1) rejData(:,2)]); %reject data for ICA.
            end

            set(0, 'CurrentFigure', Fig1);
            topoplot([],EEG_rej.chanlocs, 'style', 'blank',  'electrodes', 'numpoint', 'chaninfo', EEG_rej.chaninfo);
            set(0, 'CurrentFigure', Fig2);
            plot_spec(EEG_rej.data',EEG_rej.srate, 45), title('original power spectrum');
            yLimit_spec = get(gca,'YLim');
                        
            %Calculating for each component for how much variance it accounts in
            %the data
            disp('Calculating variance accounted for')
            pvar = NaN(length(EEG_rej.icaweights),1);
            for n_pvar = 1:length(EEG_rej.icaweights)
                pvar(n_pvar) =  pvaf_PS(EEG_rej, n_pvar);
            end
            
            [bad_cmp_acitive,bad_cmp_passive] = select_components(EEG_rej.icawinv./std(EEG_rej.icawinv),EEG_rej.icaact,EEG_rej.srate,EEG_rej.chanlocs, pvar,'fig',{Fig3},'data',EEG_rej.data);
                        
            noisy_ds = input('NOISY? YES=''y'' - NO=''n'': ','s');
            if noisy_ds == 'y'
                fileID = fopen(fullfile(probDir,'EEG_noisy_subjects.txt'),'a+'); %do not overwrite
                fprintf(fileID,'\n %s %s %12s',datestr(datetime),char(name));
                fclose(fileID);
                
                fileID = fopen(fullfile(textDir,[char(name),'.txt']),'a+'); %do not overwrite
                fprintf(fileID,'\n %s %12s', datestr(datetime),';Noisy power spectrum. ');
                fclose(fileID);
            end
            
            fileID = fopen(fullfile(probDir, 'EEG_rejected_ICA_components.txt'),'a+'); %do not overwrite
            fprintf(fileID,'\n %s %s %12s %s',datestr(datetime),char(name),';Number of rejected components: ', num2str(numel(bad_cmp_acitive)));
            fclose(fileID);
            
            fileID = fopen(fullfile(textDir,[char(name),'.txt']),'a+'); %do not overwrite
            fprintf(fileID,'\n %s %12s %s', datestr(datetime),';Number of rejected components: ', num2str(numel(bad_cmp_acitive)));
            fclose(fileID);
            rejcomp = find([bad_cmp_acitive]);
            text = ['The number of rejected active components: ', num2str(length(rejcomp))];
            disp(text)
            
            %% clean
            EEG_a = EEG;
            EEG_a_rej = EEG_rej;
            
            
            EEG_a_p = EEG;
            EEG_a_p_rej = EEG_rej;
            if numel(bad_cmp_acitive) %if any active component is selected
                EEG_a = pop_subcomp(EEG_a, bad_cmp_acitive, 0);
                EEG_a_rej = pop_subcomp(EEG_a_rej, bad_cmp_acitive, 0);
            end
            if numel(bad_cmp_passive) %if any passive component is selected
                bad_cmp_a_p = [bad_cmp_acitive;bad_cmp_passive];
                EEG_a_p = pop_subcomp(EEG_a_p, bad_cmp_a_p, 0); % EEG_a_p
                EEG_a_p_rej = pop_subcomp(EEG_a_p_rej, bad_cmp_a_p, 0);
            end
            
            % Plot time series
            
            EEG = EEG_a; % EEG is now the whole EEG dataset with active components only
            
            %% plotting the spectrum of cleaned with active and passive components
            
            set(0, 'CurrentFigure', Fig4);
            subplot(1,2,1), plot_spec(EEG_a_rej.data',EEG_a_rej.srate,45), title('cleared with active components'), ylim(yLimit_spec)
            if numel(bad_cmp_passive)
                subplot(1,2,2), plot_spec(EEG_a_p_rej.data',EEG_rej.srate,45), title('cleared with active and passive components'),ylim(yLimit_spec)
            else
                subplot(1,2,2), plot_spec(EEG_a_p_rej.data',EEG_rej.srate,45), title('ORIGINAL SPECTRUM'),ylim(yLimit_spec)
            end
            %
            
            clear ind_F_temp X_spec_a X_spec_a_p F
            clear EEG_a EEG_a_p EEG_a_rej EEG_a_p_rej
            
            EEG = eeg_checkset(EEG);
            cntu = input('SAVE? NO=''n'' - YES=any: ','s');
            
            if cntu == 'n'
                %close all
            else
                EEG.reject.icarejedcomp = bad_cmp_acitive;
                EEG.reject.icamusclecomp = bad_cmp_passive;
                EEG = pop_saveset( EEG, 'filename',[char(name), '_prep'],'filepath',saveDir);
                %close all
                
            end
 
            clearvars -except dataDir saveDir  probDir cntu textDir namesICA isa Fig1 Fig2 Fig3 Fig4 name s_name files
            
        end
     end
end
