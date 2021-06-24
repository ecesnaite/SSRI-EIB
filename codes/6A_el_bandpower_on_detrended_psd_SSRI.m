% (C) Elena Cesnaite, email: e.cesnaite@gmail.com, page: https://www.researchgate.net/profile/Elena-Cesnaite

% This code was created to analyze data described in a paper:'One-week escitalopram intake shifts excitation-inhibition balance 
% in healthy female brain: implications for individual cortical responsivity to SSRIs' by Zsido & Molloy et al. 
% This code calculates alpha band parameters on detrended psd
% Last updated 22.06.2021

clear all, close all

psd_all=load('') %original PSD

dataDir = ''
saveDir = ''
files= dir([dataDir,'0*'])

f_orig = [0:0.25:45]; %original sampling frequency
fs = [1:0.25:40]; %sampling frequency based on which FOOOF estimates were made between 1 - 40 Hz

for p = 1:length(psd_all.psd)
    
    id = psd_all.psd(p).ID
    original = psd_all.psd(strcmp({psd_all.psd.ID},id));
    
    data = load([dataDir, strcat(id, '.mat')]); %1 - background params (offset and exponent (slope)); 2 - peak params (central freq, ampl, bandwi); 3 - r_squared; 5 - gaussian_params
    for ch = 1:27 % 27 channels
        rsq = data.results{1,1}{ch,3};
        
        if rsq >= 0.8
            aperiodic = data.results{1,1}{ch,1};
            offset =aperiodic(1);
            slope = aperiodic(2);
            
            logL = offset - slope*log10(fs);
            L = 10.^logL;
         
            % subtract the slope from the original psd %
            detrend = original.spect(ch,f_orig>=1 & f_orig<=40) - L;
            
            %find peaks that exceed 0.05 uV
            [pow,ifreq,width,~,bounds] = findpeaks_adjusted_new_SSRI(detrend,fs,'MinPeakProminence',0.05,'Annotate','extents','WidthReference','halfheight');% min peak for power threshold
            
            any_alp = (ifreq >= 7 & ifreq <= 13);
            if any(any_alp)
                if sum(any_alp)==1 % if one alpha peak is present
                    
                    freq = ifreq(any_alp);
                    range = bounds(any_alp,:);
                    
                    if diff(range) > 6      %whenever the width of the peak is more than 3Hz, set it to 3
                        range = [freq-3, freq+3]; 
                    end
                    
                    alpha_broad_power = sum(detrend(fs>=round(range(1),2) & fs<=round(range(2),2)))*0.25;
                    
                    alpha_pow(ch) = double(alpha_broad_power);
                    iaf(ch) = freq;
                    slp(ch) = slope;
                    off(ch) = offset;
                    
                    clearvars any_alp freq range alpha_broad_power pow ifreq width detrend L logL bounds aperiodic offset slope
                    
                elseif sum(any_alp)>1 % if more than one peak in alpha range is present
                    
                    [val alp_loc] = max(pow(any_alp)); %out of two alpha peaks which one is higher in amplitude
                    loc = find(any_alp); % find two alpha peak locations out of every peak
                    alp_loc = loc(alp_loc); %higher peak location out of every peak
                    any_alp(:) = 0; % set everything to 0 as if no alpha peak was found
                    any_alp(alp_loc) = 1; %set the maxima location of alpha peak to 1
                    
                    freq = ifreq(any_alp);
                    range = bounds(any_alp,:);
                    
                    if diff(range) > 6      %whenever the width of the peak is more than 3Hz, set it to 3
                        range = [freq-3, freq+3];
                    end
                    
                    alpha_broad_power = sum(detrend(fs>=round(range(1),2) & fs<=round(range(2),2)))*0.25;
                    
                    alpha_pow(ch) = double(alpha_broad_power);
                    iaf(ch) = freq;
                    slp(ch) = slope;
                    off(ch) = offset;
                    
                    clearvars val alp_loc loc any_alp freq range alpha_broad_power pow ifreq width detrend L logL bounds aperiodic offset slope
                    
                else
                    error('no peak was found')
                    alpha_pow(ch) =NaN;
                    iaf(ch) = NaN;
                    slp(ch) = slope;
                    off(ch) = offset;
                end
            else 
                alpha_pow(ch) =NaN;
                iaf(ch) = NaN;
                slp(ch) = slope;
                off(ch) = offset;
            end
        else
            alpha_pow(ch) =NaN;
            iaf(ch) = NaN;
            slp(ch) = NaN;
            off(ch) = NaN;
            clearvars rsq
            
        end
    end
    
    alpha_params(p).ID = id
    alpha_params(p).power = alpha_pow
    alpha_params(p).freq = iaf
    alpha_params(p).slope= slp
    alpha_params(p).offset= off
    
    clearvars -except p f_orig fs files saveDir dataDir psd_all alpha_params
end

save([saveDir, 'alpha_params_findpeaks_detrend_1_40_7_13_threshold_05'], 'alpha_params')

%% create matrixes for alpha power %%

alpha_broad=struct()
alpha_broad.ID = {alpha_params.ID};

for i = 1:length(alpha_params)
  broadpow(i,:) = alpha_params(i).power;
  iaf(i,:) = alpha_params(i).freq;
  slope(i,:) = alpha_params(i).slope;
  offset(i,:) = alpha_params(i).offset;
end

save([saveDir, 'alpha_pow_mtx_05'], 'broadpow')
save([saveDir, 'alpha_freq_mtx_05'], 'iaf')
save([saveDir, 'slope_mtx_05'], 'slope')
save([saveDir, 'offset_mtx_05'], 'offset')

%% Negative values could appeared due to a small alpha peak that is estimated under the PSD curve. Check if there are any
[row channel] = find(broadpow<0)%NONE
