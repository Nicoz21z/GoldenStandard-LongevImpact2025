clc,
clear all, 
close all

%% CONFIGURATION

%subfolder= 'Controllers3days\'
subfolder = 'Fallers3days\';
%Use DIR to find all data in the subfolder, automatically scanning
%both 'co' and 'fl' records
dataset = dir([subfolder '*.hea']);
%Initialize an empty structure to save all data
dati3days_hackathon = struct();

walk_tot=[];
step_tot=[];
median_bout_tot=[];
%-----------IMPORT DATA -----------------
% scelta = input('number pls: ');
for scelta=1:length(dataset)
    name_full = dataset(scelta).name;
    % Remove hea
    recordName = erase(name_full,'.hea');
    % Create complete path to pass to rdsamp ('dati_3days\CO001')
    fullPath = [subfolder recordName];
    % Check that the associated metadata file (.hea) exists
    heaFile = [fullPath '.hea'];

    if exist(heaFile, 'file')
        % If the .hea file exists, rdsamp reads and merges .dat and .hea.
        [segnale, Fs, tm] = rdsamp(fullPath);

        % Save the data in the 'hackathon_data' structure 
        dati3days_hackathon.(recordName).Segnale = segnale;
        dati3days_hackathon.(recordName).Fs = Fs;
        dati3days_hackathon.(recordName).Tempo = tm;
        disp(['Importato con successo: ' recordName]);
    else
        error(['File HEA per ' recordName ' non trovato.']);
    end

    %----------SELECT SIGNAL--------------
    acc_V_g = dati3days_hackathon.(recordName).Segnale(:,1);
    acc_AP_g = dati3days_hackathon.(recordName).Segnale(:,2);
    acc_ML_g = dati3days_hackathon.(recordName).Segnale(:,3);
    %convert the signal from g to m/s^2
    acc_V = acc_V_g.*9.81;
    acc_AP = acc_AP_g.*9.81;
    acc_ML = acc_ML_g.*9.81;
    
    %----------------FILTERING-------------------
    %acceleration filter with HP adn LP
    fc = 0.1; 
    [b,a] = butter(4, fc/(Fs/2), 'high');
    acc_V_f = filtfilt(b,a,acc_V);
    acc_AP_f = filtfilt(b,a,acc_AP);
    acc_ML_f = filtfilt(b,a,acc_ML);
    
    fc1 = 15;
    [b,a] = butter(4, fc1/(Fs/2), 'low');
    acc_V_filt = filtfilt(b,a,acc_V_f);
    acc_AP_filt = filtfilt(b,a,acc_AP_f);
    acc_ML_filt = filtfilt(b,a,acc_ML_f);
    
    %-------------MOE-NILSSEN CORRECTION----------
    %correct the tilt and the gravitational component
    [Acc_APcorr, Acc_MLcorr, Acc_Vcorr] = algo_Moe_Nilssen(acc_AP_filt, acc_ML_filt, acc_V_filt, 'tiltAndNoG');
    
    %------------NORMALIZATION-----------------------
    %normalization of the signal
    Acc_APcorr = (Acc_APcorr-mean(Acc_APcorr))/std(Acc_APcorr);
    Acc_Vcorr  = (Acc_Vcorr-mean(Acc_Vcorr))/std(Acc_Vcorr);
    Acc_MLcorr = (Acc_MLcorr-mean(Acc_MLcorr))/std(Acc_MLcorr);
    
    %----------------WINDOWING-------------------
    %----1. Segment the signal using a window of 1 sec with an overlap of 50%
    window_sec = 1;
    N = Fs*window_sec; 
    overlap = 0.5;
    step = round(N*(1-overlap));
    num_windows = floor((length(Acc_Vcorr) - N) / step) + 1;
    
    SMA_values = zeros(num_windows,1);
    energy_values = zeros(num_windows,1);
    ratio_values = zeros(num_windows,1);
    window_indices = zeros(num_windows,2);
    barLength=50;
    %----2. For each window we extract the locomotorian index(SMA, 
    for i = 1:num_windows
        idx_start = (i-1)*step+1;
        idx_end   = idx_start+N-1;
        accV_win  = Acc_Vcorr(idx_start:idx_end);
        accAP_win = Acc_APcorr(idx_start:idx_end);
        accML_win = Acc_MLcorr(idx_start:idx_end);
    
        % ----2.1 SMA
        SMA_values(i) = mean(abs(accV_win)+abs(accAP_win)+abs(accML_win));
    
        %-----2.2 Energy ratio (AP axis) with the PSD to estimate the peak
        %frequency
        [pxxAP,f] = pwelch(accAP_win,hanning(N),step,[],Fs);
        % integration in the band of 0.5-3 Hz 
        p_AP_locomotionband = bandpower(pxxAP,f,[0.5 3],'psd');
        durata = N/Fs;
        energy_values(i) = p_AP_locomotionband*durata;
    
        window_indices(i,:) = [idx_start, idx_end];
        perc = i/num_windows;
    
        nHash = round(perc*barLength);
    
        barStr = ['[' repmat('#',1,nHash) repmat(' ',1,barLength-nHash) ']'];
        
        % stampa sulla stessa riga
    
        fprintf('\r%s %3.0f%%', barStr, perc*100);
    end
    
    %-------- Setting the THRESHOLDS-------------
    SMA_threshold   = 0.135*9.81;  
    E_threshold     = 0.05*9.81;
    %--------- CHECK THE CONIDTION WITH AN OR---------
    cond_SMA = (SMA_values > SMA_threshold); 
    cond_E   = (energy_values > E_threshold) ;
    correct_mask = cond_SMA | cond_E;
    
    %------------LOCOMOTION MASK---------
    locomotion_mask = zeros(length(Acc_APcorr),1);
    for i=1:num_windows
        if correct_mask(i)
            locomotion_mask(window_indices(i,1):window_indices(i,2)) = 1;
        end
    end
    
    %----------FILTER BOUTS > 60s-----------
    %consider only the bouts>60 s
    minBoutSamples = 60*Fs;
    starts = find(diff([0;locomotion_mask])==1);
    stops  = find(diff([locomotion_mask;0])==-1);
    
    for i=1:length(starts)
        if (stops(i)-starts(i)) < minBoutSamples
            locomotion_mask(starts(i):stops(i)) = 0;
        end
    end
    
    %% --------------WALKING PERCENTAGE------------
    total_walk_samples = sum(locomotion_mask);
    perc_walk = total_walk_samples / length(Acc_APcorr) * 100;
    disp(perc_walk);
    %% ------------- STEP COUNT-------------------
    
    %---------1. Extract only the wlaking portion (mask=1)
    accAP_walk = Acc_APcorr(locomotion_mask==1);
    accV_walk  = Acc_Vcorr(locomotion_mask==1);
    
    %---------2. total duration of the walking 
    T_walking = length(accAP_walk)/Fs;   % in seconds
    
    %---------3. Calculate PSD on the AP axis (more informative for the frequency of the steps)
    [pxx,f] = pwelch(accAP_walk,hanning(4*Fs),[],[],Fs);
    
    %---------4. Consider the locomotion range 0.5–3 Hz
    valid_idx = (f>=0.5 & f<=3);
    [~,idx_dom] = max(pxx(valid_idx));
    f_step = f(valid_idx);
    f_step_dom = f_step(idx_dom);   % dominant frequency ≈ frequency of the steps
    
    %--------5. total number of steps
    N_steps = round(f_step_dom * T_walking);
    
    %% ----------- STRIDE DURATION --------------
    % all_stride = [];
    % for i = 1:length(starts)
    %     if locomotion_mask(starts(i))==0, continue; end
    %     accBout = Acc_APcorr(starts(i):stops(i));
    %     %-----1. autocorrelation shifts and compare the signla with itslef ad
    %     %different time lags--> walking is cyclic signla so has 
    %     % peaks at lags =  periodicity of the signal
    %     %periodicity of walking and dominant peaks are related to the stride
    %     %and step
    %     [r,lags] = xcorr(accBout,'coeff');
    %     lags_sec = lags/Fs;
    % 
    %     mask = (lags_sec > 0.5 & lags_sec < 2.0); % --> filter for realistic peak 
    %     r_pos   = r(mask);
    %     lags_pos = lags_sec(mask);
    % 
    %     if ~isempty(r_pos)
    %         [~,peakIdx] = max(r_pos);
    %         stride_duration = lags_pos(peakIdx);
    %         all_stride(end+1) = stride_duration;
    %     end
    % end
    %% WALKING BOUT DURATION (per replicare ~112 s del paper)
    minBoutSamples = 60*Fs;   % filtro: solo bout >= 60 s
    bout_durations = [];

    for i = 1:length(starts)
        dur_samples = stops(i) - starts(i) + 1;   % +1 per includere ultimo campione
        if dur_samples >= minBoutSamples
           bout_durations(end+1) = dur_samples / Fs; 
        end
    end

% Media e mediana
median_bout_duration = median(bout_durations); %less affected by the outliers

    
    walk_tot(scelta)= perc_walk;
    step_tot(scelta)=N_steps;
    median_bout_tot(scelta)=median_bout_duration;
   
end
