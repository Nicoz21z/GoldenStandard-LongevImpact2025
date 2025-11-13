clc;
clear all;
close all;

%% CONFIGURATION
oldWarnState = warning; 
warning('off','all');
subfolder = 'LabWalks/'; 
%Use DIR to find all data in the subfolder, automatically scanning
%both 'co' and 'fl' records
dataset=dir([subfolder '*_base.hea']); 
%Initialize an empty structure to save all data
datiHackathon = struct(); 

%Loop to import into rdsamp
for i = 1:length(dataset)
    
    name_full=dataset(i).name;
    % Remove hea
    recordName=erase(name_full,'.hea');
    % Create complete path to pass to rdsamp ('LabWalks\c0001_base')
    fullPath =[subfolder recordName];
    % Check that the associated metadata file (.hea) exists
    heaFile = [fullPath '.hea'];
    
    if exist(heaFile, 'file')
        
        % If the .hea file exists, rdsamp reads and merges .dat and .hea.
        [segnale, Fs, tm] = rdsamp(fullPath);
        
        % Save the data in the 'hackathon_data' structure 
        dati_hackathon.(recordName).Segnale = segnale;
        dati_hackathon.(recordName).Fs = Fs;
        dati_hackathon.(recordName).Tempo = tm;
        
        disp(['Importato con successo: ' recordName]);
        
    else
        warning(['Saltato: File HEA per ' recordName ' non trovato.']);
    end
end 
save('Dati_Hackathon_Complessivi.mat', 'dati_hackathon');
%% ---------------- LOADING DATA ----------------
load('Dati_Hackathon_Complessivi.mat');
fields = fieldnames(dati_hackathon);

%% ---------------- SUBJECT PARAMETERS ----------------
h  = 1.64;   sensor_h  = 0.57*h;   % Controller (sensor high)
h1 = 1.61;   sensor_h1 = 0.57*h1;  % Faller
   %% ---------------- FILTERS PARAMETERS ----------------
%acceleration filter
fs = dati_hackathon.co001_base.Fs; 
fc_low = 0.5;
fc_high = 15; 
wn = [fc_low fc_high]/(fs/2);
order=4;
[b, a] = butter(order, wn, 'bandpass');         
%velocity filter 
fc_low_pos = 0.5;
fc_high_pos = 10; 
wn = [fc_low_pos fc_high_pos]/(fs/2);
order=4;
[b_pos, a_pos] = butter(order, wn, 'bandpass');  

%% ========================================================================
%                          NON FALLERS 
% ========================================================================

for i = 1:numel(fields)-35

    Name = fields{i};
    fprintf('Analyzing %s (controller)...\n', Name);

   % --- 1. extract 3 axis of acceleration ---
    Acc_V_segment  = dati_hackathon.(Name).Segnale(:,1);
    Acc_AP_segment = dati_hackathon.(Name).Segnale(:,2);
    Acc_ML_segment = dati_hackathon.(Name).Segnale(:,3);

   % --- 2. Attitude and gravity correction (Moe-Nilssen) ---
    [accH_AP, accH_ML, accV] = algo_Moe_Nilssen(Acc_AP_segment, Acc_ML_segment, Acc_V_segment, 'tiltAndNoG');
    g = 9.81;
    accH_AP=accH_AP*g;
    accH_ML=accH_ML*g;
    accV=accV*g;

  % --- 3. Calculate PSD to estimate the pitch frequency ---
    [Pxx, f] = pwelch(accV, [], [], [], fs);
    idx_range = (f >= 0.5 & f <= 3);
    [~, idx_peak] = max(Pxx(idx_range));
    f_band = f(idx_range);
    f_peak = f_band(idx_peak); % dominant frequency (passi/sec)

 % --- 4. acceleration filtering ---
    accV_final  = filtfilt(b, a, accV);
   
 % --- 5. Obtain vertical variation of the center of mass ---    
    vel = cumtrapz(1/fs, accV_final);  % integration -> velocity
    vel = filtfilt(b_pos, a_pos, vel);
    disp = cumtrapz(1/fs, vel);        % integration -> displacement
    hdisp = max(disp) - min(disp);     %vertical amplitude 

  % --- 6. Step Length (inverse pendulum) ---
    StepLength(i) = 2 * sqrt(2 * sensor_h * hdisp - hdisp^2);

  % --- 7. Gait speed ---
    gait_speed(i) = f_peak * StepLength(i);

    fprintf('f_step = %.2f Hz -> Step = %.2f m -> Speed = %.2f m/s\n', ...
        f_peak, StepLength(i), gait_speed(i));
end

%% ========================================================================
%                           FALLERS 
% ========================================================================
% filters parameters acc 
fs = dati_hackathon.co001_base.Fs; 
fc_low = 0.5;
fc_high = 15; 
wn = [fc_low fc_high]/(fs/2);
order=4;
[b, a] = butter(order, wn, 'bandpass');        
% velocity filter 
fc_low_pos = 0.5;
fc_high_pos = 10; 
wn = [fc_low_pos fc_high_pos]/(fs/2);
order=4;
[b_pos, a_pos] = butter(order, wn, 'bandpass');  

for i = numel(fields)-35:numel(fields)
 Name = fields{i};
    fprintf('Analyzing %s (controller)...\n', Name);

    % --- 1. extract 3 axis of acceleration ---
    Acc_V_segment  = dati_hackathon.(Name).Segnale(:,1);
    Acc_AP_segment = dati_hackathon.(Name).Segnale(:,2);
    Acc_ML_segment = dati_hackathon.(Name).Segnale(:,3);

   % --- 2. Attitude and gravity correction (Moe-Nilssen) ---
    [accH_AP, accH_ML, accV] = algo_Moe_Nilssen(Acc_AP_segment, Acc_ML_segment, Acc_V_segment, 'tiltAndNoG');
    g = 9.81;
    accH_AP=accH_AP*g;
    accH_ML=accH_ML*g;
    accV=accV*g;

    % --- 3. Calculate PSD to estimate the peak frequency ---
    [Pxx, f] = pwelch(accV, [], [], [], fs);
    idx_range = (f >= 0.5 & f <= 3);
    [~, idx_peak] = max(Pxx(idx_range));
    f_band = f(idx_range);
    f_peak = f_band(idx_peak); % dominant frequency (passi/sec)

    % --- 4. acceleration filtering ---
    accV_final  = filtfilt(b, a, accV);
   
    % --- 5. Obtain vertical variation of the center of mass ---  
    vel = cumtrapz(1/fs, accV_final);  % integration -> velocity
    vel = filtfilt(b_pos, a_pos, vel);
    disp = cumtrapz(1/fs, vel);        % integration -> displacement
    hdisp = max(disp) - min(disp);     % vertical amplitude 

    % --- 6. Step Length (inverse pendulum) ---
    StepLength1(i) = 2 * sqrt(2 * sensor_h1 * hdisp - hdisp^2);

    % --- 7. Gait speed ---
    gait_speed1(i) = f_peak * StepLength1(i);

    fprintf('f_step = %.2f Hz -> Step = %.2f m -> Speed = %.2f m/s\n', ...
        f_peak, StepLength1(i), gait_speed1(i));
end

gait_speed1 = gait_speed1(gait_speed1 ~= 0);
avg_speed=mean(gait_speed);
std_speed=std(gait_speed);

avg_speed1=mean(gait_speed1);
std_speed1=std(gait_speed1);

% boxplots
figure
subplot(121)
boxplot(gait_speed);
title('Controllers gait speed [m/s]');
subplot(122)
boxplot(gait_speed1);
title('Fallers gait speed [m/s]');

fprintf('\nControllers average speed = %.2f ± %.2f m/s\n', avg_speed, std_speed);
fprintf('\nFallers average speed = %.2f ± %.2f m/s\n', avg_speed1, std_speed1);

