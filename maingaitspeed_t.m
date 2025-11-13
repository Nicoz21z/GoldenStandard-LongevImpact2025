%% 
clc;
clear all;
close all;
%% ---------------- LOADING DATA ----------------
load('Turns.mat');
fields = fieldnames(dati_hackathon);

%% ---------------- SUBJECT PARAMETERS ----------------
h  = 1.64;   sensor_h  = 0.57*h;   % Controller (sensor high)
h1 = 1.61;   sensor_h1 = 0.57*h1;  % Faller
%% ---------------- FILTERS PARAMETERS ----------------
fs = dati_hackathon.co001_base.Fs; 
% Cut-off frequency
fc = 15; 
% 4nd order Butterworth filter
[b, a] = butter(4, fc/(fs/2)); 


features_all = [];  % final matrix
labels_all = [];    % vector labels
%% ========================================================================
%                          NON FALLERS 
% ========================================================================

for i = 1:numel(fields)-35

    Name = fields{i};


    Acc_V_filt = filtfilt(b, a, dati_hackathon.(Name).acc_V);
    Acc_AP_filt = filtfilt(b, a, dati_hackathon.(Name).acc_AP);
    Acc_ML_filt = filtfilt(b, a, dati_hackathon.(Name).acc_ML);
    [Acc_APcorr, Acc_MLcorr, Acc_Vcorr] = algo_Moe_Nilssen(Acc_AP_filt, Acc_ML_filt, Acc_V_filt, 'tiltAndNoG');

% INTEGRATING THE ACC V:

    int_acc_V = cumtrapz(1/fs,Acc_Vcorr);

% DIFFERENTIATION THE SIGNAL OBTAINED BY USING A CWT:

    S1 = -cwt(int_acc_V,10,'gaus1',1/fs);

% USING FINDPEAKS TO IDENTIFY THE MINIMA OF S1:

    %[ICs, index1] = findpeaks(-S1,'MinPeakHeight',0.10);
    [ICs, index1] = findpeaks(-S1);
    IC = index1;
    IC_t = IC./fs;

% DIFFERENTIATING S1

    S2 = -cwt(S1,10,'gaus1',1/fs);

% USING FINDPEAKS TO IDENTIFY THE MAXIMA OF S2:

    [FCs, index2] = findpeaks(S2,'MinPeakDistance',50);
    FC = index2;
    FC_t = FC./fs;


    try
        %if an error in this function occurs we skip to the next subject
        [StepTime, StanceTime, StrideTime, SwingTime, StepLength, StepVelocity]=calculateSpatioTemporalGaitCharacteristics(IC,FC,Acc_Vcorr,sensor_h,fs);

    catch ME

        warning('Errore nel soggetto %d: %s. Passo al successivo.', i, ME.message);

        % plot the signal

        % figure()
        % plot(Acc_Vcorr)
        % title(sprintf('Subject %d . Select start and stop', i))
        % xlabel('Sammples')
        % ylabel('Acc V corr')

            continue 

    end
    
    gait_speed(i) = mean(StepVelocity); %gait speed of the controller subjects 

    % mean of features 
feat = [ ...
    mean(StepTime, 'omitnan'), ...
    mean(StanceTime, 'omitnan'), ...
    mean(StrideTime, 'omitnan'), ...
    mean(SwingTime, 'omitnan'), ... 
    mean(StepLength, 'omitnan'), ...
    mean(StepVelocity, 'omitnan')  % gait speed = mean velocity
];

features_all = [features_all; feat];
labels_all = [labels_all; 0];  % 0 = control
    
end


%% ========================================================================
%                           FALLERS 
% ========================================================================

for i = numel(fields)-35:numel(fields)

    Name = fields{i};

% filtering the segmented signals
    Acc_V_filt = dati_hackathon.(Name).acc_V;
    Acc_AP_filt = dati_hackathon.(Name).acc_AP;
    Acc_ML_filt =  dati_hackathon.(Name).acc_ML;

    [Acc_APcorr, Acc_MLcorr, Acc_Vcorr] = algo_Moe_Nilssen(Acc_AP_filt, Acc_ML_filt, Acc_V_filt, 'tiltAndNoG');

% INTEGRATING THE ACC V:

    int_acc_V = cumtrapz(1/fs,Acc_Vcorr);

% DIFFERENTIATION THE SIGNAL OBTAINED BY USING A CWT:

    S1 = -cwt(int_acc_V,10,'gaus1',1/fs);

% USING FINDPEAKS TO IDENTIFY THE MINIMA OF S1:

    [ICs, index1] = findpeaks(-S1,'MinPeakHeight',0.10);
    IC = index1;
    IC_t = IC./fs;

% DIFFERENTIATING S1

    S2 = -cwt(S1,10,'gaus1',1/fs);

% USING FINDPEAKS TO IDENTIFY THE MAXIMA OF S2:

    [FCs, index2] = findpeaks(S2,'MinPeakDistance',200);
    FC = index2;
    FC_t = FC./fs;


    try
        % if an error in this function occurs we skip to the next subject
        [StepTime1, StanceTime1, StrideTime1, SwingTime1, StepLength1, StepVelocity1]=calculateSpatioTemporalGaitCharacteristics(IC,FC,Acc_Vcorr,sensor_h1,fs);

    catch ME

        warning('Errore nel soggetto %d: %s. Passo al successivo.', i, ME.message);

  
        continue
    end

    
    gait_speed_1(i) = mean(StepVelocity1); %gait speed of the faller subjects

    feat = [ ...
    mean(StepTime1, 'omitnan'), ...
    mean(StanceTime1, 'omitnan'), ...
    mean(StrideTime1, 'omitnan'), ...
    mean(SwingTime1, 'omitnan'), ...
    mean(StepLength1, 'omitnan'), ...
    mean(StepVelocity1, 'omitnan')  % gait speed
];

features_all = [features_all; feat];
labels_all = [labels_all; 1];  % 1 = faller


end
% We remove the zeros in the giat speed vector due to the possible errors in the previous function
gait_speed = gait_speed(gait_speed~=0);
gait_speed_1 = gait_speed_1(gait_speed_1~=0);


%we calculate mean and standard deviation
average_speed = mean(gait_speed);
std_gait_speed = std(gait_speed);


average_speed_1 = mean(gait_speed_1);
std_gait_speed_1 = std(gait_speed_1);

% boxplots
figure
subplot(121)
boxplot(gait_speed);
title('Controllers gait speed [m/s]');
subplot(122)
boxplot(gait_speed_1);
title('Fallers gait speed [m/s]');

%no outliers
gait_speed_1na = [gait_speed_1(1:14) gait_speed_1(16:29)];
average_speed_1na = mean(gait_speed_1na);
std_gait_speed_1na = std(gait_speed_1na);
%%
figure
subplot(121)
boxplot(gait_speed);
title('Controllers gait speed [m/s]');
subplot(122)
boxplot(gait_speed_1na);
title('Fallers gait speed, no outliers [m/s]');

% Return gate speed in the command window:

fprintf('Controllers average speed = %.2f ± %.2f m/s\n', average_speed, std_gait_speed);
fprintf('Fallers average speed = %.2f ± %.2f m/s\n', average_speed_1, std_gait_speed_1);
fprintf('Fallers average speed without outliers = %.2f ± %.2f m/s\n', average_speed_1na, std_gait_speed_1na);

%% =============================== ML MODEL =========================================
feature_names = {'StepTime','StanceTime','StrideTime','SwingTime','StepLength','GaitSpeed'};
T = array2table(features_all, 'VariableNames', feature_names);
T.Label = labels_all;  % aggiungi la colonna label

% label in the first column 
T = movevars(T, 'Label', 'Before', 1);
writetable(T, 'GaitFeaturesTable.csv');

% Index controller (label = 0)
idx_ctrl = find(T.Label == 0);

% Index faller (label = 1)
idx_fall = find(T.Label == 1);

% Training: first 30 controller and first 23 faller
train_idx = [idx_ctrl(1:30); idx_fall(1:23)];

% Testing: last 8 controller e last 6 faller 
test_idx = [idx_ctrl(end-7:end); idx_fall(end-5:end)];

% Table for training and testing the ML classification model 
T_train = T(train_idx, :);
T_test  = T(test_idx, :);

% save in a file 
writetable(T_train, 'GaitFeaturesTable_train.csv');
writetable(T_test,  'GaitFeaturesTable_test.csv');
