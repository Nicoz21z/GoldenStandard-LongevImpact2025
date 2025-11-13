%%
clc;
clear all;
close all;

%% CONFIGURATION
subfolder = 'LabWalks/';
%Use DIR to find all data in the subfolder, automatically scanning
%both 'co' and 'fl' records
dataset=dir([subfolder '*_base.hea']); 
%Initialize an empty structure to save all data
dati_hackathon = struct(); 

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
        [signal, Fs, tm] = rdsamp(fullPath);
        
        % Save the data in the 'hackathon_data' structure.
        dati_hackathon.(recordName).Signal = signal;
        dati_hackathon.(recordName).Fs = Fs;
        dati_hackathon.(recordName).Time = tm;
        
    else
        warning(['Skipped: File HEA per ' recordName ' not found.']);
    end
end 
save('Dati_Hackathon_Complessivi2.mat', 'dati_hackathon');

%% ========= cut all the subjects signals and save them in a new struct ===============

thresholdDegrees = 45; %only turns with degrees equal or over this threshold are identified
thresholdLPfiltering = 1.5; %Hz
thresholdPeakDS = 15; %degrees/second
thresholdCrossingDS = 5; %degrees/second
thresholdIntraTurnS = 0.05; %seconds
thresholdTurnDurationLowS = 0.5;
thresholdTurnDurationHighS = 10;


for i = 1:length(dataset)

    Record = dataset(i).name;
    Name = erase(Record, '.hea');

    acc_V = dati_hackathon.(Name).Signal(:,1);
    acc_AP = dati_hackathon.(Name).Signal(:,2);
    acc_ML = dati_hackathon.(Name).Signal(:,3);
    gyr_V = dati_hackathon.(Name).Signal(:,4);


    acc_V = acc_V.*9.81;
    acc_AP = acc_AP.*9.81;
    acc_ML = acc_ML.*9.81;

    t = dati_hackathon.(Name).Time;

    % figure
    % plot(gyr_V);
    % xlabel('Samples');
    % ylabel('gyr');

    [totalTurns,reportAllTurns,reportAverageAndStdTurns] = findTurns(gyr_V, Fs, thresholdDegrees, thresholdLPfiltering, thresholdPeakDS,thresholdCrossingDS,thresholdIntraTurnS,thresholdTurnDurationLowS, thresholdTurnDurationHighS );

    acc_V_Clean = acc_V; % copy of original signal
    acc_ML_Clean = acc_ML;
    acc_AP_Clean = acc_AP;

   
    mask = true(length(acc_V), 1);

    for k = 1:size(totalTurns,1)
        turnInterval = totalTurns(k,1):totalTurns(k,2);
        mask(turnInterval) = false;
    end

    acc_V_Clean = acc_V(mask);
    acc_ML_Clean = acc_ML(mask);
    acc_AP_Clean = acc_AP(mask);
    t_Clean = t(mask);


    dati_hackathon.(Name).acc_V = acc_V_Clean;
    dati_hackathon.(Name).acc_ML = acc_ML_Clean;
    dati_hackathon.(Name).acc_AP = acc_AP_Clean;
    dati_hackathon.(Name).T = t_Clean;

    figure

    subplot(2,1,1)
    plot(acc_V_Clean)
    title('cleaned signal')
    subplot(2,1,2)
    plot(acc_V)
    title('original signal')
    linkaxes


end


save('Turns.mat','dati_hackathon')