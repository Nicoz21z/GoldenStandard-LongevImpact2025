function [ totalTurns,reportAllTurns,reportAverageAndStdTurns] = findTurns(gyroVertical, FS, thresholdDegrees, thresholdLPfiltering, thresholdPeakDS,thresholdCrossingDS,thresholdIntraTurnS,thresholdTurnDurationLowS, thresholdTurnDurationHighS )
% This function detectes all turns over a certain thresholds of degrees (thresholdDegrees)
%it assumes the vertical axis is aligned with the gravity
%gyro should be on the back (orr waist) sensor. 
%if no turn is detected then empty vectors are returned
totalTurns=[];
gyroVertical= lowPassFilter(gyroVertical,FS, thresholdLPfiltering);
if size(gyroVertical,2)>1
    error('should be a column vector')
end

dataAbs=abs(gyroVertical);

[PKS,LOCS]=findpeaks(dataAbs, 'MinPeakHeight',thresholdPeakDS);

%For each maximum (peak)
for i=1:size(LOCS,1)
    
    startTurn=[]; 
    endTurn=[];
    for j=1:LOCS(i)-1 %go backwards to find the crossing before the peak
        if dataAbs(LOCS(i)-j)<thresholdCrossingDS
            startTurn=LOCS(i)-j;
            break
        end
    end
    for j=1:size(dataAbs,1)-LOCS(i) %go forward to find the crossing after the peak
        if dataAbs(LOCS(i)+j)<thresholdCrossingDS
            endTurn=LOCS(i)+j;
            break
        end
    end
    directionTurn=sign(gyroVertical(LOCS(i))); %we need this for the following check on the direction of the turn
    if isempty(startTurn)||isempty(endTurn) %because you reached the beginning or end of the signal and you did not find those crossings
        %do nothing, do not add the turn
    else
        totalTurns=[totalTurns;[startTurn endTurn directionTurn ]];
    end
end
%since there can be more than one maximum between the same 5deg/s crossing,
%then use unique
totalTurns=unique(totalTurns,'rows');


timeAxis=(1:size(gyroVertical,1))/FS;
% figure;
% subplot(4,1,1)
% plot(timeAxis,gyroVertical,'g');
% xlabel('s')
% ylabel('deg/s')
% hold on;
% plot(timeAxis(LOCS),gyroVertical(LOCS),'o');
% for i=1:size(totalTurns,1) 
%     turnInterval=totalTurns(i,1):totalTurns(i,2);
% 
%     if totalTurns(i,3)==1 %left
%         plot(timeAxis(turnInterval),gyroVertical(turnInterval),'r');
%     else%right
%       plot(timeAxis(turnInterval),gyroVertical(turnInterval),'b');
%     end
% end

% title(['Step 1 - find maxima and 5 deg/s crossing, number of turns: ' int2str(size(totalTurns,1))]);
% legend('V gyro', 'maxima','right turn','left turn');


indexesToDelete=[];
%for each turn
for i=2:size(totalTurns,1)% from 2 because there must be a previous turn in order to do this combination
    intraTurnDurationS=(totalTurns(i,1)-totalTurns(i-1,2)+1)/FS; %(startTurn of this turn-endTurn of previous turn+1)/FS
    if intraTurnDurationS<thresholdIntraTurnS && totalTurns(i,3)==totalTurns(i-1,3) %intra turn duration less than threshold and previous turn is in the same direction
        totalTurns(i,:)=[totalTurns(i-1,1) totalTurns(i,2) totalTurns(i,3)];    
        indexesToDelete=[indexesToDelete;i-1]; %it does not exist anymore, it is combined in the next
    end
end
%now we need to eliminate those rows/turns
totalTurns(indexesToDelete,:)=[]; %delete



% subplot(4,1,2)
% plot(timeAxis,gyroVertical,'g');
% xlabel('s')
% ylabel('deg/s')
% hold on;
% for i=1:size(totalTurns,1) 
%     turnInterval=totalTurns(i,1):totalTurns(i,2);
%     if totalTurns(i,3)==1 %left
%     plot(timeAxis(turnInterval),gyroVertical(turnInterval),'r');
%     else
%       plot(timeAxis(turnInterval),gyroVertical(turnInterval),'b');
%     end
% end
% title(['Step 2 - combine near turns in the same direction, number of turns: ' int2str(size(totalTurns,1))]);


%for each remaining turn
indexesTurnsToDelete=[];
for i=1:size(totalTurns,1)
    turnDurationS=(totalTurns(i,2)-totalTurns(i,1))/FS;
    if turnDurationS<thresholdTurnDurationLowS || turnDurationS>thresholdTurnDurationHighS
        indexesTurnsToDelete=[indexesTurnsToDelete;i];
    end
end
totalTurns(indexesTurnsToDelete,:)=[];


% subplot(4,1,3)
% plot(timeAxis,gyroVertical,'g');
% xlabel('s')
% ylabel('deg/s')
% hold on;
% for i=1:size(totalTurns,1) 
%     turnInterval=totalTurns(i,1):totalTurns(i,2);
%    if totalTurns(i,3)==1 %left
%     plot(timeAxis(turnInterval),gyroVertical(turnInterval),'r');
%     else
%       plot(timeAxis(turnInterval),gyroVertical(turnInterval),'b');
%     end
% end
% title(['Step 3 - delete too short or too long turns, number of turns: ' int2str(size(totalTurns,1))]);



indexesTurnsToDelete=[];
degrees=zeros(size(totalTurns,1),1);
%for each remaining turn calculate turn angle
for i=1:size(totalTurns,1)
    degrees(i)=getTurningDegreesDS(gyroVertical(totalTurns(i,1):totalTurns(i,2)),FS);
    if abs(degrees(i))<thresholdDegrees
        indexesTurnsToDelete=[indexesTurnsToDelete;i];
   
    end
end
%add degrees of each turn
totalTurns=[totalTurns degrees];
%delete turns with too few degrees
totalTurns(indexesTurnsToDelete,:)=[];

% subplot(4,1,4)
% plot(timeAxis,gyroVertical,'g');
% xlabel('s')
% ylabel('deg/s')
% hold on;
% for i=1:size(totalTurns,1) 
%     turnInterval=totalTurns(i,1):totalTurns(i,2);
%     if totalTurns(i,3)==1 %left
%     plot(timeAxis(turnInterval),gyroVertical(turnInterval),'r');
%     else
%       plot(timeAxis(turnInterval),gyroVertical(turnInterval),'b');
%     end
% end
% title(['Step 4 - keep only turns with sufficient degrees, number of turns: ' int2str(size(totalTurns,1))]);


%linkaxes;%this is to synchronize the zoom on signals across the subplots
    

%now let's add a report with, for each turn, duration, turning angle (abs),
%and peak velocity

reportAllTurns=[];
for i=1:size(totalTurns,1) 
    reportAllTurns=[reportAllTurns; [(totalTurns(i,2)-totalTurns(i,1)+1)/FS totalTurns(i,4) max(dataAbs(totalTurns(i,1):totalTurns(i,2)))]];
   
end

reportAverageAndStdTurns=[mean(reportAllTurns);std(reportAllTurns)];

    