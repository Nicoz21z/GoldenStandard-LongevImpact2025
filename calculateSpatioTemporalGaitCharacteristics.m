function [StepTime, StanceTime, StrideTime, SwingTime, StepLength, StepVelocity]=calculateSpatioTemporalGaitCharacteristics(IC,FC,accV,w_height,Fs)


% everything should start with an IC
if FC(1)<IC(1)
    FC(1)=[];
end

% n_IC = n_FC
min_len=min([length(FC);length(IC)]);
if length(FC)~=length(IC)
    FC=FC(1:min_len);
    IC=IC(1:min_len);
end

%check IC->FC->IC
for i=1:min_len
    diff=FC(i)-IC(i);
    if diff<0
        FC(i)=[];
        if i<min_len && IC(i+1)-IC(i)>20 
            IC(i)=[];
        end
    end
end


% n_IC = n_FC
if length(FC)~=length(IC)
    min_len=min([length(FC);length(IC)]);
    FC=FC(1:min_len);
    IC=IC(1:min_len);
end

% Plot again
figure
plot(accV,'b');
hold on;
plot(IC,accV(IC),'r*');
plot(FC,accV(FC),'g*');
xlabel('samples')
legend('acc','IC','FC');

for i=1:length(IC)-1
    % Calculating temporal gait characteristics from the IC/FC events
    StepTime(i)=(IC(i+1)-IC(i))/Fs;
    if i+1<=length(FC)%In case the last event present is an IC
        StanceTime(i)=(FC(i+1)-IC(i))/Fs;
    else
        StanceTime(i)=NaN;
    end
    if i+2>length(IC)-1 %ICs are finished so you cannot calculate stride (but you still can calculate one last step)
        StrideTime(i)=NaN;
    else
        StrideTime(i)=(IC(i+2)-IC(i))/Fs;
    end
    
    % Integration and spatio-temporal estimations
      StepAcc=detrend(accV(IC(i):IC(i+1)));
      hvel=cumtrapz(1/Fs,StepAcc); %you must use the sampling period to obtain the velocity and then the displacement
      h=cumtrapz(1/Fs, hvel);
  
      h=max(h)-min(h);
    
    StepLength(i)=2*sqrt(2*(w_height)*h-h^2);
    StepVelocity(i)=StepLength(i)/StepTime(i);
end
     
SwingTime = StrideTime-StanceTime;
