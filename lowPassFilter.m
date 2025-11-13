function [data] = lowPassFilter(data,FS, frequencyForFilt)
%low pass filter of data (vector)
[C,D] = butter(1,frequencyForFilt/(FS/2));
data = filtfilt(C,D,data); 
end

