function [accH_AP, accH_ML, accV]=algo_Moe_Nilssen(accAP,accML,accV,tilt_or_tiltAndNoG)
%Acc in m/s^2
%two versions: 'tilt' only corrects orientation, 'tiltAndNoG' also clears g
%acceleration
%it provides as output the corrected accelerations
%NOTE. before using this, it should be verified that the vertical axis has a positive value near to +1g for
%still standing posture. Otherwise, if it is close to -1g, please first correct
%by multiplying by -1 the vertical values. 
%In this case, if possible, it should be understood how the sensor is flipped so flip appropriately also the AP or the ML axis.

g=9.81;

teta_ap = asin(mean(accAP)/g);
teta_ml = asin(mean(accML)/g);

accH_AP=(accAP.*cos(teta_ap))-(accV.*sin(teta_ap));
accV1=(accAP.*sin(teta_ap))+( accV.*cos(teta_ap));
accH_ML=(accML.*cos(teta_ml))-( accV1.*sin(teta_ml));

if strcmp(tilt_or_tiltAndNoG,'tilt')
    accV = (accML.*sin(teta_ml)) + (accV1.*cos(teta_ml));
elseif strcmp(tilt_or_tiltAndNoG,'tiltAndNoG')
  accV = (accML.*sin(teta_ml)) + (accV1.*cos(teta_ml))-g;  
else
    error('wrong input string');
end
