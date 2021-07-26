function [eta_EL]=mappa_eta_EL(w,torque,varargin)
if nargin==3
if w<140
    w=w/(1.25^(1-1));
    torque=torque.^(1.25^(1-1));
elseif 140<=w<140*(1+1.25^1)
    w=w/(1.25^(2-1));
    torque=torque.^(1.25^(2-1));
elseif 140*(1+1.25^1)<=w<140*(1+1.25^1+1.25^2)
    w=w/(1.25^(3-1));
    torque=torque.^(1.25^(3-1));
elseif 140*(1+1.25^1+1.25^2)<=w<140*(1+1.25^1+1.25^2+1.25^3)
    w=w/(1.25^(4-1));
    torque=torque.^(1.25^(4-1));
elseif w>140*(1+1.25^1+1.25^2+1.25^3)
    w=w/(1.25^(5-1));
    torque=torque.^(1.25^(5-1)); % Calcolo la coppia per diverse velocità
end
else
end

Kc=0.40;    %copper losses
Ki=8e-3;    %iron losses
Kw=0.2e-5;  %windage losses
C=300;      %constant losses
Creg=0;

nel=numel(torque);
WheelPower=w*abs(torque);
TotPowerMot=WheelPower+w.^3*Kw+w*Ki+torque.^2*Kc+C;
for jj=1:nel
        eta_EL(jj)=WheelPower(jj)/TotPowerMot(jj);
end


