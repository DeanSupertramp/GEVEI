function [eta_ICE]=mappa_eta_ICE(w,torque,varargin)
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
    torque=torque.^(1.25^(5-1));
end
else
end
eta_ICE=-(0.8e-8)*(w-3630/60*2*pi).^2 - (3.2e-6).*(torque-160).^2+0.41;
