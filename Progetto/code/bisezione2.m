%% Bisezione
% calcolo lambda ottimale data la diff tra SoC_finale e Soc_target
    f_p=find(f_SOC_target>0);
    L_sup=L1_0_vec(f_p(1));     % Valori di lambda inf e sup 
    L_inf=L1_0_vec(f_p(1)-1);   % vicini a quello ottimale del caso precedente 

    
    iterazioni=10;
for indice_bisez=1:iterazioni
    L1=ones(1,length(t))*L_sup;
    ciclo_tempo;
    SOC_fu=SOC(end)-SOC_target;
    
    L1=ones(1,length(t))*L_inf;
    ciclo_tempo;
    SOC_fl=SOC(end)-SOC_target;
    
    L_m=0.5*(L_sup+L_inf);  % lambda medio
    L1=ones(1,length(t))*L_m;
    ciclo_tempo;            % Ripeto il ciclo temporale per questi 3 valori di L
    SOC_fm=SOC(end)-SOC_target;
    
    fm=SOC_fu*SOC_fl;
    fu=SOC_fu*SOC_fm;
    fl=SOC_fm*SOC_fl;
    if fl>0         % Lower
       L_inf=L_m;
    elseif fu>0     % Upper
       L_sup=L_m;
    else 
       L0=L_m;      % Medio
       break
    end

end  

if(fm~=0)
    L0=0.5*(L_sup+L_inf);
    L1=ones(1,length(t))*L0;
    ciclo_tempo;    
    SOC_fm=SOC(end);   
end

[min_val, indice]=min([abs(SOC_fu-SOC_target),abs(SOC_fl-SOC_target),abs(SOC_fm-SOC_target)]);
if(indice==1)
    L0=L_sup;
elseif(indice==2)
    L0=L_inf;
else
    L0=L_m; 
end

L1=ones(1,length(t))*L0;
ciclo_tempo; 