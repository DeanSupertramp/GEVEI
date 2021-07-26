%%
clc
close all
clear all

%% Caricamento del ciclo di guida standard
load('WLTC_full.mat');
load('FUDS.mat');
t=FUDS(:,1)';
v=FUDS(:,2)';
%t=[0:1:(length(WLTC_full(:,1))+length(FUDS(:,1))-1)]';
%v=[WLTC_full(:,2);FUDS(:,2)]';

Kc=0.40;    %copper losses
Ki=8e-3;    %iron losses
Kw=0.2e-5;  %windage losses
C=300;      %constant losses

%% Parametri del veicolo 
Cd=0.26;                    %coefficiente aereodinamico
mu_r0=0.024;                %rolling resistance
G = 1;                      %rapporto %gtr(itr)
G0=4.113;
r = 0.32;                   %raggio ruota [m]
m = 1370;                   % kg con 2 passeggeri da 70kg l'uno
m_eff = 1421;               % massa efficace considerando il momento di inerzia del motore
wc_el = 1300/(60/(2*pi));   %rad s^-1, w critico motore elettrico
Tmax_el = 65;               %N*m coppia massima motore elettrico
Pmax_el = wc_el*Tmax_el;    %W

vc_ic_rpm=3300;             %rad s^-1, w critico motore a combustione
wc_ic=vc_ic_rpm/(60/(2*pi));
Tmax_ic=155;                % N*m coppia massima motore a combustione
Pmax_ic=wc_ic*Tmax_ic;

A = 2.33;                   % area, m^2
eff_mech = 0.96;            % efficienza della trasmissione eta meccanico
v_max_r = 225/3.6;          % m/s velocità massima veicolo
g = 9.81;                   % m/s^2
ro = 1.225;                 % [Kg/m^3] Densità dell'aria
dt = t(2)-t(1);
alpha=0;                    %angolo di pendenza [gradi]

Q0=6500*3.6;                % [C] capacità iniziale pacco batteria
V0=248;                     % tensione iniziale pacco batteria
E0=5.8e6;                   % [J] Energia iniziale pacco batteria
R0=0.12;                    % [Ohm] Resistenza interna batteria
SOC0=0.85;                  % SOC iniziale
SOC_target=0.75;
m1=33;                      % pendenza della curva di scarica linearizzata
k=1.01;                     % coefficiente di Peukert
PbattMax=31000;             % [W]
IbattMax=130;               % [A]


%% Inizializzazione vettori
Ftr=zeros(1,length(t));             %forza di trazione
Ftr_p=zeros(1,length(t));           %forza di trazione positiva
Ftr_n=zeros(1,length(t));           %forza di trazione negativa
Ltr_p=zeros(1,length(t));           %lavoro di trazione positivo
Ltr_n=zeros(1,length(t));           %lavoro di trazione negativo
P=zeros(1,length(t));               %Potenza meccanica sulle ruote
P_p=zeros(1,length(t));             %Potenza meccanica positiva
P_n=zeros(1,length(t));             %Potenza meccanica negativa
P_motore=zeros(1,length(t));        %Potenza meccanica del motore
P_eff_p=zeros(1,length(t));         %Potenza meccanica trasmessa dal motore
P_eff_n=zeros(1,length(t));         %Potenza meccanica trasmessa al motore
T=zeros(1,length(t));               %Coppia
T_p=zeros(1,length(t));             %Coppia positiva
T_n=zeros(1,length(t));             %Coppia negativa
P_batt=zeros(1,length(t));          %Potenza richiesta alla batteria
P_batt_p=zeros(1,length(t));        %Potenza erogata dalla batteria
P_batt_n=zeros(1,length(t));        %Potenza fornita alla batteria
eff_motore_el=zeros(1,length(t));   %Efficienza del motore elettrico
eff_gen=zeros(1,length(t));         %Efficienza del generatore
eff_motore_ic=zeros(1,length(t));
w=zeros(1,length(t));               %velocità angolare di rotazione del motore
I_p=zeros(1,length(t));             %corrente erogata dalla batteria
I_n=zeros(1,length(t));             %corrente fornita alla batteria
Q=zeros(1,length(t));               %capacità della batteria
SOC=zeros(1,length(t));             %SOC della batteria
D_SOC_p=zeros(1,length(t));         %D_SOC della batteria
D_SOC_n=zeros(1,length(t));         %D_SOC della batteria
V=zeros(1,length(t));               %tensione della batteria

Tel_0=0;
Tel_n=ones(1,length(t))*Tel_0;
Tel_p=ones(1,length(t))*Tel_0;
P_ic=zeros(1,length(t));
%Calcolo vettore accelerazione, con a(1)=0
a=[0,diff(v)/dt]; 

%% Calcolo valori all'istante t=0
Ftr_max=Tmax_el*G*eff_mech/r;
vc = wc_el*r/G;
Rg=m*g*sind(alpha);                 %forza di arrampicata in salita

x(1)=0;
Ftr(1)=mu_r0*m*g*cosd(alpha)+Rg+(0.5*A*Cd*ro*(v(1)^2)+(m_eff*a(1)));
SOC(1)=SOC0;
V(1)=V0;
eff_motore_el(1)=0;
eff_gen(1)=0;
eff_tot(1)=0;

P_p(1)=Ftr(1)*v(1);
P_eff_p(1)=P_p(1)/eff_mech;
w(1)=v(1)/r*G0;
T_p(1)=0;
Tic_0=(T_p(1)/G-Tel_p(1));
Tic=ones(1,length(t))*Tic_0;

P_pk_p(1)=P_eff_p(1); 
P_pk_n(1)=0;

L1_0_vec=[0:1:7];               %L1_0_vec=[0:0.5:5];
L1=ones(1,length(t));

%% SIMULAZIONE
for j=1:length(L1_0_vec)        % Shooting Method sui valori iniziali di ?
    L1_0=L1_0_vec(j);
    L1(1)=L1_0;
    L1(2)=L1_0;
    for i = 2:(length(t))       % calcolo traiettorie temporali, evoluzione nel tempo del sistema

        x(i) = x(i-1) + dt*(v(i)+v(i-1))/2;
        mu_r=mu_r0*(1+v(i)/100);
        Rr0=mu_r*m*g;
        Rr=Rr0*cosd(alpha);     % resistenza al rotolamento
        Ftr(i)=Rr+Rg+(0.5*A*Cd*ro*(v(i)^2)+(m_eff*a(i)));   % forza di trazione

        Ftr_p(i)=max(0,Ftr(i));
        Ftr_n(i)=abs(min(0,Ftr(i)));

        P_p(i)=Ftr_p(i)*v(i);
        P_n(i)=Ftr_n(i)*v(i);

        P_eff_p(i)=P_p(i)/eff_mech;         % P_eff=potenza effettiva da dividere tra IC e EL
        P_eff_n(i)=P_n(i)*eff_mech;

        w(i)=v(i)/r*G0;                     % w motore
        v_rpm=(60/(2*pi))*w(i);

        if w(i)==0
            T_p(i)=0;
            T_n(i)=0;
            Tmax_el_w=0;
        else
            T_p(i)=P_eff_p(i)/w(i);
            T_n(i)=P_eff_n(i)/w(i);
            if(w(i)>wc_el)
                Tmax_el_w=Pmax_el/w(i);
            else
                Tmax_el_w=Tmax_el;
            end
        end    
%% Calcolo valore ottimo di coppia del motore elettrico al tempo t attraverso la minim. dell'Hamiltoniana
        if(Ftr(i)>0)    % caso trazione
           if((T_p(i)/(G*G0))>Tmax_el_w)
               Tstop=Tmax_el_w;
           else
               Tstop=T_p(i)/(G*G0);
           end
           T_opt=linspace(0,Tstop,1000);
           T_ice=T_p(i)/(G*G0)-T_opt;     
           eff_motore_ic_vec=mappa_eta_ICE(w(i),T_ice,E0);
           Pfuel_vec=(T_ice*w(i))./eff_motore_ic_vec;
           eff_motore_el_vec=mappa_eta_EL_mot(w(i),T_opt,E0);
           P_batt_p_vec=T_opt*w(i)./eff_motore_el_vec;
           I_p_vec=P_batt_p_vec./V(i-1);
           D_SOC_p_vec=((I_p_vec).^k)*dt/Q0;
           SOC_vec=SOC(i-1)-D_SOC_p_vec;
           V_vec=215+m1.*(SOC_vec);
           P_pk_p_vec=I_p_vec.^k .* V_vec;
           SOCp_vec=-(1/Q0)*(V_vec/(2*R0)-sqrt((V_vec/(2*R0)).^2-((P_pk_p_vec)/R0)));
           Pech_vec=-V_vec*Q0.*SOCp_vec;
           w_SOC=5*(0.9-SOC(i-1));
           H_vec=Pfuel_vec+(L1(i)+w_SOC).*Pech_vec;          
           [Hmin,min_index]=min(H_vec);
           Tel_p(i)=T_opt(min_index);
           Tel_n(i)=0;
        else        % caso frenata rigenerativa
           if(T_n(i)>Tmax_el_w)
               Tstop=-Tmax_el_w;
           else
               Tstop=-T_n(i);
           end
           T_opt=linspace(Tstop,0,1000);
           eff_motore_el_vec=mappa_eta_EL_gen(w(i),T_opt,E0);
           P_batt_n_vec=T_opt*w(i)./eff_motore_el_vec;
           I_n_vec=-P_batt_n_vec./V(i-1);
           D_SOC_n_vec=((I_n_vec).^(1/k))*dt/Q0;
           SOC_vec=SOC(i-1)+D_SOC_n_vec; 
           SOC_vec(find(SOC_vec>1))=1;
           V_vec=215+m1.*(SOC_vec);
           P_pk_n_vec=I_n_vec.^(1/k) .* V_vec;
           SOCp_vec=-(1/Q0)*(V_vec/(2*R0)-sqrt((V_vec/(2*R0)).^2-((-P_pk_n_vec)/R0)));
           Pech_vec=-V_vec*Q0.*SOCp_vec;
           w_SOC=5*(0.9-SOC(i-1));
           H_vec=(L1(i)+w_SOC).*Pech_vec;
           [Hmin,min_index]=min(H_vec);
           Tel_p(i)=0;
           Tel_n(i)=-T_opt(min_index);
        end
%%        
        if w(i)==0
            T_p(i)=0;
            T_n(i)=0;
            Tic(i)=0;
            Tel_p(i)=0;
            Tel_n(i)=0;
            eff_motore_el(i)=0;
            eff_gen(i)=0;
        else
            Tic(i)=(T_p(i)/(G*G0))-Tel_p(i);    
            eff_motore_el(i)=mappa_eta_EL_mot(w(i),Tel_p(i),E0);
            eff_gen(i)=mappa_eta_EL_gen(w(i),Tel_n(i),E0);
            eff_motore_ic(i)=mappa_eta_ICE(w(i),Tic(i),E0);
            P_ic(i)=Tic(i)*w(i)/eff_motore_ic(i);
        end
%% Aggiornamento parametri sistema (SoC e lambda)        
        P_batt_p(i)=Tel_p(i)*w(i)/eff_motore_el(i);
        if(P_batt_p(i)>PbattMax)
            P_batt_p(i)=PbattMax;
        end
        if(isnan(P_batt_p(i)))
            P_batt_p(i)=0;
        end
        P_batt_n(i)=Tel_n(i)*w(i)*eff_gen(i);   % calcolo la potenza fornita alla batteria
        if(isnan(P_batt_n(i)))
            P_batt_n(i)=0;
        end

        I_p(i)=P_batt_p(i)/V(i-1);              % corrente erogata dalla batteria
        if(I_p(i)>IbattMax)
            I_p(i)=IbattMax;
        end
        I_n(i)=P_batt_n(i)/V(i-1);              % corrente fornita alla batteria

        D_SOC_p(i)=((I_p(i))^k)*dt/Q0;
        D_SOC_n(i)=(I_n(i)^(1/k))*dt/Q0;        % Con frenata rigenerativa
        %D_SOC_n(i)=0;                          % senza frenata rigenerativa
        SOC(i)=SOC(i-1)-D_SOC_p(i)+D_SOC_n(i);     
        if(SOC(i)>1)
            SOC(i)=1;
        end
        Q(i)=Q0*SOC(i);                         % carica all'istante i
        V(i)=215+m1*(SOC(i));                   % tensione nominale all'istante i
        P_pk_n(i)=I_n(i).^(1/k) .* V(i);
        P_pk_p(i)=I_p(i).^(k) .* V(i);
        dSOCp(i)=((5e-4)/(35e3))*((P_pk_p(i)-P_pk_n(i))/1e3);
        L1(i+1)=L1(i)-dt*L1(i)*dSOCp(i);

    end 
% ciclo temporale concluso, si calcola la funzione costo
% e lo scostamento tra SoC finale e quello target
    SOC_end(j)=SOC(end);
    f_SOC_target(j)=SOC_end(j)-SOC_target;
    J(j)=dt*sum(P_ic);
    figure(1);
    plot(t,SOC*100,'color',[0, 0.4470, 0.7410]);
    hold on;
end

figure(1);
plot0=plot(t,ones(1,length(t))*SOC_target*100, '--k');
xlabel('Tempo [s]')
ylabel('SOC [%]')
title("State Of Charge")

%% %%%%%%%%%%%%%%BISEZIONE%%%%%%%%%%%%%%%%%%%%
bisezione2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

plot1=plot(t,SOC*100,'r'); 
legend([plot0 plot1],"SoC target","Soluzione ottimale", 'Location','southwest')

figure(2)
subplot(2,1,1)
plot(L1_0_vec,J)
ylabel('J')
title("J(\lambda_0)")
subplot(2,1,2)
plot(L1_0_vec,f_SOC_target, linspace(0,L0,10), zeros(1,10),'--k')
lim=get(gca,'YLim');
line([L0 L0],[lim(1),0],'Color','black','LineStyle','--')
ylabel('f\_SOC\_target')
xlabel("\lambda_0")
title("SOC(t_F)-SOC_{TARGET}")

figure()
subplot(3,1,1)
plot(t,v)
ylabel('v(t)')
title("Ciclo di guida - Velocità [m/s]")
subplot(3,1,2)
plot(t,a)
ylabel('a(t)')
title("Accelerazione [m/s^2]")
subplot(3,1,3)
plot(t,x/1000)
xlabel('Tempo [s]')
ylabel('x(t)')
title("Distanza percorsa [Km]")


%% MAPPA DI EFFICIENZA MOTORE A COMBUSTIONE
n_punti=100;
T_val=[1:Tmax_ic/n_punti:Tmax_ic]';
w_val=linspace(0,628,n_punti)';
[w_g,T_g]=meshgrid(w_val,T_val);        %griglia dei punti T,w
v_rpm1=(60/(2*pi)).*w_g;
eta_ic=-(0.8e-8)*(v_rpm1-3630).^2 - (3.2e-6)*(T_g-160).^2 +0.41;
eta_ic(eta_ic<0)=0.1;
max_p=[w_val,zeros(size(w_val)),zeros(size(w_val))]; % è una matrice con 3 colonne:
for jk=1:size(w_val)
    if w_val(jk)<wc_ic
        max_p(jk,2)=Tmax_ic;            % La coppia è massima fino ad w_c
    else
        max_p(jk,2)=Pmax_ic/w_val(jk);  % Superata w_c, la potenza è costante (pari a P_max) e la coppia è inversamente proporzionale ad w
    end
    max_p(jk,3)=-(0.8e-8)*(w_val(jk)/(60/(2*pi))-3630).^2 - (3.2e-6)*(max_p(jk,2)-160).^2 +0.41;
end
eta_values=[.25:.01:0.41];
figure()
p1=surface(w_g,T_g,T_g*0,eta_ic,'edgecolor','none');
set(gca,'view',[0 90])
colorbar
hold on
p2=plot3(max_p(:,1),max_p(:,2),max_p(:,3),'k');
p3=plot3(max_p(:,1),max_p(:,2).*eff_mech,max_p(:,3),'b');
plot3(w,Tic,eff_motore_ic,'ob')
contour3(w_g,T_g,eta_ic,eta_values,'r','ShowText','on');
legend([p1, p2, p3],'efficienza','curva limite ideale','curva limite reale')
ylabel('Coppia')
xlabel('\omega_{IC}', 'fontsize',16)
title("Mappa di efficienza del motore a combustione")


%% MAPPA DI EFFICIENZA MOTORE ELETTRICO
n_punti=100;
T_val=[1:Tmax_el/n_punti:Tmax_el]';
w_val=linspace(0,628,n_punti)';
[w_g,T_g]=meshgrid(w_val,T_val);                            % griglia dei punti T,w
OutPower=w_g.*T_g;                                          % potenza disponibile alla meccanica
TotPower=OutPower+w_g.^3*Kw+w_g*Ki+T_g.^2*Kc+C;             % potenza totale erogata
eta_motor=OutPower./TotPower;                               % efficienza del motore

max_p=[w_val,zeros(size(w_val)),zeros(size(w_val))];        % è una matrice con 3 colonne:
for jk=1:size(w_val)
    if w_val(jk)<wc_el
        max_p(jk,2)=Tmax_el;                                % La coppia è massima fino ad w_c
    else
        max_p(jk,2)=Pmax_el/w_val(jk);                      % Superata w_c, la potenza è costante (pari a P_max) e la coppia è inversamente proporzionale ad w
    end
    OP=w_val(jk).*max_p(jk,2);
    TP=OP+w_val(jk).^3*Kw+w_val(jk)*Ki+max_p(jk,2).^2*Kc+C;
    max_p(jk,3)=OP./TP;
end

eta_values=[.72:.04:1];
figure()
p1=surface(w_g,T_g,T_g*0,eta_motor,'edgecolor','none');
surface(w_g,-T_g,T_g*0,eta_motor,'edgecolor','none')
set(gca,'view',[0 90])
colorbar
hold on
plot3(max_p(:,1),max_p(:,2),max_p(:,3),'k')%
p3=plot3(max_p(:,1),max_p(:,2).*eff_mech,max_p(:,3),'b');
p2=plot3(max_p(:,1),-max_p(:,2),max_p(:,3),'k');
plot3(max_p(:,1),-max_p(:,2).*eff_mech,max_p(:,3),'b')
plot3(w,Tel_p,eff_motore_el,'ob')
plot3(w,-Tel_n,eff_gen,'ob')
contour3(w_g,T_g,eta_motor,eta_values,'r','ShowText','on');
contour3(w_g,-T_g,eta_motor,eta_values,'r','ShowText','on');
legend([p1, p2, p3],'efficienza','curva limite ideale','curva limite reale')
ylabel('Coppia')
xlabel('\omega', 'fontsize',16)
title("Mappa di efficienza del motore elettrico")



L_motore_p=cumsum(P_ic).*dt;
L_motore_el=cumsum(P_batt_p).*dt;
L_motore_el_n=cumsum(P_batt_n).*dt;
figure()
plot(t,L_motore_p./3.6e6,t,L_motore_el./3.6e6,t,L_motore_el_n./3.6e6)
xlabel('tempo [s]')
ylabel('Lavoro [kWh]')
title('Lavoro meccanico')
%%decommentare solo se SOC(0)=SOC_target
%hold on
%load('L_motore_IC.mat')
%plot(t, L_motore_IC/3.6e6)
%legend('Ibrido - motore IC','Ibrido - propulsione elettrica','Ibrido - frenata rigenerativa', 'Motore IC', 'Location','northwest')