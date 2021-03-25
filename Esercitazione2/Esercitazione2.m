close all
clear all

% *********************** 1. Caricare ciclo di guida standard ***********************

load('cicliguida/FUDS.mat');
t = FUDS(:,1);
v = FUDS(:,2); %[m/s]
dt = t(2)-t(1);
plot(t,v);
title("Ciclo di guida standard FUDS");
xlabel("Tempo [s]")
ylabel("Velocit� istantanea [m/s]")

% --------------------- Costanti ---------------------
Cd=0.19; %coefficiente aereodinamico
mu_rr=0.0048; %rolling resistance
G = 11; %rapporto
r = 0.3; %[m] - raggio ruota 
m = 1540; % kg con 2 passeggeri da 70kg l'uno
m_eff = 1560; % [kg] - massa efficace considerando il momento di inerzia del motore
Tmax = 140; %[N*m] - coppia massima
wc = 733; %[rad s^-1] - w critico
A = 1.8; %[m^2] - area, 
eff = 0.95; % efficienza
giri_max = 12000; %[rpm] - giri al minuto
g = 9.8; %[m/s^2]
ro = 1.225; % [Kg/m^3] Densit� dell'aria
%dt = 0.1;
alpha=0; % [rad] pendenza
Rr0 = mu_rr*m*g;
Rr = Rr0 * cosd(alpha);
Rg = m*g*sind(alpha);
Ra0 = 0.5*ro*A*Cd;

% --------------------- Definisco e inizializzo i vettori ---------------------
x = zeros(1,length(t));
a = zeros(1,length(t));
F_m = zeros(1,length(t)); % [N]
F_m_pos = zeros(1,length(t)); % [N]
F_m_neg = zeros(1,length(t)); % [N]
P_m = zeros(1,length(t)); % [W] potenza motrice
P_m_pos = zeros(1,length(t)); % [W]potenza motrice
P_m_neg = zeros(1,length(t)); % [W]potenza motrice
L_m = zeros(1,length(t)); % [Wh]Lavoro della forza motrice
L_m_pos = zeros(1,length(t)); % [Wh] Lavoro della forza motrice
L_m_neg = zeros(1,length(t)); % [Wh] Lavoro della forza motrice
Ra_m = zeros(1,length(t)); % Resistenza aereodinamica
Rr_m = ones(1,length(t))*Rr;


% --------------------- Loop ---------------------
for i = 2:length(t)
% *********************** 2. Calcolare accelerazione e distanza percorsa  ***********************
    x(i) = x(i-1) + dt*(v(i)+v(i-1))/2; %formula del trapezio (integrale)
    a(i) = (v(i)-v(i-1))./dt;
    Ra_m(i) = Ra0*v(i-1).^2;
% *********************** 3. Calcolare forza di trazione istante per istante ***********************
    F_m(i-1) = Ra_m(i) + Rr + Rg + m_eff*a(i-1);
    F_m_pos(i-1) = max(F_m(i-1),0); % siccome il lavoro mi � utile per il calcolo dei consumi, uso solo la parte positiva del lavoro del motore
    F_m_neg(i-1) = min(F_m(i-1),0); 
% *********************** 4. Calcolare il lavoro della forza di trazione (positivo) per stimare i consumi  ***********************
    L_m(i) = L_m(i-1)+(x(i)+x(i-1)-x(i-1))*(F_m(i-1)+F_m(i))/2;
    L_m_pos(i) = L_m_pos(i-1)+(x(i)+x(i-1)-x(i-1))*(F_m_pos(i-1)+F_m_pos(i))/2;
% *********************** 5. Calcolare il lavoro della forza di trazione (negativo) per stimare energia recuperabile dalla frenata rigenerativa  ***********************
    L_m_neg(i) = L_m_neg(i-1)+(x(i)+x(i-1)-x(i-1))*(F_m_neg(i-1)+F_m_neg(i))/2;

    P_m(i) = F_m(i-1).*v(i-1);
    P_m_pos(i) = F_m_pos(i-1).*v(i-1);
    P_m_neg(i) = F_m_neg(i-1).*v(i-1);

end


% --------------------- Plot ---------------------
figure(2)
subplot(1,2,1);
plot(t,a)
title("Accelerazione")
xlabel("Tempo [s]")
ylabel("Accelerazione [m/s^2]")
subplot(1,2,2);
plot(t,x);
title("Distanza percorsa")
xlabel("Tempo [s]")
ylabel("Distanza [m]")

figure(3)
subplot(2,1,1);
plot(t,F_m)
title("Forza di trazione")
xlabel("Tempo [s]")
ylabel("Forza [N]")
subplot(2,1,2);
plot(t,F_m_pos)
title("Forza di trazione")
xlabel("Tempo [s]")
ylabel("Forza [N]")
hold on
plot(t,F_m_neg)


figure(4)
subplot(2,1,1);
plot(t,P_m/1000)
title("Potenza di trazione")
xlabel("Tempo [s]")
ylabel("Potenza")
subplot(2,1,2);
plot(t,P_m_pos/1000)
title("Potenza di trazione")
xlabel("Tempo [s]")
ylabel("Potenza [kW]")
hold on
plot(t,-P_m_neg/1000)
title("Potenza di trazione")
xlabel("Tempo [s]")
ylabel("Potenza [kW]")

figure(5)
subplot(2,1,1);
plot(t,L_m/1000)
title("Lavoro complessivo")
xlabel("Tempo [s]")
ylabel("Lavoro [kWh]")
subplot(2,1,2);
plot(t,L_m_pos/1000)
title("Lavoro")
xlabel("Tempo [s]")
ylabel("Lavoro [kWh]")
hold on
plot(t,-L_m_neg/1000)

%calcolo consumo del motore (lavoro motore), lo spazio percorso,
%teniamo conto dell'efficienza,
%per la forza di trazione, la moltiplichiamo per la massa efficace,
%sommiamo i termini resistivi, da cui poi ricaviamo potenza, lavoro etc
