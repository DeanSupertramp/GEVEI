close all
clear 

% Costanti
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
ro = 1.225; % [Kg/m^3] Densità dell'aria
dt = 0.1;
alpha=0;

w_max_m = 2*pi*giri_max/60; % [giri/sec] motore
v_max_r = w_max_m*r/G; % velocità max ruota
vc = wc*r/G; % velocità critica
Pmax = Tmax*wc;
Ftr0=Tmax*G*eff/r; % forza di trazione, supponiamo sempre massima

%definisco e inizializzo i vettori 
t = 0:dt:30;

v = zeros(1,length(t));

a = zeros(1,length(t));
a(1) = (Ftr0-(mu_rr*m*g)-(0.5*A*Cd*ro*v(1)^2))/m_eff;

Ftr = ones(1,length(t)).*Ftr0;

x = zeros(1,length(t));

P = zeros(1,length(t));
P(1)=Ftr(1)*v(1);

L = zeros(1,length(t));

% loop
for i = 2:length(t)
    v(i) = v(i-1)+dt*a(i-1);
    if v(i)>= v_max_r
        v(i)=v_max_r;
    end
    x(i) = x(i-1) + dt*(v(i)+v(i-1))/2; %formula del trapezio (integrale)
    if (v(i)<=vc)
        a(i) = (Ftr0-(mu_rr*m*g)-(0.5*A*Cd*ro*v(i)^2))/m_eff;
    else
        Ftr(i)=eff*Pmax/v(i);
        a(i) = (Ftr(i)-(mu_rr*m*g)-(0.5*A*Cd*ro*v(i)^2))/m_eff;
    end
    P(i)=Ftr(i)*v(i);
end

L=cumsum(P).*dt; %lavoro = potenza * tempo

Ek = 0.5*m_eff*v.^2; % energia cinetica

Rr = mu_rr*m*g; % resistenza al rotolamento
Pr = v.*Rr; % potenza associata alla resistenza al rotolamento
Lr = cumsum(Pr).*dt; % lavoro associato alla resistenza al rotolamento

Ra = 0.5*ro*A*Cd*v.^2;
Pa = v.*Ra;
La = cumsum(Pa).*dt;

figure(1)
subplot(1,2,1);
plot(t,v)
title("Velocità")
xlabel("Tempo [s]")
ylabel("Velocità [m/s]")
grid on;
%figure(2)
subplot(1,2,2);
plot(t,a)
title("Accelerazione")
xlabel("Tempo [s]")
ylabel("Accelerazione [m/s^2]")
grid on;

figure(3)
plot(t,Ek,t,L,t,Lr,t,La)
title("Energia cinetica, Lavoro, Lavoro rotolamento, Lavoro aereodinamico")
xlabel("Tempo [s]")
ylabel("Ek, L, Lr, La")

figure(4)
plot(t,x)
title("Spazio percorso")
xlabel("Tempo [s]")
ylabel("Distanza [m]")
grid on;

