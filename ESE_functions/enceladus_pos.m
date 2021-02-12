% Enceladus coe (ref Saturn) , data ref 1/1/2028 alle 9:59:00
function [r, v, saturn_r] = enceladus_pos(year, month, day, hour, min, sec)

% pos saturno nella data stabilita
[saturn_coe, saturn_r, saturn_v, ~] =...
                        planet_elements_and_svMOD(6,year,month,day, hour, min, sec);

% differenza di secondi dalla data di riferimento
j_oggi = J0(year, month, day);
secondi_oggi = hour*60*60 + min*60 + sec;
j_ref = J0(2028,1,1);
secondi_ref = 9*60*60 + 59*60;
delta_secondi = secondi_ref - secondi_oggi;
delta_giorni = j_oggi - j_ref;
delta_secondi = delta_giorni*24*60*60 + delta_secondi;


mu = 3.793074669799999e7; % saturn

% coe enceladus
e = 0;
RA = saturn_coe(3);
incl = saturn_coe(4) + 0.009;
w = 0;
orbital_period = 7 + 53*60 + 32*60*60;  %orbital_period = 32 h 53 min 07 s
a = 238000; %Km
velocita = 2*pi*a/orbital_period;
h = a*velocita;

omega = 2*pi/orbital_period;

TA = 96.4*pi/180 - omega*delta_secondi;

coe_enc = [h e RA incl w TA];

[r, v] = sv_from_coe(coe_enc,mu);
r = r + saturn_r;
v = v + saturn_v;

%test per valutare semiasse maggiore effettivo
% TA = 0;
% coe_enc = [h e RA incl w TA];
% [r1, v1] = sv_from_coe(coe_enc,mu);
% TA = 180*pi/180;
% coe_enc = [h e RA incl w TA];
% [r2, v2] = sv_from_coe(coe_enc,mu);
% r5 = norm(r2-r1)/2;





% %% trovare TA del 13.5.2028 alle 7:50:00    ----> TA =   214.62*pi/180              
% figure(12)
% hold on
% hold on
% plot(0,0,'xr')
% plot(saturn_r(1),saturn_r(2),'ob')
% A = [0 saturn_r(1)]; 
% B = [0 saturn_r(2)]; 
% line(A,B)
% 
% TA = 214.62*pi/180;
% 
% coe_enc = [h e RA incl w TA];
% 
% [r1, v1] = sv_from_coe(coe_enc,mu);
% 
% 
% r1 =r1 + saturn_r;
% 
% plot(r1(1),r1(2),'or')
% 
% xlim([ 1.1368*10^9, 1.1374*10^9])
% ylim([7.848*10^8, 7.851*10^8])
% 
% hold off



 