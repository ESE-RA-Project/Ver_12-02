% trovare TA di riferimento

% % %% trovare TA del 13.5.2028 alle 7:50:00    ----> TA = 101*pi/180;   
% (old 214.62*pi/180 , non so perchè sia cambiato)
% 
% % pos saturno nella data stabilita
% [saturn_coe, saturn_r, saturn_v, ~] =...
%                         planet_elements_and_svMOD(6,2028,5,13, 7, 50, 0);
% 
% mu = 3.793074669799999e7; % saturn
% 
% % coe enceladus
% e = 0;
% RA = saturn_coe(3);
% incl = saturn_coe(4) + 0.009;
% 
% 
% 
% w = 0;
% orbital_period = 7 + 53*60 + 32*60*60;  %orbital_period = 32 h 53 min 07 s
% a = 238000; %Km
% velocita = 2*pi*a/orbital_period;
% 
% h = a*velocita;
% 
% TA = 101*pi/180;
% 
% coe_enc = [h e RA incl w TA];
% 
% [r, v] = sv_from_coe(coe_enc,mu);
% r = r + saturn_r;
% v = v + saturn_v;
%            
% figure(12)
% hold on
% hold on
% plot(0,0,'xr')
% plot(saturn_r(1),saturn_r(2),'ob')
% A = [0 saturn_r(1)]; 
% B = [0 saturn_r(2)]; 
% line(A,B)
% 
% plot(r(1),r(2),'or')
% 
% xlim([ 1.1368*10^9, 1.1374*10^9])
% ylim([7.848*10^8, 7.851*10^8])
% 
% hold off


% %% trovare TA del 1.1.2028 alle 9:59:00    ----> TA = 96.4*pi/180;

% pos saturno nella data stabilita
[saturn_coe, saturn_r, saturn_v, ~] =...
                        planet_elements_and_svMOD(6,2028,1,1, 9, 59, 0);

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

TA = 96.4*pi/180;

coe_enc = [h e RA incl w TA];

[r, v] = sv_from_coe(coe_enc,mu);
r = r + saturn_r;
v = v + saturn_v;
           
figure(12)
hold on
hold on
plot(0,0,'xr')
plot(saturn_r(1),saturn_r(2),'ob')
A = [0 saturn_r(1)]; 
B = [0 saturn_r(2)]; 
line(A,B)

plot(r(1),r(2),'or')

xlim([ 1.201*10^9, 1.203*10^9])
ylim([6.909*10^8, 6.919*10^8])

hold off