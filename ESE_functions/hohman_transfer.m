mu = 1.327565122000000e+11; %[km^3/s^2]

% nasa's eye
dist_t_s = 92966632.4; % mile
dist_s_sat = 861215270.7; % mile

dist_t_s = dist_t_s*1.60934; %Km
dist_s_sat = dist_s_sat*1.60934; %Km

% wiki
% dist_t_s = 152000000; %Km
% dist_s_sat = 1429000000; % Km

orb_t = 200; % Km
orb_sat = 400; % Km

R_t = 3671; %Km
R_sat = 58232; %Km

r1 = dist_t_s + orb_t + R_t;
r2 = dist_s_sat + orb_sat + R_sat;

v_init = sqrt(mu/r1);
v_trans_p = sqrt((2*mu*r2)/(r1*(r1+r2)));

Delta_v1 = v_trans_p - v_init;

v_trans_a = sqrt((2*mu*r1)/(r2*(r1+r2)));
v_fin = sqrt(mu/r2);

Delta_v2 = v_fin - v_trans_a;

Delta_t = pi*sqrt(((r1+r2)^3)/(8*mu));

a = (r1 + r2)/2;
e = (r2 - r1)/(r1+r2);

Delta_v = Delta_v1 + Delta_v2;

%% plot 
c=e*a;
b=sqrt(a*a-c*c);
h=ellipse(a,b,0,c,0)
hold on
orb_init = ellipse(r1,r1,0,0,0,'k')
orb_fin = ellipse(r2,r2,0,0,0,'k')
plot(0,0,'r+') %sole

