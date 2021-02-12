function [r,v,Delta_t] = orbit_hohmann(Saturn_pos,Saturn_coe, Enceladus_pos,Enceladus_coe)

dist = Enceladus_pos - Saturn_pos;

G = 6.6742e-20; %[km^3/kg/s^2]

masses = 10^24 * [0.330104
                  4.86732
                  5.97219
                  0.641693
                  1898.13
                  568.319
                  86.8103
                  102.410
                  0.01309
                  8.6*1e-5
                  1989100]; %[kg]

radii = [2439.7
         6051.8 
         6371
         3389.5
         69911
         58232
         25362
         24622
         1151
         249.9
         695508]; %[km] 
     
distances = [57909227
             108209475
             149598262
             227943824
             778340821
             1426666422
             2870658186
             4498396441
             1426904442
             413690250
                     0];%[km]
                 
dist = Enceladus_pos - Saturn_pos;
% Sat_Enc=2.38*10^5; %Km
Sat_Enc=norm(dist); %Km
mu_sat = G*masses(6); %Saturn
orb_park_dep=400;
orb_park_arr=100;

Enceladus_SOI_Sat = (masses(10)/masses(6))^(2/5)*Sat_Enc; %[km]


r1=orb_park_dep + radii(6);
r2=orb_park_arr + radii(10) + Sat_Enc ;

v_init = sqrt(mu_sat/r1);
v_trans_p = sqrt((2*mu_sat*r2)/(r1*(r1+r2)));

Delta_v1 = v_trans_p - v_init;

v_trans_a = sqrt((2*mu_sat*r1)/(r2*(r1+r2)));
v_fin = sqrt(mu_sat/r2);

Delta_v2 = v_fin - v_trans_a;

Delta_t = pi*sqrt(((r1+r2)^3)/(8*mu_sat));

a = (r1 + r2)/2;
e = (r2 - r1)/(r1+r2);

Delta_v = Delta_v1 + Delta_v2;


%%

vers_dist = dist/norm(dist);
ref = [1 0 0];
theta = acos(vers_dist*ref')*180/pi;
theta = theta + pi;


% coe hohmann

h = r1*v_trans_p;
ecc = e;
RA = Enceladus_coe(3) + 206.87*pi/180;
RA = 320.45*pi/180;
incl = Saturn_coe(4) + 0.009;
% RA = theta;
% incl = 0;
w = 0;
TA = 0;

coe_h = [h ecc RA incl w TA];

[r0, v0] = sv_from_coe(coe_h,mu_sat);
i=0;
r = zeros(300,3);
v= zeros(300,3);
for TA = 0:0.6:180
    
    i=i+1;
    TA = TA*pi/180;
    coe_h = [h ecc RA incl w TA];
    [r(i,:), v(i,:)] = sv_from_coe(coe_h,mu_sat);
    
end
r = r + Saturn_pos;

hold on
plot3(r(:,1),r(:,2),r(:,3),'r')