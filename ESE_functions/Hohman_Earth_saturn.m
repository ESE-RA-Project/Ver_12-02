%% Hohmann transfer from Earth to Saturn.
%% Initialization.

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
                 
colors = ["g"          %green
          "m"          %magenta
          "b"          %blue
          "r"          %red
          "#A2142F"    %darker red
          "#7E2F8E"    %purple
          "#4DBEEE"    %darker cyan
          "c"          %(bright) cyan
         "#D95319"    %orange
          "#77AC30"    %darker green
          "#D95319"];  %orange, not visible due to Sun orbit dimensions
                 

%% Interplanetary data.
%Distances from the Sun.
Earth_to_Sun = distances(3);%0.14960e9 [km]
Saturn_to_Sun = distances(6);%1.42667e9 [km]

%Radii
Saturn_radius = radii(6); % 2.499e+02 [km]
Earth_radius = radii(3); % 6371 [km]

%Masses
Earth_mass = masses(3); %5.97219e24; %[kg]
Saturn_mass= masses(6); %5.68319e26
Sun_mass = masses(11);%1.989e30; %[kg]

%SOI: (m_planet/m_Sun)^(2/5) * distance_from_Sun
Earth_SOI = (Earth_mass/Sun_mass)^(2/5)*Earth_to_Sun;%9.24522e5; %[km]
Saturn_SOI = (Saturn_mass/Sun_mass)^(2/5)*Saturn_to_Sun; %5.4538e7 %[km]

[xx,yy,zz] = sphere(10);

mu_sun = G*masses(11); %1.327565122000000e+11 %[km^3/s^2]

%% Orbit parameters
Eorb_park_dep=200; %[km]
Sorb_park_arr=400; %[km]

%Inner orbit radius.
r1E=Eorb_park_dep + Earth_radius + Earth_to_Sun; %  149604833 [km]
%Outer orbit radius.
r2S=Sorb_park_arr + Saturn_radius + Saturn_to_Sun ; % 1.426725054000000e+09 [km]

%Velocities close to Earth.
Ve_init = sqrt(mu_sun/r1E); %  29.788943805379283 [km/s]
Ve_trans_p = sqrt((2*mu_sun*r2S)/(r1E*(r1E+r2S))); % 40.078982921885647 [km/s]

%First thrust.
Delta_v1E = Ve_trans_p - Ve_init; % 10.290039116506364 [km/s]


%Velocities close to Saturn.
Vs_trans_a = sqrt((2*mu_sun*r1E)/(r2S*(r1E+r2S))); % 4.202638434103334 [km/s]
Vs_fin = sqrt(mu_sun/r2S); %   9.646233561216789 [km/s]

%Second Thrust.
Delta_v2S = Vs_fin - Vs_trans_a; %5.443595127113455 [km/s]

%Transfer period.
Delta_tE2S = pi*sqrt(((r1E+r2S)^3)/(8*mu_sun)); %1.907864938411930e+08 [s]
Delta_tE2Syear = Delta_tE2S /(3600*24*365); %6.049800033016012 [years]
fprintf('La durata del trasferimento orbitale vale 6 anni 18 giorni e 6 ore. \n')

%Transfer orbit (ellipse) parameters.
a2 = (r1E + r2S)/2; % 7.881649435000000e+08 [km]
e2 = (r2S - r1E)/(r1E+r2S); % 0.810185882747270

%Total thrust.
Delta_vES = Delta_v1E + Delta_v2S; %15.733634243619818 [km/s]
fprintf('Il delta(v) totale secondo Hohmann (trascurando le SOI dei pianeti) vale 15.734 km/s.\n ')

%% Transfer orbit plot.
c2=e2*a2;
b2=sqrt(a2*a2-c2*c2);

if exist('figure2') == 0
    figure()
else
    figure2()
end

%Heliocentric frame.
body_sphere_MY(11,[0,0,0]);
hold on
plot3(0,0,0,'ro')

xlabel('x')
ylabel('y')
zlabel('z')

%Transfer orbit.
h=ellipse(a2,b2,0,c2,0); 
% Transfer orbit coe.   
origin_coe=coe_from_sv([-r1E, 0 ,0], [0 , -Ve_trans_p, 0], mu_sun);

%Earth.
body_sphere_MY(3,[-Earth_to_Sun,0,0]);
plot3(-Earth_to_Sun,0,0,'bo');
%Earth orbit.
orb_earth = ellipse(r1E,r1E,0,0,0,'k');
%Initial parking orbit.
Eparkorb_init = ellipse(Eorb_park_dep+Earth_radius,Eorb_park_dep+Earth_radius,0,-Earth_to_Sun, 0);
%Initial spacecraft position.
plot3(-r1E,0,0,'k^');



%Saturn.
body_sphere_MY(6,[Saturn_to_Sun,0,0]);
plot3(Saturn_to_Sun,0,0,'yo');
%Saturn orbit.
orb_saturn = ellipse(r2S,r2S,0,0,0,'k');
%Final parking orbit.
Sparkorb_final = ellipse(Sorb_park_arr,Sorb_park_arr,0,0,Saturn_to_Sun + Saturn_radius);
%Final spacecraft position.
plot3(r2S,0,0,'^');
%Saturn SOI.
surface(Saturn_to_Sun+Saturn_SOI*xx, Saturn_SOI*yy,...
        Saturn_SOI*zz,'FaceColor','none','EdgeColor',colors(6))
%Earth SOI.    
surface(-Earth_to_Sun+Earth_SOI*xx, Earth_SOI*yy,...
        Earth_SOI*zz,'FaceColor','none','EdgeColor',colors(3))



%% Zoom on Earth.
if exist('figure2') == 0
    figure()
else
    figure2()
end

hold on

xlim([-149611204, -149592091])
ylim([-9.5565e+03, 9.5565e+03])
zlim([-9.5565e+03, 9.5565e+03])

xlabel('x')
ylabel('y')
zlabel('z')
view(-10,45)

%Earth.
body_sphere_MY(3,[-Earth_to_Sun,0,0]);
%Earth orbit.
orb_earth = ellipse(r1E,r1E,0,0,0,'k');
% % %Initial parking orbit.
% % Eparkorb_init = ellipse(Eorb_park_dep,Eorb_park_dep,0,0,-Earth_to_Sun - Earth_radius);
%Initial parking orbit.
Eparkorb_init = ellipse(Eorb_park_dep+Earth_radius,Eorb_park_dep+Earth_radius,0,-Earth_to_Sun, 0);

%Initial spacecraft position.
plot3(-r1E,0,0,'k^');
%Transfer orbit.
h=ellipse(a2,b2,0,c2,0); 
%Earth SOI.
surface(-Earth_to_Sun+Earth_SOI*xx, Earth_SOI*yy,...
        Earth_SOI*zz,'FaceColor','none','EdgeColor',colors(3))
                
%Escape from Earth.
orbit_esc=[h.XData(1),h.YData(1),0;h.XData(2),h.YData(2),0];
[esc_hyp, delta_vesc] = Hohmann_esc_hyp(orbit_esc, Eorb_park_dep, origin_coe, [0 , -Ve_trans_p, 0]);
%% Exiting Earth SOI.
if exist('figure2') == 0
    figure()
else
    figure2()
end

hold on
xlim([-1.509850453107824e+08, -1.482114786892176e+08])
ylim([-1.386783310782394e+06, 1.386783310782394e+06])
zlim([-1.386783310782394e+06, 1.386783310782394e+06])

grid
xlabel('x')
ylabel('y')
zlabel('z')
view(-10,45)

body_sphere_MY(3,[-Earth_to_Sun,0,0]);
%Initial spacecraft position.
plot3(-r1E,0,0,'k^');
%Transfer orbit.
h=ellipse(a2,b2,0,c2,0); 
%Earth orbit.
orb_earth = ellipse(r1E,r1E,0,0,0,'k');
% %Initial parking orbit.
% Eparkorb_init = ellipse(Eorb_park_dep+Earth_radius,Eorb_park_dep+Earth_radius,0,0,-Earth_to_Sun );

%Earth SOI.
surface(-Earth_to_Sun+Earth_SOI*xx, Earth_SOI*yy,...
        Earth_SOI*zz,'FaceColor','none','EdgeColor',colors(3))
    
%Escape from Earth.
orbit_esc=[h.XData(1),h.YData(1),0;h.XData(2),h.YData(2),0];
[esc_hyp, delta_vesc] = Hohmann_esc_hyp(orbit_esc, Eorb_park_dep, origin_coe, [0 , -Ve_trans_p, 0]);



%% Zoom on Saturn.

if exist('figure2') == 0
    figure()
else
    figure2()
end

hold on
view(10,10)
grid

xlim([1.426491726000000e+09, 1.426841118000000e+09])
ylim([-174696, 174696])
zlim([-174696, 174696])


%Saturn.
body_sphere_MY(6,[Saturn_to_Sun,0,0]);
%Saturn orbit.
orb_saturn = ellipse(r2S,r2S,0,0,0,'k');
%Final parking orbit.
Sparkorb_final = ellipse(Sorb_park_arr + Saturn_radius,Sorb_park_arr + Saturn_radius,0,Saturn_to_Sun ,0);
%Final spacecraft position.
plot3(r2S,0,0,'^');
%Saturn SOI.
surface(Saturn_to_Sun+Saturn_SOI*xx, Saturn_SOI*yy,...
        Saturn_SOI*zz,'FaceColor','none','EdgeColor',colors(6))

%Transfer orbit.
h=ellipse(a2,b2,0,c2,0); 

%Capture on Saturn.
orbit_cap=[h.XData(end-1),h.YData(end-1),0;h.XData(end),h.YData(end),0];
[cap_hyp, delta_vcap] = Hohmann_cap_hyp(orbit_cap, Sorb_park_arr, origin_coe, [0, Vs_trans_a,0]);


%% Exiting Saturn SOI.
if exist('figure2') == 0
    figure()
else
    figure2()
end

hold on
view(-10,45)

xlim([1.344859901299280e+09, 1.508472942700720e+09])
ylim([-8.180652070072033e+07, 8.180652070072033e+07])
zlim([-8.180652070072033e+07, 8.180652070072033e+07])


%Saturn.
body_sphere_MY(6,[Saturn_to_Sun,0,0]);
%Saturn orbit.
orb_saturn = ellipse(r2S,r2S,0,0,0,'k');
%Final parking orbit.
Sparkorb_final = ellipse(Sorb_park_arr + Saturn_radius,Sorb_park_arr + Saturn_radius,0,Saturn_to_Sun ,0);
%Final spacecraft position.
plot3(r2S,0,0,'^');
%Saturn SOI.
surface(Saturn_to_Sun+Saturn_SOI*xx, Saturn_SOI*yy,...
        Saturn_SOI*zz,'FaceColor','none','EdgeColor',colors(6))
%Transfer orbit.
h=ellipse(a2,b2,0,c2,0); 

%Capture on Saturn.
orbit_cap=[h.XData(end-1),h.YData(end-1),0;h.XData(end),h.YData(end),0];
[cap_hyp, delta_vcap] = Hohmann_cap_hyp(orbit_cap, Sorb_park_arr, origin_coe, [0, Vs_trans_a,0]);

%Effective Delta_v
Eff_delta_v= delta_vesc + delta_vcap;
fprintf('Il delta(v) effettivo (considerando le fasi iperboliche di fuga e di cattura) vale 18.229 km/s.\n ')

