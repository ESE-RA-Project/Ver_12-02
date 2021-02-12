function [traj, delta_v] = ...
  Hohmann_cap_hyp(orbit, park_r, origin_coe, v_in)
% CAPTURE_HYP(goal_id, orbit, arr_time, park_r, park_i, origin_coe, v_in)
%   computes the trajectory the spacecraft will follow
%   to enter the sphere of influence of object GOAL_ID at time
%   ARR_TIME and to position itself into a circular parking orbit of 
%   radius PARK_R.
%
%   It uses the last two points of the arrival interplanetary orbit 
%   (computed via the patched conics method) stored in ORBIT to generate
%   a trajectory with orbital elements ORIGIN_COE that will allow it
%   to reach velocity V_IN at the entrance of the body's SOI.
%
%   [traj, delta_v] = CAPTURE_HYP(...) returns the computed capture
%   trajectory TRAJ and the needed change of velocity DELTA_V.
%
%   Options for the Sun are not contemplated since that would be
%   the general case of an interplanetary trajectory.
%
%   goal_id  - identifier of the destination planet:
%                    1 = Mercury
%                    2 = Venus
%                    3 = Earth
%                    4 = Mars
%                    5 = Jupiter
%                    6 = Saturn
%                    7 = Uranus
%                    8 = Neptune
%                    9 = Pluto
%                   10 = Enceladus
%                   11=Sun
%
%   orbit    - last two points of the interplanetary trajectory
%              computed via the patched conics method
%
%   arr_time - array specifying time of arrival with elements 
%                    (in this order):
%                     year         - range: 1901 - 2099
%                     month        - range: 1 - 12
%                     day          - range: 1 - 31
%                     hour         - range: 0 - 23
%                     minute       - range: 0 - 60
%                     second       - range: 0 - 60
%
%   park_r   - radius of the circular parking orbit around the body
%
%   origin_coe - classical orbital elements of the origin interplanetary
%                 orbit:
%                h    = angular momentum (km^2/s)
%                e    = eccentricity
%                RA   = right ascension of the ascending
%                       node (rad)
%                incl = inclination of the orbit (rad)
%                w    = argument of perigee (rad)
%                TA   = true anomaly (rad)
%                a    = semimajor axis (km)
%
%   v_in - velocity of the spacecraft at the entry into the SOI
%          of the planet


    %% Data
    mu = 1.327565122000000e+11;

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

G = 6.6742e-20; %[km^3/kg/s^2]
    
    %% Input data
    %SOI: (m_planet/m_Sun)^(2/5) * distance_from_Sun
    pl_SOI = (masses(6)/masses(11))^(2/5)...
        * distances(6); %[km]
    
    pl_mu = G * masses(6); %[km^3/s^2]
    pl_radius = radii(6); %[km]
    
    %% Needed variables

    v_arr = [0 ,sqrt(mu /distances(6)), 0];
    pl_r0= [distances(6), 0 , 0];
    vinf = norm(v_arr-v_in); % [km/s^2] , v-infinity

    rp = pl_radius + park_r; % [km] perigee
    e = 1 + rp*vinf^2/pl_mu; % eccentricity
    a = rp/(e-1); % [km] semi-major axis
    
    % target radius for the right hyperbola (Fig. 8.13-8.14 Curtis)
    Delta = rp*sqrt(1+2*pl_mu/(rp*vinf^2)); % [km] aiming radius
    
    %Angle between arrival and departure branch of the hyperbola
    half_delta = asin(1/e);
    
    v_hyp = sqrt(vinf^2+2*pl_mu/rp);%sqrt(-2*pl_mu/rp + pl_mu/a);
    %vc = sqrt(2*pl_mu/rp);? Eq 8.59 with circular parking orbit
    vc = sqrt(pl_mu/rp); % velocity of parking orbit

    h = Delta*vinf; % [km^2/s] specific angular momentum
    RA = origin_coe(3); %[rad] right ascension of ascending node
    incl = origin_coe(4)+pi; %[rad] inclination
    w = 0;%origin_coe(5);
    TA = 0;%origin_coe(6);

    n = sqrt(pl_mu/a^3); % mean motion
    
    %% Trajectory computation
    hyp = zeros(100*24*3600/60,3);
    
    %Angle and direction of the orbit at the entering point
    in_dir = (orbit(1,1:3)-orbit(2,1:3))'; %exit vector: (1,1:3)<-(2,1:3)
    in_angle = deg2rad(atan2d_0_360(in_dir(2),in_dir(1))) + 6.9 *pi/180;

    %Initial point
    coe = [h, e, RA, incl, w, TA];
    [rprova,~] = sv_from_coe(coe, pl_mu);
    alpha = deg2rad(atan2d_0_360(rprova(2),rprova(1)));
    
    %Desired characteristics to align with the interplanetary orbit
    xi_des = in_angle - pi - half_delta;
    alpha_des = pi/2 + xi_des;
    w_des = alpha_des - alpha;
    
    counter = 1;
    for t=0:60:100*24*3600%+ 60*100
        %Using Kepler's method to compute the point
        M = n*t;
        F = kepler_H(e,M);
        cosf = (cosh(F)-e)/(1-e*cosh(F)); %(Eq. 3.41b) Curtis
        f = acos(cosf);
        
       
        coe = [h, e, RA, incl, w_des+pi,f];
        [r,~] = sv_from_coe(coe, pl_mu);
        
        %Sometimes the above algorithm produces NaN elements because of the
        %kepler_H function (it seems not to be able to deal with high
        %numbers)
        if(any(isnan(r)))
            if(hyp(2,:) ~= [0 0 0]) %from the third point onward
                diff = hyp(counter-1,:)-hyp(counter-2,:);
                diff = 60*norm(v_in)*diff/norm(diff);
                point = hyp(counter-1,:)' + diff';
            else %first two points
                
                    coe = [h, e, RA, incl, w_des+pi,-t/6]; %w aveva +pi
               
                [peri,~] = sv_from_coe(coe, pl_mu);
                point = pl_r0' + peri';
            end
        else %if kepler_H returned a valid result
             point = pl_r0' + r';
        end
        
        %Adding the point to the list
        hyp(counter,:) = point';
        counter = counter+1;
        
        %To stop computing when the SOI is exited
        if (all(hyp(1,:) ~= [0 0 0]) && norm(point'-hyp(1,:))>=  pl_SOI) %5.453768046714688e+07
            break;
        end
    end
    %Getting rid of unused elements in the array
    hyp = hyp(1:counter-2, 1:3);
    
    %% Hyperbola plot
    plot3(hyp(:,1),hyp(:,2),hyp(:,3),'m-')
    % view Saturn
%     hold on
%     xlim([1.36e9 1.5e9])
%     ylim([-6e7 6e7])
%     hold off

%view earth
%  hold on
%     xlim([-1.4962e8 -1.4958e8])
%     ylim([-2e4 2e4])
%     hold off
    
    %% Output arguments
    traj = flip(hyp);
    delta_v = v_hyp - vc;
end    