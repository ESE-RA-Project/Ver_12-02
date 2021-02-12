function [traj, delta_v] = escape_hyp_MY(orbit, park_r, goal_coe, v_out)
% ESCAPE_HYP(planet_id,goal_id,orbit,dep_time,park_r,park_i,goal_coe,v_out)
%   computes the trajectory the spacecraft will follow
%   to escape the sphere of influece of the object OBJ_ID AT
%   DEP_TIME.
%
%   It uses the first points of the departing interplanetary orbit 
%   (computed via the patched conics method) stored in ORBIT and the
%   parking orbit radius PARK_R and inclination PARK_I to generate
%   a trajectory with orbital elements GOAL_COE that will allow it
%   to reach velocity V_OUT at the end of the travel.
%
%   [traj, delta_v] = ESCAPE_HYP(...) returns the computed escape
%   trajectory TRAJ and the needed change of velocity DELTA_V.
%
%   Options for the Sun are not contemplated since that would be
%   the general case of an interplanetary trajectory.
%
%   obj_id  - identifier of the origin planet:
%                            1 = Mercury
%                            2 = Venus
%                            3 = Earth
%                            4 = Mars
%                            5 = Jupiter
%                            6 = Saturn
%                            7 = Uranus
%                            8 = Neptune
%                            9 = Pluto
%                           10 = Enceladus
%                           11 = Sun
%
%   orbit    - first two points of the interplanetary trajectory
%              computed via the patched conics method
%
%   dep_time - array specifying time of departure with elements 
%                    (in this order):
%                     year         - range: 1901 - 2099
%                     month        - range: 1 - 12
%                     day          - range: 1 - 31
%                     hour         - range: 0 - 23
%                     minute       - range: 0 - 60
%                     second       - range: 0 - 60
%
%   park_r   - radius of the circular parking orbit around origin planet
%
%   goal_coe - classical orbital elements of the target interplanetary
%              orbit:
%                h    = angular momentum (km^2/s)
%                e    = eccentricity
%                RA   = right ascension of the ascending
%                       node (rad)
%                incl = inclination of the orbit (rad)
%                w    = argument of perigee (rad)
%                TA   = true anomaly (rad)
%                a    = semimajor axis (km)
%
%   v_out - escape velocity from the object SOI

    %% Argument validation
%     validateattributes(goal_coe,{'double'},{'size',[1 7]})

    %% Constants
    mu= 1.327565122000000e+11;

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
                 413690250];%[km]

    G = 6.6742e-20; %[km^3/kg/s^2]
    
    %% Input data
    %SOI: (m_planet/m_Sun)^(2/5) * distance_from_Sun
    pl_SOI = (masses(3)/masses(11))^(2/5)...
        * distances(3); %[km]
    
    pl_mu = G * masses(3); %[km^3/s^2]
    pl_radius = radii(3); %[km]
    
    park_i = goal_coe(4);

    %% Needed variables
    pl_r0= [-distances(3), 0 , 0];
    v_dep = [0,-sqrt(mu /distances(3)),0];

    %v-infinity of the departure hyperbola
    vinf = norm(v_out - v_dep); %[km/s]

    %Hyperbola characteristics
    rp = pl_radius + park_r; %[km], periapsis
    e = 1 + rp*vinf^2/pl_mu; %eccentricity
    a = rp/(e-1); %[km], semi-major axis
%     b = a*sqrt(e^2-1); %[km], semi-minor axis

    %Velocity at hyperbola periapsis
    vp = sqrt(vinf^2 + 2*pl_mu/rp);

    %Angle between hyperbola center and exiting-branch
%     beta = acos(1/e);
    half_delta = asin(1/e);

    %Hyperbola orbital elements
    h = rp*vp;
    RA = goal_coe(3);
    incl = goal_coe(4)+2*pi;
    w = goal_coe(5);
    TA = goal_coe(6);

    %Mean motion
    n = sqrt(pl_mu/a^3);
    
    %Gravitational parameter of the origin planet
    mu_dep = masses(3) * G; %[km^3/s^2]
    %Velocity of the spacecraft after burn
    v_b = sqrt(vinf^2 + 2*mu_dep/rp); %[km/s]
    %Velocity of circular parking orbit around origin planet
    v_park = sqrt(mu_dep/rp);
    
    %% Trajectory computation
    rr = zeros(100*24*3600/60,3);
    out_dir = (orbit(2,1:3)-orbit(1,1:3))'; %exit vector: (2,1:3)<-(1,1:3)
    out_angle = deg2rad(atan2d_0_360(out_dir(2),out_dir(1)))+0.185*pi/180;
%     up_angle = deg2rad(-atan2d_0_360(out_dir(3),-out_dir(1)));

%     p = -a*(1-e^2); %semilatum [km]
%     tran = acos(1/e-rp/p); %true anomaly at perigee
    
%     coe = [h, e, RA, incl, 0, 0];
    coe = [h, e, RA, incl, w, TA];
    [rprova,~] = sv_from_coe(coe, pl_mu);
    alpha = deg2rad(atan2d_0_360(rprova(2),rprova(1)));

    xi_des = out_angle - pi - half_delta;
    alpha_des = pi/2 + xi_des;
    w_des = alpha_des - alpha;        
    counter = 1;
    for t=0:60:100*24*3600%ceil(pl_SOI/norm(v_out))
        M = n*t; %Hyperbolic mean anomaly
        F = kepler_H(e,M); %Hyperbolic eccentric anomaly
        cosf = (e-cosh(F))/(e*cosh(F)-1);
        f = acos(cosf); %True anomaly
        coe = [h, e, RA, incl, w_des, f]; %TA = f
        [r,~] = sv_from_coe(coe, pl_mu); %spacecraft position
        if(any(isnan(r)))
            if(rr(2,:) ~= [0 0 0])
                diff = rr(counter-1,:)-rr(counter-2,:);
                diff = 60*norm(v_out)*diff/norm(diff);
                point = rr(counter-1,:)' + diff';
            else
                coe = [h, e, RA, incl, w_des, t/6];
                [peri,~] = sv_from_coe(coe, pl_mu);
                point = pl_r0' + peri';
            end
        else
             point = pl_r0' + r';
        end
        rr(counter, :) = point';
        counter = counter + 1;
        %rr = cat(1,rr,point');
        if all(rr(1,:) ~= [0 0 0]) && norm(point'-rr(1,:))>= pl_SOI
            break;
        end
    end
    rr = rr(1:counter-2, 1:3);
    %% Deleted because useless, I think
    %Angle of orientation of escape trajectory, to be aligned with
    %the escape velocity vector
%     out_dir = Rotz(goal_coe(3))'*Rotx(park_i)'*...
%         (orbit(2,1:3)-orbit(1,1:3))'; %exit vector: (2,1:3)<-(1,1:3)
%     out_angle = deg2rad(atan2d_0_360(out_dir(2),out_dir(1)));
% 
%     t = 0:0.1:5;
% 
% %     Parametric hyperbola equations
%     xh_l = -a*cosh(t);
%     xh_r = a*cosh(t);
%     yh = b*sinh(t);
% 
%     hyp = [];
%     for i = 1:length(t)
%         point = pl_r0' + Rotz(goal_coe(3))*Rotx(park_i)*...
%                     Rotz(out_angle)*Rotz(beta)*...
%                     ([xh_l(i); -yh(i);0] + [-(rp-a);0;0]);
%         hyp = cat(1,hyp,point');
%         if norm(hyp(size(hyp,1),:)-hyp(1,:))>= pl_SOI
%             break;
%         end
%     end
    
    %% Hyperbola plot
    plot3(rr(:,1),rr(:,2),rr(:,3),'m-')%plot3(hyp(:,1),hyp(:,2),hyp(:,3),'m-')
    
    %% Output arguments
    traj = rr;%hyp;
    delta_v = v_b - v_park;
end