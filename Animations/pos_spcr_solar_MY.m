% This script calculates a vector containing each row position of
% Dawn_spacecraft with respect to solar system

% init
pos_spcr = zeros(n_days, 3);


%% Earth to Venus
% ev_days = datenum([2022 12 2]) - datenum([2022 8 2]);
ev_days = J0(2022, 12, 2) - J0(2022, 8 ,2);

% [body_pos1, sp_v1, body_posf1, sp_vf1,tof1, orb_elem1] = ...
%                 gen_orbitMOD2(3,2,[data0.year data0.month data0.day 0 0 0],[data1.year data1.month data1.day 0 0 0],0);
%             
% Ev_orbit = intpl_orbit(tof1,Earth_r0,sp_v1);

ev_pos = Ev_orbit(:,1:3);

q = floor(size(ev_pos,1)/ev_days); % rate of samples for each day
counter=0;

for i = 1:ev_days
	pos_spcr(i,:) = ev_pos( i*q ,:);
%     counter=counter+1;
end

%% Venus to Earth
ve_days = J0(2023, 5, 2)- J0(2022 ,12 ,2);

%[body_pos2, sp_v2, body_posf2, sp_vf2,tof2, orb_elem2] = ...
%                gen_orbitMOD2(2,3,[data1.year data1.month data1.day,0,0,0],[data2.year data2.month data2.day,0,0,0],0);
%Ve_orbit = intpl_orbit(tof2,Venus_r1,sp_v2);

ve_pos = Ve_orbit(:,1:3);

q = size(ve_pos,1)/ve_days; % rate of samples for each day


for i = 1:ve_days
	pos_spcr(i+ev_days,:) = ve_pos(floor(q * i),:);
%     counter=counter+1;
end


%% Earth to Saturn
es_days = J0(2028 ,2 ,2)- J0(2023, 5, 2);
es_pos = Es_orbit(:,1:3);

q = size(es_pos,1)/es_days; % rate of samples for each day

for i = 1:es_days
	pos_spcr(i+ev_days+ve_days,:) = es_pos( floor(q*i ),:) ;
%     counter=counter+1;
end

%% Saturn Park orbit (2028/5/28 - 02/02/28)
s_days = J0(2028, 5, 28)- J0(2028, 2, 2);

for i = (ev_days + ve_days +es_days+ 1):(ev_days + ve_days +es_days+ s_days)
	
	year_now	= time_vector(i,1);
	month_now	= time_vector(i,2);
	day_now		= time_vector(i,3);
        counter=counter+1;
	
	[~, tmp_Saturn_now, ~, ~] = planet_elements_and_svMOD(6, ...
								year_now, month_now, day_now, 0, 0, 0);
	
	pos_spcr(i,:) = tmp_Saturn_now;
    
end

% Saturn to Enceladus 

se_pos = orb_SE(:,1:3)
se_days=1;


for i=0:1
pos_spcr(i+1+ev_days+ve_days+es_days+s_days,:) = se_pos( end-1+i,:) ;
end


%tmp
pos_spcr(end,:) = pos_spcr(end-1,:);

%% tests
% %sanity check
% (ev_days + ve_days + es_days + s_days +se_days) - (n_days-1);


%% days

day_left_earth	= 1;
day_venus       = 1 + ev_days;
day_earth       = 1 + ev_days + ve_days;
day_saturn      = 1 + ev_days + ve_days + es_days;
day_left_saturn = 1 + ev_days + ve_days + es_days + s_days;
day_enceladus   = 1 + ev_days + ve_days + es_days + s_days + se_days;
