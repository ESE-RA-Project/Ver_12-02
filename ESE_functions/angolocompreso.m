function [teta, alpha_x,beta_y,gamma_z, R]=angolocompreso(v1,v2)
%v1 e v2 espressi come vettori riga
radtodeg=180/pi;
vers_v1=v1/norm(v1);
% i1=vers_v1(1)*[1 0 0]';
% j1=vers_v1(2)*[0 1 0]';
% k1=vers_v1(3)*[0 0 1]';

vers_v2=v2/norm(v2);
% i2=vers_v2(1)*[1 0 0]';
% j2=vers_v2(2)*[0 1 0]';
% k2=vers_v2(3)*[0 0 1]';

cos_teta=vers_v1*vers_v2';
teta=radtodeg*acos(cos_teta);

R=vers_v2'*vers_v1;

%inversione yaw-pitch-roll
beta_y=radtodeg*atan2(-R(3,1),sqrt(R(3,2)^2 + R(3,3)^2));
alpha_x=radtodeg*atan2(R(3,2),R(3,3));
gamma_z=radtodeg*atan2(R(2,1),R(1,1));

end

% alpha_x,beta_y,gamma_z