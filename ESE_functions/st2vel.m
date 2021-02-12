%%%%Velocità da r e t%%%%
function orb_vel=st2vel(orb_r, t, v_init)
vel=zeros(length(t), 3);
vel(1, :)=v_init;

%integration
for i=1:length(t)-1
    vel(i+1, :)=(orb_r(i+1, :)-orb_r(i, :))/(t(i+1)-t(i));
end

orb_vel= vel;
end