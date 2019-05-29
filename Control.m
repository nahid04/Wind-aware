clear all;
clc;
syms v_com r_com x_com y_com z_com r x y z yaw_com yaw r_dot_com

r_com = transpose([x_com y_com z_com]);

r = transpose ([x y z]);
c1= 0.5;


if c1 > 0
v_com = c1*[r_com - r] + r_dot_com

else v_com = r_dot_com
    
end



