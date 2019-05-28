syms pv dpv dpsv r ga gaw psw hw ew dpsw esw psv v h ev yaw pitch roll dv x y z m tau_x tau_y tau_z ctx cty ctz wb Rb wbx wby wbz vbx vby vbz cdx cdy cdz vb n dcd dct e3 fd dw sv Xw g esv esw Re wx wy wz ct vx vy vz ft w p q r ew j j1 j2 j3 i wm sw om n jm tau tau_d fi
% In general -------------------------------------------------------

Rb =  [cos(yaw)*cos(pitch) cos(yaw)*sin(roll)*sin(pitch)-cos(roll)*sin(yaw) sin(roll)*sin(yaw)+cos(roll)*cos(yaw)*sin(pitch); ...
     cos(pitch)*sin(yaw) cos(roll)*cos(yaw)+sin(roll)*sin(yaw)*sin(pitch) cos(roll)*sin(yaw)*sin(pitch)-cos(yaw)*sin(roll); ...
     -sin(pitch) cos(pitch)*sin(roll) cos(roll)*cos(pitch)];
ft = symsum(fi , i, 1, n);
cd = [cdx 0 0 ; 0  cdy 0 ; 0 0 cdz];
%r = transpose([x y z]);
vb = transpose ([vbx vby vbz]);
n = norm (vb);
 
v = transpose([vx vy vz]);
e3 = transpose([0 0 1]);
g = transpose([0 0 -g]) ;
% diff(r) = v
% velocity----------------------------------------------------------------------------
sv = (-1/m) * ( Rb )  * (cd ) * (vb ) * (n);
dv = ( (1/m ) * ft * Rb * e3 )+ g + sv 
v = int(dv)

% estimated velocity----------------------------------------------------------------------------
ga = 0.5;
h = 0.5;
dpsv = ga * ev;
psv  = int (dpsv);
dpv = ( (1/m ) * ft * Rb * e3 )+ g + psv + h * ev;
pv = int(dpv)
% general equation---------------------------------

j = diag([j1 j2 j3]);
w = transpose([p q r]); 
bw = transpose([-q p 0]); 
Xw = [0 -w(3) w(2) ; w(3) 0 -w(1) ; -w(2) w(1) 0 ] ;
wm = symsum(om * ((-1)^ i) , i, 1, n); % wrong, right it after finishing;

jm = 1;
tau = transpose ([tau_x tau_y tau_z]);
%dw = ((cross ( -w ,  j .* w ) ) ./ j )  + ((jm * (wm .* ew))./ j ) + ( tau ./ j) + (tau_d ./ j)
k = norm (wb);
ct = [ctx 0 0 ; 0  cty 0 ; 0 0 ctz];
wb = transpose ([wbx wby wbz]);
% rotational velocity-------------------------------
sw = (-inv(j) * ct * wb * k);
dw = (cross (( -(inv(j)) * w ) , (j * w) ))  + jm * wm * (inv(j)) * bw + (inv(j)) * tau + sw ;
w = int(dw);

% estimated rotational velocity ----------------------------------------
hw = 0.5;
ew = 0.5;
dpsw = gaw * ew;
psw = int (dpsw)
dpw = (cross (( -(inv(j)) * w ) , (j * w) ))  + jm * wm * (inv(j)) * bw + (inv(j)) * tau + psw + hw*ew ;
pw = int(dpw)
