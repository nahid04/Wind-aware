clear all;
clc;
syms ac pwbcx psbwx dpwcb f tf u ux0 uy0 uz0  uxf uyf uzf dpwb Re tw jerk psvb df  ctx dx pwbcy dotv T dotw dpsv dpsw psbwy cty dy pwbcz psbwz ctz dz vabx vaby vabz bx by bz pv dpv dpsv r ga gaw psw pwbx pwby pwbz psvby psvbx psvbz cdy cdz cdx ax ay az hw ew dpsw esw psv v h ev yaw pitch roll dv x y z m tau_x tau_y tau_z ctx cty ctz wb Rb wbx wby wbz vbx vby vbz cdx cdy cdz vb n dcd dct e3 fd dw sv Xw g esv esw Re wx wy wz ct vx vy vz ft w p q r ew j j1 j2 j3 i wm sw om n jm tau tau_d fi
syms performance_index gx gy gz dx dy dz aox vox pox alphax betax gammax pfx vfx afx dpx dax dvx copx cosx  aoy voy poy alphay betay gammay pfy vfy afy dpy day dvy copy cosy aoz voz poz alphaz betaz gammaz pfz vfz afz dpz daz dvz copz cosz
syms gset next_r p0 p1 p2 p3 p4 inverse_diagonal Vmax maximum_value
% In general -------------------------------------------------------

Rb =  [cos(yaw)*cos(pitch) cos(yaw)*sin(roll)*sin(pitch)-cos(roll)*sin(yaw) sin(roll)*sin(yaw)+cos(roll)*cos(yaw)*sin(pitch); ...
     cos(pitch)*sin(yaw) cos(roll)*cos(yaw)+sin(roll)*sin(yaw)*sin(pitch) cos(roll)*sin(yaw)*sin(pitch)-cos(yaw)*sin(roll); ...
     -sin(pitch) cos(pitch)*sin(roll) cos(roll)*cos(pitch)];
 
 
Re  = [ cos(yaw)*cos(pitch) cos(pitch)*sin(yaw) -sin(pitch); cos(yaw)*sin(roll)*sin(pitch)-cos(roll)*sin(yaw) cos(roll)*cos(yaw)+sin(roll)*sin(yaw)*sin(pitch) cos(pitch)*sin(roll); ...
    sin(roll)*sin(yaw)+cos(roll)*cos(yaw)*sin(pitch) cos(roll)*sin(yaw)*sin(pitch)-cos(yaw)*sin(roll) cos(roll)*cos(pitch)]

ft = symsum(fi , i, 1, n);
cd = [cdx 0 0 ; 0  cdy 0 ; 0 0 cdz];
r = transpose([x y z]);
vb = transpose ([vbx vby vbz]);
n = norm (vb);
 fd = - Rb * ([vabx * bx * cdx; vaby * by * cdy ; vabz * bz * cdz] )
p = transpose([vx vy vz]);
e3 = transpose([0 0 1]);
g = transpose([0 0 -g]) ;
% diff(r) = v
% predicted sv ----------------------------------------------------------------------------
psv = (-1/m) * ( Rb )  * (cd ) * (vb ) * (n)


% general equation---------------------------------

psvb = Re * psv;
psvbx                                                                  = psvb (1);
psvby = psvb (2);
psvbz = psvb (3);

j = diag([j1 j2 j3]);
w = transpose([p q r]); 
bw = transpose([-q p 0]); 
Xw = [0 -w(3) w(2) ; w(3) 0 -w(1) ; -w(2) w(1) 0 ] ;
wm = symsum(om * ((-1)^ i) , i, 1, n); % wrong, right it after finishing;
bx= norm (vabx)
by= norm (vaby)
bz= norm (vabz)
jm = 1;
tau = transpose ([tau_x tau_y tau_z]);
%dw = ((cross ( -w ,  j .* w ) ) ./ j )  + ((jm * (wm .* ew))./ j ) + ( tau ./ j) + (tau_d ./ j)
k = norm (wb);
ct = [ctx 0 0 ; 0  cty 0 ; 0 0 ctz];
wb = transpose ([wbx wby wbz]);
ax = norm (psvbx);
ay = norm (psvby);
az = norm (psvbz);
% predicted rotational velocity-------------------------------
psw = (-inv(j) * ct * wb * k)

% predicted cd---------
cdx   = (m * ax) / (bx)^2 ; 
cdy  = (m * ay) / (by)^2; 
cdz   = (m * az) / (bz)^2; 

% Predicted Wb-------------------------

pwbx = vbx - ( sign (psvbx)) * sqrt ((m/cdx) * ax)
pwby = vby - ( sign (psvby)) * sqrt ((m/cdy) * ay)
pwbz = vbz - ( sign (psvbz)) * sqrt ((m/cdz) * az)

% Predicted Wc---------------------
dx = norm (psbwx);
dy = norm (psbwy);
dz = norm (psbwz);

pwbcx = p - sign(psbwx) * sqrt ((j1/ctx) * dx)
pwbcy = q - sign(psbwy) * sqrt ((j2/cty) * dy)
pwbcz = r - sign(psbwz) * sqrt ((j3/ctz) * dz)

% dot of psv-----------------------------------------

dpsv = ((-1/m) * ( Rb ) * Xw  * (cd ) * (vb ) * (n)) - (2 * 1/m * Rb * cd * (diff (vb)) * (n) )
dpsw = ( - (2 * inv (j)) * ct * (diff (wb)) * (k) )

%dot of pwb -----------------------------------------

dpwb = ((dpsv) +  ((1/m) * ( Rb ) * Xw  * (cd ) * (vb ) * (n)) + (2 * 1/m * Rb * cd * (diff (vb)) * (n) ) ./  (2 * 1/m * Rb * cd * (n) ))

%dot of pwcb -----------------------------------------

dpwcb = ((psw) + (2 * inv (j)) * ct * (diff (wb)) * (k) ) ./  ( (2 * inv (j)) * ct * (k) )

%----------------------------------------------------------------------------------------------------------
%Acceleration-----------------------------------------

dotv = ( (1/m ) * ft * Rb * e3 )+ g + psv

%Angular rate change--------------------------------

dotw = (cross (( -(inv(j)) * w ) , (j * w) ))  + jm * wm * (inv(j)) * bw + (inv(j)) * tau + psw 


%thrust------------------------------------------------
ac  = diff(diff(r))

T = ac - g - psv

ft = norm (T) * m


% jerk + differential of mass normalised
% thrust---------------------------------------------

jerk =  diff (diff(diff(r)));
df = diff (ft/m); 
jerk = df * Rb * e3 + (ft/m) * Rb * Xw * e3 + psv

df = transpose (e3) * Re * (jerk - psv) 

% Angular rate -------------------------
tw = [wx -wy 0];
tw = (ft/m) * [1 0 0 ; 0 1 0 ; 0 0 0] * Re * (jerk-psv)

% axis motion primitives in terms of jerk --------------------

u = transpose([ux0 uy0 uz0 uxf uyf uzf]);
f = u.^2;
performance_index =  int (f , 0, tf);


%{
optimal cost function gx gy gz and optimal path dx dy dz
tf =3;
pox= 0;
vox =0;
aox =0;
pfx = 2;
vfx = 1;
afx = 0;
poy= 0;
voy =0;
aoy =0;
pfy = 2;
vfy = 1;
afy = 0;
poz= 0;
voz =0;
aoz =0;
pfz = 2;
vfz = 1;
afz = 0;

%-- x coordinate
dpx = pfx-pox-(vox*tf)-(1/2*aox*tf^2);
dvx = vfx-vox-aox*tf;
dax = afx-aox;
copx = [dpx; dvx; dax];

alphax = ((1/tf^5) * (720*dpx-360*tf*dvx+60*tf^2*dax)) ;
betax = ((1/tf^5) * (-360*tf*dpx+168*tf^2*dvx-24*tf^3*dax));
gammax = ((1/tf^5) * (60*tf^2*dpx-24*tf^3*dvx+3*tf^4*dax));

cosx= [alphax; betax; gammax] ;

dx = [((alphax/120) * tf^5) + ((betax/24) * tf^4) + ((gammax/6) * tf^3) + ((aox/2)*tf^2) + (vox*tf) + pox;
    ((alphax/24) * tf^4) + ((betax/6) * tf^3) + ((gammax/2) * tf^2) + (aox*tf) + vox;
    ((alphax/6) * tf^3) + ((betax/2) * tf^2) + ((gammax) * tf) + aox];


gx = gammax^2 + betax*gammax*tf + (1/3)*betax^2*tf^2+(1/3)*alphax*gammax*tf^2+(1/4)*alphax*betax*tf^3+(1/20)*alphax^2*tf^4;

%--- y coordinate

dpy = pfy-poy-(voy*tf)-(1/2*aoy*tf^2);
dvy = vfy-voy-aoy*tf;
day = afy-aoy;
copy = [dpy; dvy; day];

alphay = ((1/tf^5) * (720*dpy-360*tf*dvy+60*tf^2*day)) ;
betay = ((1/tf^5) * (-360*tf*dpy+168*tf^2*dvy-24*tf^3*day));
gammay = ((1/tf^5) * (60*tf^2*dpy-24*tf^3*dvy+3*tf^4*day));

cosy= [alphay; betay; gammay] ;

dy = [((alphay/120) * tf^5) + ((betay/24) * tf^4) + ((gammay/6) * tf^3) + ((aoy/2)*tf^2) + (voy*tf) + poy;
    ((alphay/24) * tf^4) + ((betay/6) * tf^3) + ((gammay/2) * tf^2) + (aoy*tf) + voy;
    ((alphay/6) * tf^3) + ((betay/2) * tf^2) + ((gammay) * tf) + aoy];


gy = gammay^2 + betay*gammay*tf + (1/3)*betay^2*tf^2+(1/3)*alphay*gammay*tf^2+(1/4)*alphay*betay*tf^3+(1/20)*alphay^2*tf^4;

%--- z coordinate

dpz = pfz-poz-(voz*tf)-(1/2*aoz*tf^2);
dvz = vfz-voz-aoz*tf;
daz = afz-aoz;
copz = [dpz; dvz; daz];

alphaz = ((1/tf^5) * (720*dpz-360*tf*dvz+60*tf^2*daz)) ;
betaz = ((1/tf^5) * (-360*tf*dpz+168*tf^2*dvz-24*tf^3*daz));
gammaz = ((1/tf^5) * (60*tf^2*dpz-24*tf^3*dvz+3*tf^4*daz));

cosz= [alphaz; betaz; gammaz]  ;

dz = [((alphaz/120) * tf^5) + ((betaz/24) * tf^4) + ((gammaz/6) * tf^3) + ((aoz/2)*tf^2) + (voz*tf) + poz;
    ((alphaz/24) * tf^4) + ((betaz/6) * tf^3) + ((gammaz/2) * tf^2) + (aoz*tf) + voz;
    ((alphaz/6) * tf^3) + ((betaz/2) * tf^2) + ((gammaz) * tf) + aoz];


gz = gammaz^2 + betaz*gammaz*tf + (1/3)*betaz^2*tf^2+(1/3)*alphaz*gammaz*tf^2+(1/4)*alphaz*betaz*tf^3+(1/20)*alphaz^2*tf^4;

next_r = [dx(1) dy(1) dz(1)]';
gset= [gx; gy; gz];

%}

%---Random waypoints

p =  floor((rand (8,3))*100);


    dx = sym('dx',[3 1])

%-------------------- Vmax

%inverse_diagonal = [1/(cd*(ft/m)*(e3) + (Re*g)) 0 0; 0 (1/(cd*(ft/m)*(e3) + (Re*g))) 0; 0 0 (1/(cd*(ft/m)*(e3) + (Re*g)))]
%maximum_value = ([[1 1 1]' * (abs(inverse_diagonal))])
%(Vmax)^2 = max ([[1 1 1]' * (abs(inverse_diagonal))])
Vmax = 10;
syms tgo 
% distant measurement---------------
i = 0;
for n = 2:8
    i = i+1;
    d(i) = norm (p(n,:) - p(n-1,:));
    tgo = d / Vmax
    tf = tgo(i);
pox(i)= p(n-1,1);
vox(i) =1;
aox(i) =1;
pfx(i) = p(n,1);
vfx(i) = 1;
afx(i)= 0;
poy(i)= p(n-1,2);
voy(i) =1;
aoy(i) =1;
pfy(i) = p(n,2);
vfy(i) = 1;
afy(i) = 0;
poz(i)= p(n-1,3);
voz(i) =1;
aoz(i) =1;
pfz(i) = p(n,3);
vfz(i) = 1;
afz(i) = 0;

%-- x coordinate
dpx(i) = pfx(i)-pox(i)-(vox(i)*tf)-(1/2*aox(i)*tf^2);
dvx(i) = vfx(i)-vox(i)-aox(i)*tf;
dax(i) = afx(i)-aox(i);
%copx(i) = [dpx(i) dvx(i) dax(i)];

alphax(i) = ((1/tf^5) * (720*dpx(i)-360*tf*dvx(i)+60*tf^2*dax(i))); 
betax(i) = ((1/tf^5) * (-360*tf*dpx(i)+168*tf^2*dvx(i)-24*tf^3*dax(i)));
gammax(i) = ((1/tf^5) * (60*tf^2*dpx(i)-24*tf^3*dvx(i)+3*tf^4*dax(i)));

%cosx(i)= [alphax(i); betax(i); gammax(i)];
dxx(i) = [((alphax(i)/120) * tf^5) + ((betax(i)/24) * tf^4) + ((gammax(i)/6) * tf^3) + ((aox(i)/2)*tf^2) + (vox(i)*tf) + pox(i)];
dxv(i) = [   ((alphax(i)/24) * tf^4) + ((betax(i)/6) * tf^3) + ((gammax(i)/2) * tf^2) + (aox(i)*tf) + vox(i)];
 dxa(i) = [  ((alphax(i)/6) * tf^3) + ((betax(i)/2) * tf^2) + ((gammax(i)) * tf) + aox(i)];


%gxx(i) = gammax(i)^2 + betax(i)*gammax(i)*tf + (1/3)*betax(i)^2*tf^2+(1/3)*alphax(i)*gammax(i)*tf^2+(1/4)*alphax(i)*betax(i)*tf^3+(1/20)*alphax(i)^2*tf^4;


%--- y coordinate

dpy(i) = pfy(i)-poy(i)-(voy(i)*tf)-(1/2*aoy(i)*tf^2);
dvy(i) = vfy(i)-voy(i)-aoy(i)*tf;
day(i) = afy(i)-aoy(i);
%copy(i) = [dpy(i); dvy(i); day(i)];

alphay(i) = ((1/tf^5) * (720*dpy(i)-360*tf*dvy(i)+60*tf^2*day(i))); 
betay(i) = ((1/tf^5) * (-360*tf*dpy(i)+168*tf^2*dvy(i)-24*tf^3*day(i)));
gammay(i) = ((1/tf^5) * (60*tf^2*dpy(i)-24*tf^3*dvy(i)+3*tf^4*day(i)));

%cosy(i)= [alphay(i); betay(i); gammay(i)]; 

dyx(i) = [((alphay(i)/120) * tf^5) + ((betay(i)/24) * tf^4) + ((gammay(i)/6) * tf^3) + ((aoy(i)/2)*tf^2) + (voy(i)*tf) + poy(i)];
 dyv(i) = [((alphay(i)/24) * tf^4) + ((betay(i)/6) * tf^3) + ((gammay(i)/2) * tf^2) + (aoy(i)*tf) + voy(i)];
 dya(i) =  [((alphay(i)/6) * tf^3) + ((betay(i)/2) * tf^2) + ((gammay(i)) * tf) + aoy(i)];


%gy(i) = gammay(i)^2 + betay(i)*gammay(i)*tf + (1/3)*betay(i)^2*tf^2+(1/3)*alphay(i)*gammay(i)*tf^2+(1/4)*alphay(i)*betay(i)*tf^3+(1/20)*alphay(i)^2*tf^4

%--- z coordinate

dpz(i) = pfz(i)-poz(i)-(voz(i)*tf)-(1/2*aoz(i)*tf^2);
dvz(i) = vfz(i)-voz(i)-aoz(i)*tf;
daz(i) = afz(i)-aoz(i);
%copz(i) = [dpz(i); dvz(i); daz(i)]

alphaz(i) = ((1/tf^5) * (720*dpz(i)-360*tf*dvz(i)+60*tf^2*daz(i))) 
betaz(i) = ((1/tf^5) * (-360*tf*dpz(i)+168*tf^2*dvz(i)-24*tf^3*daz(i)))
gammaz(i) = ((1/tf^5) * (60*tf^2*dpz(i)-24*tf^3*dvz(i)+3*tf^4*daz(i)))

%cosz(i)= [alphaz(i); betaz(i); gammaz(i)]  

dzx(i) = [((alphaz(i)/120) * tf^5) + ((betaz(i)/24) * tf^4) + ((gammaz(i)/6) * tf^3) + ((aoz(i)/2)*tf^2) + (voz(i)*tf) + poz(i)];
dzv(i) =    [((alphaz(i)/24) * tf^4) + ((betaz(i)/6) * tf^3) + ((gammaz(i)/2) * tf^2) + (aoz(i)*tf) + voz(i)];
 dza(i) =  [ ((alphaz(i)/6) * tf^3) + ((betaz(i)/2) * tf^2) + ((gammaz(i)) * tf) + aoz(i)];


%gz(i) = gammaz(i)^2 + betaz(i)*gammaz(i)*tf + (1/3)*betaz(i)^2*tf^2+(1/3)*alphaz(i)*gammaz(i)*tf^2+(1/4)*alphaz(i)*betaz(i)*tf^3+(1/20)*alphaz(i)^2*tf^4

next_r(i,1:3) = [dxx(i) dyx(i) dzx(i)]
next_v(i,1:3) = [dxv(i) dyv(i) dzv(i)];
next_a(i,1:3) = [dxa(i) dya(i) dza(i)];
%plot3(next_r,'o')
figure(1)
plot3(p(:,1),p(:,2),p(:,3),'-')
figure(2)
 plot3(next_r(:,1),next_r(:,2),next_r(:,3),'-')
 
%gset= [gx; gy; gz];
   end
%p =  floor((rand (15,3))*100);

