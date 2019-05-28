syms r v x y z m Rb dcd dct e3 fd dw Xw g esv esw Re wx wy wz ct vx vy vz ft w p q r ew j j1 j2 j3 i wm om n jm tau tau_d fi

r = transpose([x y z])
diff(r) = v
v = transpose([vx vy vz])
e3 = transpose([0 0 1])
g = transpose([0 0 -g])
 dcd = -0.5 * ([vx.*vx;vy.*vy;vz.*vz]) .* Re .* esv 
 dct = 0.5 * ([wx.*wx;wy.*wy;wz.*wz]) .* Re .* esw
 cd = int (dcd)
 ct = int (dct)
fd = - Rb * cd .* ([vx.*vx;vy.*vy;vz.*vz] )
tau_d = - ct .* ([wx.*wx;wy.*wy;wz.*wz] )
diff(v) = (Rb .* e3 .* ft)/m + fd ./ m + g

j = diag([j1 j2 j3])
w = transpose([p q r]) 
ew = transpose([-q p 0]) 
Xw = -transpose (w)  
wm = symsum(om * ((-1)^ i) , i, 1, n) % wrong, right it after finishing
ft = symsum(fi , i, 1, n)
jm = 1
dw = ((cross ( -w ,  j .* w ) ) ./ j )  + ((jm * (wm .* ew))./ j ) + ( tau ./ j) + (tau_d ./ j)
