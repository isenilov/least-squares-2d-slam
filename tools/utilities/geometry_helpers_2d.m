% computes the pose 2d pose vector v from an homogeneous transform A
% A:[ R t ] 3x3 homogeneous transformation matrix, r translation vector
% v: [x,y,theta]  2D pose vector
function v=t2v(A)
	v(1:2, 1)=A(1:2,3);
	v(3,1)=atan2(A(2,1),A(1,1));
end

% computes the homogeneous transform matrix A of the pose vector v
% A:[ R t ] 3x3 homogeneous transformation matrix, r translation vector
% v: [x,y,theta]  2D pose vector
function A=v2t(v)
  	c=cos(v(3));
  	s=sin(v(3));
	A=[c, -s, v(1) ;
	s,  c, v(2) ;
	0   0  1  ];
end


% normalizes and angle between -pi and pi
% th: input angle
% o: output angle
function o = normalizeAngle(th)
	o = atan2(sin(th),cos(th));
end

################### ADDED FOR LEAST SQUARES 2D (work in progress) ###############
#rotation matrix for 2d
function R=R2d(rot)
 c=cos(rot);
 s=sin(rot);
 R= [c, -s;
     s, c];
endfunction

#derivative of rotation matrix for 2d
function R=R_prime(rot)
 dc=-sin(rot); #derivative of cos(rot(x)
 ds=cos(rot);  #derivative of sin(rot(x)
 R= [dc, -ds;
     ds, dc];
endfunction


function v=flattenIsometryByColumns(T)
v=zeros(6,1);
v(1:4)=reshape(T(1:2,1:2),4,1);
v(5:6)=T(1:2,3);
endfunction

function T=unflattenIsometryByColumns(v)
  T=eye(3);
  T(1:2,1:2)=reshape(v(1:4),2,2);
  T(1:2,3)=v(5:6);
endfunction


#derivative of rotation matrix in 0
global  R0=[0, -1;
	           1, 0];

############ Taken from 3d, not yet rewrited #################

#{
function S=skew(v)
  S=[0,    -v(3), v(2);
     v(3),  0,    -v(1);
     -v(2), v(1), 0];
endfunction

function v=flattenIsometry(T)
v=zeros(12,1);
v(1:9)=reshape(T(1:3,1:3)',9,1);
v(10:12)=T(1:3,4);
endfunction

function T=unflattenIsometry(v)
  T=eye(4);
  T(1:3,1:3)=reshape(v(1:9),3,3)';
  T(1:3,4)=v(10:12);
endfunction

function M=flatTransformationMatrix(v)
  T=unflattenIsometry(v);
  R=T(1:3,1:3);
  t=T(1:3,4);
  M=eye(12);
  M(1:3,1:3)=R';
  M(4:6,4:6)=R';
  M(7:9,7:9)=R';
  M(10,1:3)=t';
  M(11,4:6)=t';
  M(12,7:9)=t';
endfunction;
#}
