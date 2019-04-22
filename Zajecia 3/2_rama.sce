clc;clear;cls;
function [Ke,fe]=beam2e(ex,ey,ep,eq);
[nargout,nargin] = argn(0);
ieee(1)
// 
//    INPUT:  ex = [x1 x2]
//
//            eq = [qx qy]       distributed loads, local directions
// 
//    OUTPUT: Ke : element stiffness matrix (6 x 6)
//
//            fe : element load vector (6 x 1)
  b=[ ex(2)-ex(1); ey(2)-ey(1) ];
  L=sqrt(b'*b);  n=b/L;

  E=ep(1);  A=ep(2);  I=ep(3);
 
  qx=0; qy=0;  if nargin>3; qx=eq(1); qy=eq(2); end

  Kle=[E*A/L   0            0      -E*A/L      0          0 ;
         0   12*E*I/L^3   6*E*I/L^2  0   -12*E*I/L^3  6*E*I/L^2;
         0   6*E*I/L^2    4*E*I/L    0   -6*E*I/L^2   2*E*I/L;
       -E*A/L  0            0       E*A/L      0          0 ;
         0   -12*E*I/L^3 -6*E*I/L^2  0   12*E*I/L^3  -6*E*I/L^2;
         0   6*E*I/L^2    2*E*I/L    0   -6*E*I/L^2   4*E*I/L];
   
  fle=L*[qx/2 qy/2 qy*L/12 qx/2 qy/2 -qy*L/12]';

  G=[n(1) n(2)  0    0    0   0;
    -n(2) n(1)  0    0    0   0;
      0    0    1    0    0   0;
      0    0    0   n(1) n(2) 0;
      0    0    0  -n(2) n(1) 0;
      0    0    0    0    0   1];

  Ke=G'*Kle*G;   fe=G'*fle; 
endfunction


Edof=[1 1 2 3 4 5 6;  2 10 11 12 7 8 9; 3 4 5 6 7 8 9];
K=zeros(12,12); f=zeros(12,1); //f(4)=2000;

A1=45.3e-4; A2=142.8e-4;

I1=2510e-8; I2=33090e-8; E=2.1e11;
ep1=[E A1 I1]; ep3=[E A2 I2];
ex1=[0 0]; ex2=[6 6]; ex3=[0 6]; ey1=[0 4]; ey2=[0 4]; ey3=[4 4];
eq1=[0 0]; eq2=[0 0]; eq3=[0 -10000];

Ke1=beam2e(ex1,ey1,ep1); Ke2=beam2e(ex2,ey2,ep1); 
[Ke3,fe3]=beam2e(ex3,ey3,ep3,eq3);

function [K,f] = assem(edof,K,Ke,f,fe) 

// Number of arguments in function call
[nargout,nargin] = argn(0)

// Display warning for floating point exception
ieee(1)

// K=assem(edof,K,Ke)
// [K,f]=assem(edof,K,Ke,f,fe)

// PURPOSE
//  Assemble element matrices Ke ( and fe ) into the global
//  stiffness matrix K ( and the global force vector f )
//  according to the topology matrix edof.
// 
// INPUT: edof:       dof topology matrix
//        K :         the global stiffness matrix
//        Ke:         element stiffness matrix
//        f :         the global force vector
//        fe:         element force vector
// 
// OUTPUT: K :        the new global stiffness matrix
//         f :        the new global force vector

[nie,n] = size(edof);
t = edof(:,2:n);
for i = 1:nie;
   K(t(i,:),t(i,:)) = K(t(i,:),t(i,:))+Ke;
  if nargin==5 then
    f(t(i,:))=f(t(i,:))+fe;
  end;
end;
//--------------------------end--------------------------------
endfunction

K=assem(Edof(1,:),K,Ke1); K=assem(Edof(2,:),K,Ke2);
[K,f]=assem(Edof(3,:),K,Ke3,f,fe3);

function [d,Q]=solveq(K,f,bc)
  // Number of arguments in function call
[nargout,nargin] = argn(0);

// Display warning for floating point exception
ieee(1)
// a=solveq(K,f)
// [a,Q]=solveq(K,f,bc)
//-------------------------------------------------------------
// PURPOSE
//  Solve static FE-equations considering boundary conditions.
//
// INPUT: K : global stiffness matrix, dim(K)= nd x nd
//        f : global load vector, dim(f)= nd x 1
//
//        bc : boundary condition matrix
//            dim(bc)= nbc x 2, nbc : number of b.c.'s
//
// OUTPUT:  a : solution including boundary values
//          Q : reaction force vector
//              dim(a)=dim(Q)= nd x 1, nd : number of dof's
//-------------------------------------------------------------

// LAST MODIFIED: M Ristinmaa   1993-10-06
// Copyright (c)  Division of Structural Mechanics and
//                Department of Solid Mechanics.
//                Lund Institute of Technology
//-------------------------------------------------------------
  if nargin==2 ; 
     d=K\f ; 
  elseif nargin==3;
     [nd,nd]=size(K);
     fdof=[1:nd]';
//
     dd=size(fdof);
     d=zeros(dd(1),dd(2))
     Q=zeros(dd(1),dd(2));
//
     pdof=bc(:,1);
     dp=bc(:,2);
     fdof(pdof)=[];
//
     s=K(fdof,fdof)\(f(fdof)-K(fdof,pdof)*dp);
//
     d(pdof)=dp;
     d(fdof)=s;
  end 
     Q=K*d-f;
   
endfunction

function [es,edi,eci] = beam2s(ex,ey,ep,ed,eq,n)

// Ouput variables initialisation (not found in input variables)
es=[];
edi=[];
eci=[];
// Number of arguments in function call
[nargout,nargin] = argn(0)
// Display warning for floating point exception
ieee(1)
//            ed = [u1 ... u6] element displacements
//            eq = [qx qy]     distributed loads, local directions 
//            n : number of evaluation points ( default=2 )
//          
//    OUTPUT: es = [ N1 V1 M1 ;  section forces, local directions, in 
//                   N2 V2 M2 ;  n points along the beam, dim(es)= n x 3
//                   .........]  
//           
//            edi = [ u1 v1 ;    element displacements, local directions,
//                    u2 v2 ;    in n points along the beam, dim(es)= n x 2
//                   .......]    
// 
//            eci = [ x1  ;      local x-coordinates of the evaluation 
//                    x2 ;       points, (x1=0 and xn=L)
//                    ...]
EA = ep(1)*ep(2);
EI = ep(1)*ep(3);
b = [ex(2)-ex(1);ey(2)-ey(1)]
L = sqrt(b'*b);

if max(size(ed(:,1)))>1 then
  disp("Only one row is allowed in the ed matrix !!!");
  return;
end;

qx = 0;
qy = 0;
if nargin>4 then
  qx = eq(1);
  qy = eq(2);
end;

ne = 2;
if nargin>5 then
  ne = n;
end;

C = [0,0,0,1,0,0;
     0,0,0,0,0,1;
     0,0,0,0,1,0;
     L,0,0,1,0,0;
     0,L^3,L^2,0,L,1;
     0,3*(L^2),2*L,0,1,0];

n = b/L;

G = [n(1),n(2),0,0,0,0;
     -n(2),n(1),0,0,0,0;
     0,0,1,0,0,0;
     0,0,0,n(1),n(2),0;
     0,0,0,-n(2),n(1),0;
     0,0,0,0,0,1];

M = inv(C)*(G*ed'-[0,0,0,-(qx*(L^2))/(2*EA),(qy*(L^4))/(24*EI),(qy*(L^3))/(6*EI)]');

A = [M(1),M(4)]';
B = [M(2),M(3),M(5),M(6)]';

x = (0:L/(ne-1):L)';
%v0 = size(x);
zero = zeros(%v0(1),%v0(2));
%v0 = size(x);
one = ones(%v0(1),%v0(2));

u = [x,one]*A-((x .^2)*qx)/(2*EA);
du = [one,zero]*A-(x*qx)/EA;
v = ([x .^3,x .^2,x,one]*B+((x .^4)*qy)/(24*EI));
// dv=[3*x.^2 2*x one zero]*B+(x.^3)*qy/(6*EI);
d2v = ([6*x,2*one,zero,zero]*B+((x .^2)*qy)/(2*EI));
d3v = ([6*one,zero,zero,zero]*B+(x*qy)/EI);

N = EA*du;
M = EI*d2v;
V = -EI*d3v;
es = [N,V,M];
edi = [u,v];
eci = x;
endfunction

  function [ed]=extract(edof,a)
// ed=extract(edof,a)
//-------------------------------------------------------------
// PURPOSE
//  Extract element displacements from the global displacement
//  vector according to the topology matrix edof.
//
// INPUT:   a:  the global displacement vector
//
//         edof:  topology matrix
//
// OUTPUT: ed:  element displacement matrix
//-------------------------------------------------------------

// LAST MODIFIED: M Ristinmaa 1993-08-24
// Copyright (c) 1991-94 by Division of Structural Mechanics and
//                          Department of Solid Mechanics.
//                          Lund Institute of Technology
//-------------------------------------------------------------
    [nie,n]=size(edof);
//
    t=edof(:,2:n);
//
    for i = 1:nie
        ed(i,1:(n-1))=a(t(i,:))';
    end
//
endfunction

function [magnfac]=eldia2(ex,ey,es,eci,magnfac,magnitude)
[nargout,nargin] = argn(0);
ieee(1)

// INPUT:
//	ex = [ x1 x2 ]
//	ey = [ y1 y2 ]	element node coordinates.
//
//	es = [  X1;
//		X2;
//		... ] 	vector containing the section force
//			in Nbr evaluation points along the beam.
//	
//	eci = [	x1;
//		x2;
//		... ]	local x-coordinates of the evaluation points (Nbr).
//
//	magnfac = [scalar]	magnification factor.

if ~((nargin==4)|(nargin==5)|(nargin==6))
    disp('??? Wrong number of input arguments!')
    return
 end
a=size(ex); b=size(ey);
 
 if ~((a-b)==[0 0])
    disp('??? Check size of coordinate input arguments!')
   return
 end
c=size(es);
a=size(eci); 
 if ~(c(1)==a(1)|(c(2)==1))
    disp('??? Check size of (es) or (eci) input arguments!')
   return
 end 
b=[ex(2)-ex(1);ey(2)-ey(1)];
Length=sqrt(b'*b);
n=b/Length;
if nargin==4
	// Determining magnfac.
	magnfac=max(abs(es))/(0.2*Length);
else
   disp("wykres")
end
es=es/magnfac;
if nargin>5
	// Creating a scaled bar.
	N=magnitude(1);
	L=N/magnfac;
	mag=size(magnitude);
	if mag(2)==1
		x=0;
		y=-0.5;
	else
		x=magnitude(2);
		y=magnitude(3);
	end;
	plot2d([x (x+L)],[y y],5,rect=[-2 -0.5 8 6.0]);
	plot2d([x x],[(y-L/20) (y+L/20)],5,rect=[-2 -0.5 8 6.0]);
	plot2d([(x+L/2) (x+L/2)],[(y-L/20) (y+L/20)],5,rect=[-2 -0.5 8 6.0]);
	plot2d([(x+L) (x+L)],[(y-L/20) (y+L/20)],5,rect=[-2 -0.5 8 6.0]);
	xstring((N*1.1)/magnfac+x,y,sprintf('%4.0f ', N));
end
// From local x-coordinates to global coordinates of the beam element.
An=size(eci);
Nbr=An(1);
A=zeros(Nbr,2)
A(1,1)=ex(1);
A(1,2)=ey(1);
for i=2:Nbr
	A(i,1)=A(1,1)+eci(i)*n(1);
	A(i,2)=A(1,2)+eci(i)*n(2);
end

B=A;
for i=1:Nbr
	A(i,1)=A(i,1)+es(i)*n(2);
	A(i,2)=A(i,2)-es(i)*n(1);
end
// Type is set to 0 by default.
pl=0;
if nargin==4
	// Do not plot if input arguments are missing.
	pl=3;
end
if pl==0
	// Plots striped diagrams.
	for i=1:Nbr
		plot2d([B(i,1) A(i,1)],[B(i,2) A(i,2)],5,rect=[-2 -0.5 8 6.0]);
	end
end
if pl<3
	plot2d(A(:,1),A(:,2),5,rect=[-2 -0.5 8 6.0]);
	plot2d([ex(1) A(1,1)],[ey(1) A(1,2)],5,rect=[-2 -0.5 8 6.0]);
	plot2d([ex(2) A(Nbr,1)],[ey(2) A(Nbr,2)],5,rect=[-2 -0.5 8 6.0]);
	
	plot2d(ex',ey',rect=[-2 -0.5 8 6.0]);
	
end
endfunction

bc=[1 0;2 0;3 0;10 0;11 0;12 0];
a=solveq(K,f,bc)

Ed=extract(Edof,a);
[es1,edi1,eci1]=beam2s(ex1,ey1,ep1,Ed(1,:),eq1,10);
[es2,edi2,eci2]=beam2s(ex2,ey2,ep1,Ed(2,:),eq2,10);
[es3,edi3,eci3]=beam2s(ex3,ey3,ep3,Ed(3,:),eq3,10);

subplot(2,2,1)
magnfac=eldia2(ex1,ey1,es1(:,1),eci1);
magnitude=[3e5 0.5 0];
eldia2(ex1,ey1,es1(:,1),eci1,magnfac);
eldia2(ex2,ey2,es2(:,1),eci2,magnfac);
eldia2(ex3,ey3,es3(:,1),eci3,magnfac,magnitude);

subplot(2,2,2)
magnfac=eldia2(ex3,ey3,es3(:,2),eci3);
magnitude=[3e5 0.5 0];
eldia2(ex1,ey1,es1(:,2),eci1,magnfac);
eldia2(ex2,ey2,es2(:,2),eci2,magnfac);
eldia2(ex3,ey3,es3(:,2),eci3,magnfac,magnitude);

subplot(2,2,3)
magnfac=eldia2(ex3,ey3,es3(:,3),eci3);
magnitude=[3e5 0.5 0];
eldia2(ex1,ey1,es1(:,3),eci1,magnfac);
eldia2(ex2,ey2,es2(:,3),eci2,magnfac);
eldia2(ex3,ey3,es3(:,3),eci3,magnfac,magnitude);
