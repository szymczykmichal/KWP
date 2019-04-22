clc;clear;cls;

function [Ke]=bar2e(ex,ey,ep) // Ke=bar2e(ex,ey,ep)
E=ep(1);  A=ep(2);

b=[ ex(2)-ex(1); ey(2)-ey(1) ];
L=sqrt(b'*b);

Kle=E*A/L*[ 1 -1; -1  1];
n=b'/L;
[t1,t2]=size(n)
G=[ n zeros(t1,t2); zeros(t1,t2) n];

//t = [-1 0 1 0]
//Ke=((E*A)/L)*t'*t


Ke=G'*Kle*G;
endfunction

Edof=[1 1 2 5 6;
          2 3 4 7 8;
          3 5 6 9 10;
          4 7 8 11 12;
          5 7 8 5 6;
          6 11 12 9 10;
          7 3 4 5 6;
          8 7 8 9 10;
          9 1 2 7 8;
         10 5 6 11 12];

K=zeros(12,12);
f=zeros(12,1); f(11)=0.5e6*sin(%pi/6); f(12)=-0.5e6*cos(%pi/6);

A=25.0e-4; E=2.1e11; ep=[E A];

ex1=[0 2]; ex2=[0 2]; ex3=[2 4]; ex4=[2 4]; ex5=[2 2];
ex6=[4 4]; ex7=[0 2]; ex8=[2 4]; ex9=[0 2]; ex10=[2 4];
ey1=[2 2]; ey2=[0 0]; ey3=[2 2]; ey4=[0 0]; ey5=[0 2];
ey6=[0 2]; ey7=[0 2]; ey8=[0 2]; ey9=[2 0]; ey10=[2 0];

Ke1=bar2e(ex1,ey1,ep); Ke2=bar2e(ex2,ey2,ep);
Ke3=bar2e(ex3,ey3,ep); Ke4=bar2e(ex4,ey4,ep);
Ke5=bar2e(ex5,ey5,ep); Ke6=bar2e(ex6,ey6,ep);
Ke7=bar2e(ex7,ey7,ep); Ke8=bar2e(ex8,ey8,ep);
Ke9=bar2e(ex9,ey9,ep); Ke10=bar2e(ex10,ey10,ep);



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
K=assem(Edof(3,:),K,Ke3); K=assem(Edof(4,:),K,Ke4);
K=assem(Edof(5,:),K,Ke5); K=assem(Edof(6,:),K,Ke6);
K=assem(Edof(7,:),K,Ke7); K=assem(Edof(8,:),K,Ke8);
K=assem(Edof(9,:),K,Ke9); K=assem(Edof(10,:),K,Ke10);

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
//  pomniejszamy wektor sil wezlowy o wartosci w ezlach gdzie okreslona sa warunki brzegowe
     s=K(fdof,fdof)\(f(fdof)-K(fdof,pdof)*dp);
//
     d(pdof)=dp;
     d(fdof)=s;
  end 
     Q=K*d-f;
   
endfunction

bc=[1 0;2 0;3 0;4 0];

a=solveq(K,f,bc)

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

ed1=extract(Edof(1,:),a); ed2=extract(Edof(2,:),a);
ed3=extract(Edof(3,:),a); ed4=extract(Edof(4,:),a);
ed5=extract(Edof(5,:),a); ed6=extract(Edof(6,:),a);
ed7=extract(Edof(7,:),a); ed8=extract(Edof(8,:),a);
ed9=extract(Edof(9,:),a); ed10=extract(Edof(10,:),a);

function [es]=bar2s(ex,ey,ep,ed)
// es=bar2s(ex,ey,ep,ed)
//-------------------------------------------------------------
// PURPOSE
//  Compute normal force in two dimensional bar element.
//
// INPUT:  ex = [x1 x2]
//         ey = [y1 y2]        element coordinates
//
//         ep = [E A]          E : Young's modulus
//                             A : Cross section area
//
//         ed : [u1 u2 u3 u4]  element displacement vector
//
// OUTPUT: es = [N]            element force 
//-------------------------------------------------------------

// LAST MODIFIED: K Persson    1995-08-23
// Copyright (c)  Division of Structural Mechanics and
//                Department of Solid Mechanics.
//                Lund Institute of Technology
//-------------------------------------------------------------
   E=ep(1);  A=ep(2);
 
   b=[ ex(2)-ex(1); ey(2)-ey(1) ];
   L=sqrt(b'*b);

   Kle=E*A/L*[1 -1 ; -1  1 ];

   n=b'/L;   G=[ n zeros(size(n)); zeros(size(n)) n];
   u=ed';
   N=E*A/L*[-1 1]*G*u;
   es=N;
endfunction

