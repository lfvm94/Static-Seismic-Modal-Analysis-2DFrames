function [Me]=FiniteMassBeams2D(ex,ey,A,PV,g)

% SYNTAX : [Me]=FiniteMassBeams2D(ex,ey,A,PV,g)
%---------------------------------------------------------------------
%    PURPOSE
%     Compute the mass matrix for a two dimensional beam element. 
% 
%    INPUT:  ex = [x1 x2]
%            ey = [y1 y2]       element node coordinates
%
%
%            eq = [qx qy]       distributed loads, local directions
% 
%            A                  Transversal area of element
%            PV                 Volumetric weigth of material element
%            g                  Gravity acceleration

%    OUTPUT: Me : element mass matrix (6 x 6)
%
%--------------------------------------------------------------------

% LAST MODIFIED: L.Verduzco    2021-11-14
% Copyright (c)  Faculty of Engineering
%                Autonomous University of Queretaro
%--------------------------------------------------------------------

b=[ ex(2)-ex(1); ey(2)-ey(1) ];
L=sqrt(b'*b);  n=b/L;

Mle=PV*A*L/(420*g)*[140  0     0      70     0           0;
                    0   156    22*L   0      54      -13*L;
                    0   22*L   4*L^2  0      13*L   -3*L^2;
                    70   0     0      140    0           0;
                    0   54     13*L   0      156     -22*L;
                    0   -13*L  -3*L^2 0     -22*L    4*L^2];
        
G=[n(1) n(2)  0    0    0   0;
    -n(2) n(1)  0    0    0   0;
      0    0    1    0    0   0;
      0    0    0   n(1) n(2) 0;
      0    0    0  -n(2) n(1) 0;
      0    0    0    0    0   1];
  
  Me=G'*Mle*G;