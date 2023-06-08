function [fmaxDOF,Mgl,Kgl,T,La,Egv]=SeismicModalMDOF2DFrames2...
    (coordxy,A,unitWeightEl,qbarxy,Edof,bc,E,I,ni,nf,Psacel,g,modal)
% SYNTAX : 
% [fmaxDOF,Mgl,Kgl,T,La,Egv]=SeismicModalMDOF2DFrames2...
%  (coordxy,A,unitWeightEl,qbarxy,Edof,bc,E,I,ni,nf,Psacel,g,modal)
%---------------------------------------------------------------------
%    PURPOSE
%     To compute the global stiffness matrix of a plane frame as well as
%     the global mass matrix to then call the function:
%     "ModalsMDOF2DFrames2" to determine all the structure's modals.
% 
%    INPUT:  coordxy:           Node coordinates of the structure [x,y]
%
%            A:                 Cross-sectional elements' area
%
%            E,I:               Modulus of Elasticity and Cross-sectional 
%                               inertia of the frame's elements
%
%            bc:                Boundary condition array
%
%            Psacel:            Pseudo-acceleration at the base
%
%            modal:             Mode of vibration of interest:
%                               [mode-1,mode-2,...] -> The equivalent 
%                                    inertial forces are computed with the
%                                    contribution of the given modes of 
%                                    vibration
%                               
%                               mode-i -> The equivalent inertial forces are
%                                       computed with the mode of
%                                       vibration inserted
%
%            g:                 gravity acceleration
%
%            unitWeightEl:      unit weight material of each element:
%                               Size: nbars x 1
%
%            qbarxy:            uniformly distributed loads. 
%                               Size: nbars x 2.
%                               The first column corresponds to the
%                               distributed loads in the local X' direction
%                               and the second column the the loads
%                               distributed in the local Y' direction.
%
%    OUTPUT: La :               Modal of vibration for each DOF. 
%                               Size: Nmodals x 1
%
%            Egv:               DOF's eigenvalues: NDOF x Nmodals
%
%            T :                Structure's periods for each modal      
%
%            fmaxDOF :          Equivalent DOF's forces for the modal
%                               in question
%
%            Mgl:               Global Mass matrix
%            Kgl:               Global Stiffness matrix
%
%--------------------------------------------------------------------

% LAST MODIFIED: L.Verduzco    2023-06-07
% Copyright (c)  Faculty of Engineering
%                Autonomous University of Queretaro
%--------------------------------------------------------------------
nnodes=length(coordxy(:,1)); nbars=length(E);

%% Stiffness and Mass matrices
Kgl=zeros(3*nnodes);
Mgl=zeros(3*nnodes);
for i=1:nbars      
    ex=[coordxy(ni(i),1) coordxy(nf(i),1)];
    ey=[coordxy(ni(i),2) coordxy(nf(i),2)];
    ep=[E(i) A(i) I(i)];
    
    %% Stiffness matrix
    [Kebar]=beam2e(ex,ey,ep);
    
    %% Mass matrix
    PV=unitWeightEl(i,1)-qbarxy(i,2)/A(i); % unit weight of each element
                                          % The distributed downward loads
                                          % on the BEAMS are considered.
    [Mebar]=FiniteMassBeams2D(ex,ey,A(i),PV,g);
    [Kgl]=assem(Edof(i,:),Kgl,Kebar);
    [Mgl]=assem(Edof(i,:),Mgl,Mebar);
end 

%% Modal analysis
[fmaxDOF,T,La,Egv]=ModalsMDOF2DFrames2(Mgl,Kgl,bc,Psacel,modal);

