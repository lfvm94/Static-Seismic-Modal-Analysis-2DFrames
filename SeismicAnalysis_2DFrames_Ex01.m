% SeismicAnalysis_2DFrames_Ex01
%----------------------------------------------------------------
% PURPOSE 
%    To compute the Static Modal analysis for a Reinforced Concrete
%    Plane Frame.
%
%----------------------------------------------------------------

% LAST MODIFIED: L.F.Veduzco    2023-06-07
%                Faculty of Engineering
%                Autonomous University of Queretaro
%----------------------------------------------------------------

clc 
clear all

nnodes=12;
nbars=15;

%% Materials
fpc=[300;
     300;
     300;
     300;
     300;
     300;
     300;
     300;
     300;
     300;
     300;
     300;
     300;
     300;
     300];

% Modulus of Elasticity of each element
E=14000*(fpc).^0.5;

%% Geometry
dimensions=[20 50;
            20 50;
            20 50;
            20 50;
            20 50;
            20 50;
            20 50;
            20 50;
            20 50;
            20 50;
            20 50;
            20 50;
            20 50;
            20 50;
            20 50];

A=dimensions(:,1).*dimensions(:,2);
I=1/12.*dimensions(:,1).*dimensions(:,2).^3;

% Coordinates of each node for each bar
coordxy=[0 0;
         0 300;
         0 600;
         0 900;
         600 0;
         600 300;
         600 600;
         600 900;
         1200 0;
         1200 300;
         1200 600;
         1200 900]; 
                 
%% Topology
% Node conectivity
ni=[1;2;3;4;3;2;5;6;7;8; 7; 6; 9; 10;11];
nf=[2;3;4;8;7;6;6;7;8;12;11;10;10;11;12];

% Length of each element
L=sqrt((coordxy(nf,1)-coordxy(ni,1)).^2+(coordxy(nf,2)-coordxy(ni,2)).^2);

% Topology matrix
Edof=zeros(nbars,7);
for i=1:nbars
    Edof(i,1)=i;
    Edof(i,2)=ni(i)*3-2;
    Edof(i,3)=ni(i)*3-1;
    Edof(i,4)=ni(i)*3;
    
    Edof(i,5)=nf(i)*3-2;
    Edof(i,6)=nf(i)*3-1;
    Edof(i,7)=nf(i)*3;
end

%% Boundary conditions
% Prescribed DOF - [N-DOF, displacement]
bc=[1 0;
    2 0;
    3 0;
    13 0;
    14 0;
    15 0;
    25 0;
    26 0;
    27 0];

%% Loads

type_elem=[1 "Col";
           2 "Col";
           3 "Col";
           4 "Beam";
           5 "Beam";
           6 "Beam";
           7 "Col";
           8 "Col";
           9 "Col";
           10 "Beam";
           11 "Beam";
           12 "Beam";
           13 "Col";
           14 "Col";
           15 "Col"];
       
% Distributed Loads on beams
beamsLoads=[1 -50;
          2 -50;
          3 -50;
          4 -50;
          5 -50;
          6 -50];

elem_cols=[];
elem_beams=[];
beams=0;
cols=0;
for j=1:nbars % To identify which elements are beams and which are columns
    if type_elem(j,2)=="Beam"
        beams=beams+1;
        elem_beams=[elem_beams,j];
    elseif type_elem(j,2)=="Col"
        cols=cols+1;
        elem_cols=[elem_cols,j];
    end
end
qbarxy=zeros(nbars,2);
for i=1:beams
    qbarxy(elem_beams(i),2)=1.1*(beamsLoads(i,2));
end

modal=1;

%% Seismic response spectrum

g=981; % gravity acceleration

%% Modal analysis
pvconc=0.0024; % unit weight of concrete
unitWeightElm=zeros(nbars,1)+pvconc;

% Consistent mass method
DS=1;
[Sd,fmaxDOF,Mglobal,Kglobal,T,La,Egv,Ma]=SeismicModalMDOF2DFrames2...
(coordxy,A,unitWeightElm,qbarxy,Edof,bc,E,I,ni,nf,DS,g,modal);

% Considering the equivalent seismic loads for a structural analysis
fglobal=fmaxDOF;

%% Static structural analysis with seismic forces
np=7; % number of analysis points for the mechanical elements

[Ugl1,reactions,Ex,Ey,esbarsnormal,esbarsshear,esbarsmoment]=...
StaticLinearAnalysis2DFrame(E,A,I,bc,fglobal,ni,nf,qbarxy,np,coordxy,...
1,1,6,6);

%% Plot of the modal in question and its frequency
Freq=1./T;

if length(modal)==1 % If only one modal was entered
    figure(6)
    grid on
    NoteMode=num2str(modal);
    title(strcat('Eigenmode ','- ',NoteMode))
    eldraw2(Ex,Ey,[2 3 1]); 
    Edb=extract(Edof,Egv(:,modal));
    eldisp2(Ex,Ey,Edb,[1 2 2]);
    FreqText=num2str(Freq(modal));
    NotaFreq=strcat('Freq(Hz)= ',FreqText);
    text(50,10,NotaFreq);
end
%% Plot of the seismic structural reponse for the first 10 modes
widthstruc=max(coordxy(:,1));
heightstruc=max(coordxy(:,2));

figure(7)
axis('equal')
axis off
title(strcat('Deformed structures against the seismic actions for the',...
    ' first 10 EigenModes'))
hold on

% First 5 modals at the top of the plot
for i=1:5
    
    % Modal static analysis - Consistent mass method
    [Sd,fmaxDOF,Mglobal,Kglobal,T,La,Egv,Ma]=SeismicModalMDOF2DFrames2...
    (coordxy,A,unitWeightElm,qbarxy,Edof,bc,E,I,ni,nf,DS,g,i);

    % Static structural analysis with seismic forces
    [Ugl,reactions,Ex,Ey,esbarsnormal,esbarsshear,esbarsmoment]=...
    StaticLinearAnalysis2DFrame(E,A,I,bc,fmaxDOF,ni,nf,...
    qbarxy,np,coordxy,0);

    Ext=Ex+(i-1)*(widthstruc+150);
    eldraw2(Ext,Ey,[2 3 1]);
    Edb=extract(Edof,Egv(:,i));
    eldisp2(Ext,Ey,Edb,[1 2 2]);
    FreqText=num2str(Freq(i));
    NotaFreq=strcat('Frec(Hz)= ',FreqText);

    text((widthstruc+150)*(i-1)+50,150,NotaFreq)
end

% Last 5 modals at the bottom of the plot
Eyt=Ey-(heightstruc+200);
for i=6:10
    % Modal static analysis - Consistent mass method
    [Sd,fmaxDOF,Mglobal,Kglobal,T,La,Egv,Ma]=SeismicModalMDOF2DFrames2...
    (coordxy,A,unitWeightElm,qbarxy,Edof,bc,E,I,ni,nf,DS,g,i);

    %% Static structural analysis with seismic forces
    [Ugl,reactions,Ex,Ey,esbarsnormal,esbarsshear,esbarsmoment]=...
    StaticLinearAnalysis2DFrame(E,A,I,bc,fmaxDOF,ni,nf,...
    qbarxy,np,coordxy,0);

    Ext=Ex+(i-6)*(widthstruc+150);
    eldraw2(Ext,Eyt,[2 3 1]);
    Edb=extract(Edof,Egv(:,i));
    eldisp2(Ext,Eyt,Edb,[1 2 2]);
    FreqText=num2str(Freq(i));
    NotaFreq=strcat('Frec(Hz)= ',FreqText);
    
    text((widthstruc+150)*(i-6)+50,-heightstruc,NotaFreq);
end