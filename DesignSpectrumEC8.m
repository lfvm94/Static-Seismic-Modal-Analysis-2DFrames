function [Sd,Se]=DesignSpectrumEC8(Z,gamma,beta,damp,T,q,St)
% SYNTAX : [Sd]=DesignSpectrumEC8(Z,gamma,beta,damp,T,q)
%---------------------------------------------------------------------
%    PURPOSE
%     Computes the seismic normalized acceleration from the seismic design 
%     spectrum in accordance to the Eurocode 8 EC8. 
% 
%    INPUT:  Z                  Seismic Hazard Zone - Z1,Z2,Z3
%
%
%            gamma              importance factor
% 
%            beta               
%            damp               viscous damping ratio
%            T                  Structural period corresponding to the mode
%                               of interest
%
%            q                  Behaviour factor
%
%            St:                Soil type
%
%    OUTPUT: Sd :               Design spectrum acceleration
%            Se :               Elastic spectrum acceleration
%
%--------------------------------------------------------------------

% LAST MODIFIED: L.Verduzco    2021-11-14
% Copyright (c)  Faculty of Engineering
%                Autonomous University of Queretaro
%--------------------------------------------------------------------

if Z=='Z1'
    agr=0.16; % m/s^2
elseif Z=='Z2'
    agr=0.24;
elseif Z=='Z3'
    agr=0.35;
else 
    disp('Error. The hazard zone does not exist')
end

ag=agr*gamma;
eta=sqrt(10/(5+damp)); % damping correction factor
if St=='A'
    if eta>=0.55
        S=1.0; % soil coefficient
        Tb=0.15;
        Tc=0.4;
        Td=2.5;
    else
        S=1.0; % soil coefficient
        Tb=0.05;
        Tc=0.25;
        Td=1.2;
    end
elseif St=='B'
    if eta>=0.55
        S=1.2; % soil coefficient
        Tb=0.15;
        Tc=0.5;
        Td=2.5;
        
    else
        S=1.35; % soil coefficient
        Tb=0.05;
        Tc=0.25;
        Td=1.2;
    end
elseif St=='C'
    if eta>=0.55
        S=1.15; % soil coefficient
        Tb=0.2;
        Tc=0.6;
        Td=2.5;
    else
        S=1.5; % soil coefficient
        Tb=0.1;
        Tc=0.25;
        Td=1.2;
    end
elseif St=='D'
    
    if eta>=0.55
        Tb=0.2;
        Tc=0.8;
        Td=2.5;
        S=1.35; % soil coefficient
    else
        Tb=0.1;
        Tc=0.3;
        Td=1.2;
        S=1.8; % soil coefficient
    end
    
elseif St=='E'
    if eta>=0.55
        S=1.4; % soil coefficient
        Tb=0.15;
        Tc=0.8;
        Td=2.5;
    else
        S=1.6; % soil coefficient
        Tb=0.05;
        Tc=0.25;
        Td=1.2;
    end
else 
    disp('Error. That soil type is not supported. Try another.')
end

%% Elastic spectrum acceleration
if 0<=T && T<=Tb
    Se=ag*S*(1+T/Tb*(eta*2.5-1));
elseif Tb<=T && T<=Tc
    Se=ag*S*eta*2.5;
elseif Tc<=T && T<=Td
    Se=ag*S*eta*2.5*(Tc/T);
elseif Td<=T && T<=4*S
    Se=ag*S*eta*2.5*(Tc*Td/T^2);
end

if Se<beta
    Se=beta;
end

Sd=Se/q; % design spectrum acceleration with behaviour factor
