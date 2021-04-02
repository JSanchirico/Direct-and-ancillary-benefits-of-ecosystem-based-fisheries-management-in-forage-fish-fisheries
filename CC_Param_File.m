
% Model parameters
function varargout=CC_Param_File(flag,VARFLAG, capT,Param,OA_CASE)

% parameters that do not depend on flag    
    cc_area = 302000;
    ra=0.56; % parameters for anchovy
    Ka = 500000 / cc_area; % in mt 
    rs = 0.44; % parameters for sardine
    Ks = 1500000 /cc_area;
    p=.25; % proportion of Cmax for initialization case
    
    winf.star= (7.768E-09 * 1300 ^ 3.0496) / 1000 ; % in mt
    if flag==1 && capT==0
        %Indpendent
        DC_a=0.0001; % proportion of diet that is anchovy at xa = 0.5Ka
        DC_s = 0.00001; % proportion of diet that is sardine at xs = 0.5Ks
        % parameter file for delay differential model
        
    elseif flag==1 && capT==1
        %Indpendent

        DC_a=0.0001; % proportion of diet that is anchovy at xa = 0.5Ka
        DC_s = 0.0001; % proportion of diet that is sardine at xs = 0.5Ks
        % parameter file for delay differential model
       
    
        
    elseif flag==2
        %Bottom-up, no egg predation
  
        DC_a=Param(1); % proportion of diet that is anchovy at xa = recopath levels
        DC_s=Param(2); % proportion of diet that is anchovy at xa = ecopath levels
        % parameter file for delay differential model
       
        
    end
    Fa=0; Fs=0; Fp = 0;

    if VARFLAG ==0
        amp_a = 0; amp_s = 0;
        per_a = 10;
        per_s = 10;
        
        start_a = 0;
        start_s = 90;
    else 
        amp_a = 2;
        amp_s = 2;
    end
    
    
    if VARFLAG ==1
        per_a = 10;
        per_s = 10;
        start_a = 0;
        start_s = 10;
    elseif VARFLAG ==2
        per_a = 10;
        per_s = 10;
        start_a = 0;
        start_s = 180;
    elseif VARFLAG ==3
        per_a = 10;
        per_s = 15;
        start_a = 0;
        start_s = 10;
    end
    
% reproduction parameters (all based on central and southern halibut)
% gives age 0 female recruits
  
a=1.2511E+3;% beverton holt parameter, estimated from stock assessment, assuming h = 0.8, and estimates of bzero.  I used average recruitment to get ro.  Gives female recruits.  Units of spawning biomass are in mt
b=6.0552e-04 * cc_area;% beverton holt parameter, estimated from 

% adult parameters

xazero = 431000 / cc_area;
xszero = 978484 / cc_area;
theta=0.652; % assimilation efficiency, from Holsman, Egestion is 0.2, SDA is 0.118, excretion is .03 of assimilated food so just call id 1 - 0.2 0 0.118 - 0.03 (I know, not strickly correct but close enough) = 0.652.   
Cmax = 9.27; % maximum consumption rate of predator at optimal temp for Halibut, from Holsman et al, assuming M = 0.2 (weighted average of consumption)
Y=(Cmax*p*(DC_a + DC_s - 1))/(p - 1);
alphaa=Cmax*DC_a*p/(xazero*(1-p)); % Correct!
alphas=Cmax*DC_s*p/(xszero*(1-p)); % Correct!  It assumes that xs and xa are at 0.5 ka, ka.
wr=0.0013; % weight (mt) of an age 4 female recruit (from stock assessment parameters)

Mpj=0.2;  % unknown for california halibut.  suggest 0.4, though stock assessment uses 0.2.
tjuv=4; % from stock assessment (roughly age at 50 percent maturation)

d=0.75; % from Essington et al.

wave=0.0063; % average weight, based on M = 0.2, no fishing and growth parameters listed in stock assessment
H=p*Cmax/(wave^(d-1));% # von bertalanffy parameters
vbk=(theta*H/(winf.star^(1-d))); % von bertalanffy parameters
kappa=exp(-vbk)*(1-d); % this is an approximation, seems to work OK but slows down growth a bit
 
Mp=0.2; % From Stock assessment

% parameters for anchovy and sardine model
recbars = 8.75E9/cc_area; % tuned to get biomass near 978K mt
% parameters for new sardine implementation
wrs = .01E-3; % from fitting curve to stock assessment
Ms = 0.6; % from stock assessment
winfs = 0.22E-3; % fit to stock assessment
kappas = 0.215; % fit to stock assessment

% now anchovy
recbara= 1.5E10/ cc_area; % tuned to get biomass near 431K
wra = .01E-3; % no data, assumption
Ma = 0.76; % From method, value where model worked, plus add in predation effect from chub mackerel
winfa = .097E-3;
kappaa=0.13;

% brown pelican parameters
z = 2.39;
phi = 1.1; % chosen arbitrarily, moderate density dependence
Kbp = 100; % arbitrary, can use to set level of acceptable depletion
theta_1 = 0.15; %assumed scenario from Punt et al
theta_2 = 0.2; % base case,  sensitivity analysis is 0;
theta_3 = 0.95;
Gamma_a = 0.644;% ecopath diet composition;
Gamma_s = 0.262;% ecopath diet composition 
Gamma_other = 1-Gamma_a - Gamma_s;
Mbp = 0.22; % 80 percent survival, from Punt et al.
theta_1s = 0.0;
theta_2s = 0.2; %sensitivity analysis is 0.0, 0.2 base case;
theta_3s = 0.95;


%% Economic Parameters for OA runs
% Prices $/mt
P(1)=254.2;  
P(2)=206.3; 
P(3)=7480.0;

% % Catchability coefficient
q=[0.2587    0.2012    0.1725];

% Cost parameters
% Linear 
c0=[25.3256   33.0334   36.7500]; 
% Non-Linear
c1=[4.4850    4.4850    6.5000];


% These cases correspond to what stocks are being fished. Setting q=0
% removes fishing from the pop dynamics

    if OA_CASE==1
         
         q(1)=0; 
       
    elseif OA_CASE==2
        
        q(2)=0;
        
         
    elseif OA_CASE==3
        
        q(3)=0;
        
        
    elseif OA_CASE==4
         
         q(1)=0;
         q(2)=0;
         
    elseif OA_CASE==5
         q(1)=0;
         q(3)=0;
         
    elseif OA_CASE==6
         q(3)=0;
         q(2)=0;
         
     elseif OA_CASE==7
            w=[]; % all are fished
            
    elseif OA_CASE==8 % no fishing on any stock
            q(1)=0;
            q(2)=0;
            q(3)=0;
    else
        w=[];
    
    end

varargout(1:49)={xazero, xszero, a,b,Mpj,tjuv, theta, vbk,kappa,d,wr,Y, alphaa, alphas, Cmax, Mp, amp_a, amp_s, per_a, per_s, start_a, start_s,...
    recbara, wra, Ma, winfa, kappaa, recbars, wrs, Ms, winfs, kappas,...
    z, phi, Kbp, theta_1, theta_2, theta_3, Gamma_a, Gamma_s, Gamma_other, Mbp, ...
    theta_1s, theta_2s, theta_3s, P,q,c0,c1};
