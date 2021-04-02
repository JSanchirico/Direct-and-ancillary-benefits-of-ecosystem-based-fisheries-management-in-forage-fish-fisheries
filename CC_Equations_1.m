%  Equations for Sanchirico and Essington, Direct and Ancillary
%  Benefits of EBFM in forage fish fisheries, Ecological Applications,
%  2021. Please cite the main article if you use the code/algorithms.

% Timothy E. Essington
% School of Aquatic and Fishery Sciences
% University of Washington
%
% James N. Sanchirico
% Department of Environmental Science and Policy
% University of California Davis
%

% For questions on the code, contact Jim Sanchirico at
% jsanchirico@ucdavis.edu


% Tomlab pieces
smooth=10;


% Non-negativity conditions
cbox_stocks = { 0 <= icollocate(xs)
    0 <= icollocate(xa)
    0 <= icollocate(xh) 
    0 <= icollocate(np)};

% Terminal conditions at T%
cterm = final({xa==INITCOND(1)*1.5;
               xs==INITCOND(3)*1.5
               xh==INITCOND(5)*1.5});


% Generating initial condition
cbnd = initial({....
    xa == INITCOND(1)
    na == INITCOND(2)
    xs == INITCOND(3)
    ns == INITCOND(4)
    xh == INITCOND(5)
    nh == INITCOND(6)
    np == INITCOND(7)
    });


% Time lags in population dynamics
xalag = ifThenElse(t<tau, INITCOND(1), subs(xa,t,t-tau));
xslag = ifThenElse(t<tau, INITCOND(3), subs(xs,t,t-tau));
xplag = ifThenElse(t<tau, INITCOND(5), subs(xh,t,t-tau));
nbplag = ifThenElse(t<tau, INITCOND(7), subs(np,t,t-tau));


%% Calling the equations for model

% sardine recruitment
epsst =  exp(amp_s * sin((2 * pi() * t / per_s + start_s * pi() / 180)) - amp_s^2 / 4 ); % this adds the bias correction term to make mean = recbar
recs = recbars * epsst; % sardine recruitment

% anchovy recruitment
epsat =  exp(amp_a * sin((2 * pi() * t / per_a + start_a * pi() / 180)) - amp_a^2 / 4 ); % this adds the bias correction term to make mean = recbar
reca = recbara * epsat; % anchovy recruitment

% calculate changes in winf based on consumption rate

C=Cmax*(alphaa*xa+alphas*xs+Y)/(Cmax+alphaa*xa+alphas*xs+Y);
wref=0.0063;

H=C*(wref)^(1-d);
winf=(theta*H/vbk)^(1/(1-d));

% handy to calculate recruitment here:
R=(a*xplag/(1+b*xplag))*exp(-Mpj*tjuv);

Pa = Cmax*alphaa*xh/(Cmax+Y+alphaa*xa+alphas*xs); % predation mortality rate
Ps = Cmax*alphas*xh/(Cmax+Y+alphaa*xa+alphas*xs); % predation mortality rate

% brown pelican calculations
% recruitment
d_tlag = Gamma_a * xalag / xazero + Gamma_s * xslag / xszero + Gamma_other; % relative prey abundance
psi_r = max(0, (1-theta_1 - theta_2) * theta_3 * (d_tlag - theta_1)) / ((1-theta_1)*theta_2 * (1-theta_3) + (theta_3* (1-theta_1) - theta_2) * ( d_tlag - theta_1));

% Survival
% depletion based on current biomass
d_t = Gamma_a * xa / xazero + Gamma_s * xs / xszero + Gamma_other; % relative prey abundance
psi_s = max(0, (1-theta_1s - theta_2s) * theta_3s * (d_t - theta_1s)) / ((1-theta_1s)*theta_2s * (1-theta_3s) + (theta_3s* (1-theta_1s) - theta_2s) * ( d_t - theta_1s));
psi_slag = max(0, (1-theta_1s - theta_2s) * theta_3s * (d_tlag - theta_1s)) / ((1-theta_1s)*theta_2s * (1-theta_3s) + (theta_3s* (1-theta_1s) - theta_2s) * ( d_tlag - theta_1s));

S_t = exp(-Mbp) * psi_s;
Mbp_t = - log(S_t);
S_tlag = exp(-Mbp) * psi_slag;
Mbp_tlag = -log(S_tlag);

rbp = 0.52 * exp(-4 * Mbp_tlag) * nbplag * psi_r * max(0,(1 + (phi - 1)*(1-(psi_r*nbplag/Kbp)^z)));

% difference cases with some fisheries OA and some optimal
if CASE==1
    cbox_case=[]; %cbox_sh_oa; % A optimal
   cbox_F_pos = {...
    0 <= q(1)*icollocate(e1)
     };
    
    
    e3= max(0,((P(3)*q(3)*xh-c0(3))./c1(3))); %e3=(e3>0)*e3;
    e2= max(0,((P(2)*q(2)*xs-c0(2))./c1(2))); %e2=(e2>0)*e2;
    
elseif CASE==2
    cbox_case=[]; %cbox_ah_oa; % S optimal
    cbox_F_pos = {...
    
    0 <= q(2)*icollocate(e2)
   };
    
    
    e1= max(0,((P(1)*q(1)*xa-c0(1))./c1(1))); % e1=(e1>0)*e1;
    e3= max(0,((P(3)*q(3)*xh-c0(3))./c1(3))); %e3=(e3>0)*e3; e3= ((P(3)*q(3)*xh-c0(3))./c1(3)); e3=(e3>0)*e3;
    
elseif CASE==3
    cbox_case=[]; %cbox_as_oa; % H optimal
    cbox_F_pos = {...
    
    0 <= q(3)*icollocate(e3) };
    

   e2= max(0,((P(2)*q(2)*xs-c0(2))./c1(2))); %e2=(e2>0)*e2;
   e1= max(0,((P(1)*q(1)*xa-c0(1))./c1(1))); % e1=(e1>0)*e1;
    
elseif CASE==4
    cbox_case=[];%cbox_h_oa; % A+S optimal
     cbox_F_pos = {...
    0 <= q(1)*icollocate(e1)
    0 <= q(2)*icollocate(e2)
     };

    e3= max(0,((P(3)*q(3)*xh-c0(3))./c1(3))); %e3=(e3>0)*e3;
    
elseif CASE==5
    cbox_case=[]; %cbox_s_oa; % A+H optimal
     cbox_F_pos = {...
    0 <= q(1)*icollocate(e1)
    
    0 <= q(3)*icollocate(e3) };

    
      e2= max(0,((P(2)*q(2)*xs-c0(2))./c1(2))); %e2=(e2>0)*e2;
   
elseif CASE==6
    cbox_case=[]; %cbox_a_oa; % S+H optimal
     cbox_F_pos = {...
    
    0 <= q(2)*icollocate(e2)
    0 <= q(3)*icollocate(e3) };

      e1= max(0,((P(1)*q(1)*xa-c0(1))./c1(1))); % e1=(e1>0)*e1;
else
   cbox_case=[]; 
    cbox_F_pos = {...
    0 <= q(1)*icollocate(e1)
    0 <= q(2)*icollocate(e2)
    0 <= q(3)*icollocate(e3) };
    
end


ceq = collocate({...
        dot(xa) == wra * reca + kappaa * winfa * na - (Ma + q(1)*e1 + Pa) * xa;...
        dot(xs) == wrs * recs + kappas * winfs * ns - (Ms + q(2)*e2 + Ps) * xs;...
        dot(na) == reca  - (Ma + q(1)*e1 +Pa) * na;...
        dot(ns) == recs  - (Ms + q(2)*e2 +Ps) * ns;...
        dot(xh) == wr*R+kappa*winf*nh-(Mp+q(3)*e3+kappa)*xh;...
        dot(nh) == R-Mp*nh-q(3)*e3*nh;...
        dot(np) == rbp - Mbp_t * np});


BETA = [1; 1; 1]; % not used but permits a more flexibile catch function

if CASE==1
    
    profit=P(1)*q(1)*e1*(xa^BETA(1))-c0(1)*(e1)-c1(1)*(e1^2); % A optimal
    
elseif CASE==2
    profit=P(2)*q(2)*e2*(xs^BETA(2))-c0(2)*(e2)-c1(2)*(e2^2); % S optimal
elseif CASE==3
    profit=P(3)*q(3)*e3*(xh^BETA(3))-c0(3)*(e3)-c1(3)*(e3^2); % H optimal
elseif CASE==4 % A+S optimal
    profit  = P(1)*q(1)*e1*(xa^BETA(1))-c0(1)*(e1)-c1(1)*(e1^2)+...
        P(2)*q(2)*e2*(xs^BETA(2))-c0(2)*(e2)-c1(2)*(e2^2);
    
elseif CASE==5% A+H optimal
    profit  = P(1)*q(1)*e1*(xa^BETA(1))-c0(1)*(e1)-c1(1)*(e1^2)+...
        P(3)*q(3)*e3*(xh^BETA(3))-c0(3)*(e3)-c1(3)*(e3^2);
    
elseif CASE==6 % S+H optimal
    profit  = P(2)*q(2)*e2*(xs^BETA(2))-c0(2)*(e2)-c1(2)*(e2^2)+...
        P(3)*q(3)*e3*(xh^BETA(3))-c0(3)*(e3)-c1(3)*(e3^2);
else % all optimal
    
    
    profit  = P(1)*q(1)*e1*(xa^BETA(1))-c0(1)*(e1)-c1(1)*(e1^2)+...
        P(2)*q(2)*e2*(xs^BETA(2))-c0(2)*(e2)-c1(2)*(e2^2)+...
        P(3)*q(3)*e3*(xh^BETA(3))-c0(3)*(e3)-c1(3)*(e3^2);
end


% Objective function
objective = -integrate(exp(-dis*t).*profit);



