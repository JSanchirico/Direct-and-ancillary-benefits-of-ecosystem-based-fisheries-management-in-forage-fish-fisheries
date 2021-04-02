
%  Delay Differential fixed step size solver for Sanchirico and Essington, Direct and Ancillary
%  Benefits of EBFM in forage fish fisheries, Ecological Applications,
%  2021. Please cite the main article if you use the code/algorithms.

FLAG = 2;
Param = [0.1 0.3]; % base = .1, .3
Fa = 0.5;
Fs = 2;
Fp = .4;

capT = 1;
tau = 4; % Lag length
tspan = [0 T];


[a, b] = feval(@RK_DDE_FIXEDSTEP,0, tspan(2), y0, tau, .01, tspan(2)/.01, Fa, Fs, Fp, FLAG,VARFLAG, 1, 2, capT, Param, OA_CASE);



function [tvals,yvals]=RK_DDE_FIXEDSTEP(t0,tf,y0,tau,h,n,Fa,Fs,Fp, FLAG,VARFLAG,~,~,capT, Param,OA_CASE)
%*************************************
% t0 = initial time
% tf = terminal time
% y0 = initial conditions
% tau = time delay
% h = step size
% n = steps

peek=10; 
y0=y0(:); % Converts y0 to column vector regardless of what it is

% Calls in the parameter file
[xazero, xszero, a,b,Mpj,tjuv, theta, vbk,kappa,d,wr,Y, alphaa, alphas, Cmax, Mp, amp_a, amp_s, per_a, per_s, start_a, start_s, ...
     recbara, wra, Ma, winfa, kappaa, recbars, wrs, Ms, winfs, kappas, ...
      z, phi, Kbp, theta_1, theta_2, theta_3, Gamma_a, Gamma_s, Gamma_other, Mbp,...
      theta_1s, theta_2s, theta_3s,P,q,c0,c1] = CC_Param_File(FLAG,VARFLAG,capT, Param,OA_CASE);

tc=t0; % initial start
yc=y0; % initial conditions
yp=yc; % Saving for later
tvals=tc;
ypvals=yc;
fc=feval(@derivs,tc,yc); %Starting value of F at the initial conditions
delsteps=round(tau/h); %  this is consistent with  tomlab 


for i=1:delsteps % integrate backward w/o delay to generate IC
    [tc,yp,fc]=RKstep(@derivs,tc,yp,fc,-h);
    ypvals=[y0 ypvals]; 
    tvals=[tc tvals];
end

tc=t0; yp=ypvals; fcp=feval(@derivsdel,tc,yc,yp(:,1),1);
yvals=yp;
% the system at time (t-tau) is the first column of yp
if (n>delsteps)
    yvals=yc;
    tvals=tc;
end

for j=1:n % integrate forward w/ delay
    [tc,yc,fcp]=RKdelstep(@derivsdel,tc,yc,yp,fcp,h,j);
    if mod(j,peek)==0
        yvals=[yvals yc];
        tvals=[tvals tc];
    end
    yp=[yp(:,2:delsteps+1) yc];
end


function [tnew,ynew,fnew]=RKdelstep(fname,tc,yc,yp,fcp,h,j)
ya=yp(:,1);
ya1=0.5*(yp(:,1)+yp(:,2));
ya2=yp(:,2);
k1 = h*fcp;
k2 = h*feval(fname,tc+(h/2),yc+(k1/2),ya1,j);
k3 = h*feval(fname,tc+(h/2),yc+(k2/2),ya1,j);
k4 = h*feval(fname,tc+h,yc+k3,ya2,j);
ynew = yc +(k1 + 2*k2 + 2*k3 +k4)/6;
tnew = tc+h;
fnew = feval(fname,tnew,ynew,yp(:,2),j);
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function Integrates backward at the step-size for the delay
%**************************************
function [tnew,ynew,fnew]=RKstep(fname,tc,yc,fc,h)
k1 = h*fc;
k2 = h*feval(fname,tc+(h/2),yc+(k1/2));
k3 = h*feval(fname,tc+(h/2),yc+(k2/2));
k4 = h*feval(fname,tc+h,yc+k3);
ynew = yc +(k1 + 2*k2 + 2*k3 +k4)/6;
tnew = tc+h;
fnew = feval(fname,tnew,ynew);
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is the equation file with Delay
%**************************************
function dy = derivsdel(tc,yc,ya,j)
% Differential equations function for delaydifferential
% first row is dxa/dt, second is dna/dt, third is dxs/dt, fourth is dns/dt, fifth is dxp/dt, sixth is dnp/dt
xalag = ya(1,1);
%nalag = ya(2,1);
xslag=ya(3,1);
%nalag = ya(4,1);
xplag=ya(5,1);
nbplag = ya(7,1);

xa=yc(1);
na = yc(2);
xs=yc(3);
ns = yc(4);
xp = yc(5);
np=yc(6);
nbp = yc(7);

cc_area = 302000;
t = j*h;

% sardine recruitment
epsst =  exp(amp_s * sin((2 * pi() * t / per_s + start_s * pi() / 180)) - amp_s^2 / 4 ); % this adds the bias correction term to make mean = recbar
recs = recbars * epsst; % sardine recruitment

% anchovy recruitment
epsat =  exp(amp_a * sin((2 * pi() * t / per_a + start_a * pi() / 180)) - amp_a^2 / 4 ); % this adds the bias correction term to make mean = recbar
reca = recbara * epsat; % sardine recruitment

% calculate changes in winf based on consumption rate

C=Cmax*(alphaa*xa+alphas*xs+Y)/(Cmax+alphaa*xa+alphas*xs+Y);
wref=0.0063;

H=C*(wref)^(1-d);
winf=(theta*H/vbk)^(1/(1-d));

% handy to calculate recruitment here:
R=(a*xplag/(1+b*xplag))*exp(-Mpj*tjuv);

Pa = Cmax*alphaa*xp/(Cmax+Y+alphaa*xa+alphas*xs); % predation mortality rate
Ps = Cmax*alphas*xp/(Cmax+Y+alphaa*xa+alphas*xs); % predation mortality rate
 
% brown pelican calculations
% % recruitment
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
smooth=10;

% Ensuring that fishing effort is non-negative in open access model
 e1= max(0,((P(1)*q(1)*xa-c0(1))./c1(1)));  
 e2= max(0,((P(2)*q(2)*xs-c0(2))./c1(2)));
 e3= max(0,((P(3)*q(3)*xp-c0(3))./c1(3)));
  
 haa = q(1)*e1; 
 hss = q(2)*e2; 
 hpp = q(3)*e3; 
 
 
 % Case 9 and 10 are the Catch control rules simulated
 if OA_CASE==9
     Bxs0=3.1788;
     IND=(xs<.4*Bxs0)*1+(xs>.8*Bxs0)*3; IND=(IND==0)*2;
     hss=(IND==1)*0+(IND==2)*((.5*Ms./.4)*xs/Bxs0 - 0.5 * Ms)+(IND==3)*(.5*Ms);
 end
 if OA_CASE==10
     Bcut=150000/cc_area;    
     IND=(xs<=Bcut)*1+(xs>Bcut)*2;
     catch_=(IND==1)*0+(IND==2)*.833*(.87*(xs-Bcut)*(.0465824+.06224328*epsst));
     catch1= catch_*(catch_<200000/cc_area)+(200000/cc_area)*(catch_>200000/cc_area);
     hss=catch1./xs;
 end
 
% These are the ODES of the model 
dy=zeros(7,1);
dy(1) = wra * reca + kappaa * winfa * na - (Ma + haa + Pa) * xa;
dy(2) = reca  - (Ma + haa +Pa) * na;
dy(3) = wrs * recs + kappas * winfs * ns - (Ms + hss + Ps) * xs;
dy(4) = recs  - (Ms + hss +Ps) * ns;
dy(5)= wr*R+kappa*winf*np-(Mp+hpp+kappa)*xp;
dy(6)= R-Mp*np-hpp*np;
dy(7) = rbp - Mbp_t * nbp;
end


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is the history file: Assuming constant
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 function dy=derivs(tc,yc)
 dy=[yc(1);yc(2); yc(3); yc(4); yc(5); yc(6);yc(7)]; % History is the initial conditions
 end
end
