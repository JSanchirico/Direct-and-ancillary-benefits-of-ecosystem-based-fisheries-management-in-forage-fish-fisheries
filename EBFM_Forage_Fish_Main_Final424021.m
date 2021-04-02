%  Main program file for Sanchirico and Essington, Direct and Ancillary
%  Benefits of EBFM in forage fish fisheries, Ecological Applications,
%  2021. Please cite the main article if you use the code/algorithms.

%  Please note: Code is developed to generate a solution in the paper for
%  each management scenario at a base case set of parameters. This code does not
%  recreate each figure in the paper or do the sensitivity analysis.

%  Information on how to recreate the runs is available in the paper and
%  SI. 

%  To run this code, you need to download software from Tomlab
%  tomopt.com/tomlab/ (propt and main solver package). While this is a
%  licensed product, you are able to download a free trial.

% For questions on the code, contact Jim Sanchirico at
% jsanchirico@ucdavis.edu

% To recreate the results in the paper, you will need the accompaning files:
% (1) Main program file - this file.
% (2) Parameter file
% (3) Equations file
% (4) Profit calculation
% (5) RK fixed step solver (adapted for our program)
% (6) Analysis code
% (7) Data used for intitial guesses

% All of these files are in the github repository.


% Copyrighted by Jim Sanchirico


%% *Main program file*
clear all
close all


% General parameters

% Key for parameters: 1 = Anchovy; 2=Sardine; 3=halibut; 4=pelican

% Scaling parameter for CA Current
cc_area = 302000;

% Discount rate
dis = 0.04;

% This parameter permits a more general production function than used in the paper
BETA=[1; 1; 1]; % Set to one in the paper


% Dynamic model parameters

% Time period of analysis
t0=0;  T=40; 

% Vector of #'s of collocation nodes 
NN=[5 10 45 55 75 90 100 110 120 130 140];

% Parameter used in the sensitivity analysis (diet parameters)
Param = [0.1 0.3]; 

% Anchovy biomass #, Sardine bio #, halibut bio #, pelican #
INITCOND = [400000 1000E3  750000 5000E3 1.5*11800 1.5*2.0382E6]./cc_area; 
INITCOND= [INITCOND 100]; % adding pelicans in without scaling them
y0=INITCOND; % used in the DDE code


% Main dynamic analysis

% Code for tomlab solvers -- please see Tomlab instruction manual for
% information
flag=2;
capT=1; % indicator variable that was for whether imposing the T conditions
%

% These need to be filled out for your computer. I ran the code before
% removing them for my computer and code works. 

STR=' '; % Set this for your computer
PATH=[STR ' ']; % Set this for your computer


data_params=[flag,capT,cc_area,dis,INITCOND,Param,BETA(1),BETA(2),BETA(3)];
save (strcat(PATH,'Parameters'),'data_params')

% The following is necessary when imposing the Pelican constraint on the
% optimal solution, as the constraint needs to be relevant for any given
% set of parameters. 

% Used in base case and econo sensitivity
           PelicanAVG=[95.2287;75.3698;88.332;80.415];
            
% Used in survival sensitivity case
           % PelicanAVG= [98.0313;86.76 ;90.59;88.1402];
            
% Used in recruitment sensitivity case
            %PelicanAVG= [98.0313; 67.72;86.58;75.14];


%% Running the OA model results across all four VARFLAGS, save the dynamics and save the P(T) level

% scenarios
% VARFLAG = 0; % no variability
% VARFLAG = 1; % in phase, same frequency
% VARFLAG = 2; % out of phase, same frequency
% VARFLAG = 3; % different period, same starting phase


%T=200; 
Bxs0=3.1788; %determined by running out 200 years and taking average between 100-200
for Varflag=1:4
    VARFLAG=Varflag-1;
    
    
    
    for OA_CASE=1:10  % Cycles thru the different cases of open-access fishing, no fishing, and the catch control rules.
        
        [xazero, xszero, a,b,Mpj,tjuv, theta, vbk,kappa,d,wr,Y, alphaa, alphas, Cmax, Mp, amp_a, amp_s, per_a, per_s, start_a, start_s,...
            recbara, wra, Ma, winfa, kappaa, recbars, wrs, Ms, winfs, kappas,...
            z, phi, Kbp, theta_1, theta_2, theta_3, Gamma_a, Gamma_s, Gamma_other, Mbp, ...
            theta_1s, theta_2s, theta_3s, P,q,c0,c1]=CC_Param_File(flag, VARFLAG, capT,Param,OA_CASE);
        % Calls the DDE solver for each case
        RK_DDE_FIXEDSTEP_V4JNS2
        
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
            w=[];
            
    elseif OA_CASE==8
            q(1)=0;
            q(2)=0;
            q(3)=0;
    else
        w=[];
    
     end
        
    if OA_CASE<9 % Ensures that fishing effort remains non-negative in OA runs
        
        e1oa= max(0,((P(1)*q(1)*b(1,:)-c0(1))./c1(1)));
        e2oa= max(0,((P(2)*q(2)*b(3,:)-c0(2))./c1(2)));
        e3oa= max(0,((P(3)*q(3)*b(5,:)-c0(3))./c1(3)));
          
    elseif OA_CASE==9  % For catch control rule for Sardine
        
        e1oa= max(0,((P(1)*q(1)*b(1,:)-c0(1))./c1(1))); 
        
        for i=1:length(b(3,:))
            IND=0;
            IND=(b(3,i)<.4*Bxs0)*1+(b(3,i)>.8*Bxs0)*3; IND=(IND==0)*2;
            hss=(IND==1)*0+(IND==2)*((.5*Ms./.4))*(b(3,i)./Bxs0 - .5*Ms)+(IND==3)*(.5*Ms);  
            e2oa(1,i)=hss./q(2);
            hcr_rent(i)=(P(2)*q(2)*b(3,i)*e2oa(1,i)-c0(2)*e2oa(1,i)-c1(2)*e2oa(1,i)^2);
        end
        
        e3oa= max(0,((P(3)*q(3)*b(5,:)-c0(3))./c1(3)));
        
    elseif OA_CASE==10 % For catch control rule for Sardine
        
        e1oa= max(0,((P(1)*q(1)*b(1,:)-c0(1))./c1(1))); 
        
        for i=1:length(b(3,:))
            epsst =  exp(amp_s * sin((2 * pi() * i / per_s + start_s * pi() / 180)) - amp_s^2 / 4 );
            IND=0;
            Bcut=150000/cc_area;
            IND=(b(3,i)<=Bcut)*1+(b(3,i)>Bcut)*2;
            catch_=(IND==1)*0+(IND==2)*.833*(.87*(b(3,i)-Bcut)*(.0465824+.06224328*epsst));
            catch1= catch_*(catch_<200000/cc_area)+(200000/cc_area)*(catch_>200000/cc_area);
            hss=catch1./b(3,i);            
            e2oa(1,i)=hss./q(2);
            hcr2_rent(i)=(P(2)*q(2)*b(3,i)*e2oa(1,i)-c0(2)*e2oa(1,i)-c1(2)*e2oa(1,i)^2);
            
        end
        
        e3oa= max(0,((P(3)*q(3)*b(5,:)-c0(3))./c1(3)));
    end
    
        OA_DYNAMICS{Varflag,OA_CASE}=[b; e1oa; e2oa; e3oa];
        OA_pelicanT{Varflag,OA_CASE}=b(7,end);  % Saves the Terminal Value of Pelican Population
        
    end
end

OA_CASE=7; % Resets q's for the analysis below

% save the output for graphing
save (strcat(PATH,'openaccess_pelicanT'), 'OA_pelicanT')
save (strcat(PATH,'openaccess_dynamics'), 'OA_DYNAMICS')


%% Running thru the different modeling cases : Imposing constraints on the system for the different runs


% Dynamic parameters
%VARFLAG = 0; % no variability
%VARFLAG = 1; % in phase, same frequency
%VARFLAG = 2; % out of phase, same frequency
%VARFLAG = 3; % different period, same starting phase; Cheb method works better for this case....
i=1;
 

% Loading initial guess -- this file is provided to you.
load (strcat(PATH,'InitialGuessData'))
GUESS=SolutionSavesALL;
clear SolutionSavesALL;

for PELICAN=1:5
    
    for Varflag=1:4
        
        VARFLAG=Varflag-1;
        
        for CASE=[1 2 4 5 6 7]  % Taking out case 3 from the analysis Halibut only
            
            toms t
            
            % Lag is on
            tau=4;
            options=struct;
            
            if VARFLAG<3
                Phas=tomPhase('Phas',t,t0,T,NN(8));
                options.solver='snopt'; %,[],'Cheb'); %,'fem1s')% 'fem1') %,[],'Cheb');
            else
                
                Phas=tomPhase('Phas',t,t0,T,NN(8),[],'Cheb');
                options.solver='snopt';
                
            end
            
         
            
            setPhase(Phas)
            
            if PELICAN>1
                Efix1=.01; Efix2=.01; Efix3=.01;
            else
                Efix1=.15; Efix2=.15; Efix3=.15;
            end
           
             
            if PELICAN>1
                
                if CASE<7
                    YY(CASE)=(OA_pelicanT{Varflag,CASE}-Opt_pelicanT{Varflag,CASE});
                else
                    YY(CASE)=(OA_pelicanT{Varflag,8}-Opt_pelicanT{Varflag,CASE});
                end
                
                if YY(CASE)>0
                    if PELICAN==2
                        ConsP=.05*YY(CASE)+Opt_pelicanT{Varflag,CASE};
                    elseif PELICAN==3
                        ConsP=.1*YY(CASE)+Opt_pelicanT{Varflag,CASE};
                    elseif PELICAN==4
                        ConsP=.15*YY(CASE)+Opt_pelicanT{Varflag,CASE};
                    elseif PELICAN==5
                        ConsP=.2*YY(CASE)+Opt_pelicanT{Varflag,CASE};
                    end
                elseif YY(CASE)<0&&CASE~=3
                    if PELICAN==2
                        ConsP=.05*(100-Opt_pelicanT{Varflag,CASE})+Opt_pelicanT{Varflag,CASE};
                    elseif PELICAN==3
                        ConsP=.1*(100-Opt_pelicanT{Varflag,CASE})+Opt_pelicanT{Varflag,CASE};
                    elseif PELICAN==4
                        ConsP=.15*(100-Opt_pelicanT{Varflag,CASE})+Opt_pelicanT{Varflag,CASE};
                    elseif PELICAN==5
                        ConsP=.2*(100-Opt_pelicanT{Varflag,CASE})+Opt_pelicanT{Varflag,CASE};
                    end
                    
                elseif CASE==3
                    if PELICAN==2
                        ConsP=Opt_pelicanT{Varflag,CASE};
                    elseif PELICAN==3
                        ConsP=Opt_pelicanT{Varflag,CASE};
                    elseif PELICAN==4
                        ConsP=Opt_pelicanT{Varflag,CASE};
                    elseif PELICAN==5
                        ConsP=Opt_pelicanT{Varflag,CASE};
                    end
                end
            end
            
            tomStates xs ns xa na xh nh np
            
            
            % Initial guesses after a run is finished
            x = GUESS{Varflag,PELICAN,CASE,i};
            clear x01
            
            
            [xazero, xszero, a,b,Mpj,tjuv, theta, vbk,kappa,d,wr,Y, alphaa, alphas, Cmax, Mp, amp_a, amp_s, per_a, per_s, start_a, start_s,...
                recbara, wra, Ma, winfa, kappaa, recbars, wrs, Ms, winfs, kappas,...
                z, phi, Kbp, theta_1, theta_2, theta_3, Gamma_a, Gamma_s, Gamma_other, Mbp, ...
                theta_1s, theta_2s, theta_3s, P,q,c0,c1]=CC_Param_File(flag, VARFLAG, capT,Param,7);
            
            
            if CASE==1
                
                tomControls e1
                
                % initial guess
                x01 = { icollocate({...
                    xa == INITCOND(1)/2
                    na == INITCOND(2)/2
                    xs == INITCOND(3)/2
                    ns == INITCOND(4)/2
                    xh == INITCOND(5)/2
                    nh == INITCOND(6)/2
                    np == INITCOND(7)/2
                    })
                    collocate(e1==.1)
                    };
                
                if tau>0
                    
                    x01 = {...
                        icollocate(xa == subs(xa,x))
                        icollocate(na == subs(na,x))
                        icollocate(xs == subs(xs,x))
                        icollocate(ns == subs(ns,x))
                        icollocate(xh == subs(xh,x))
                        icollocate(nh == subs(nh,x))
                        icollocate(np == subs(np,x))
                        collocate(e1==subs(e1,x))
                        };
                end
                
                cbox_E_free={ 0.0 <= q(1)*collocate(e1)
                    };
                
                cbox_E_fix={ Efix1 == icollocate(e1)};
                
                
            elseif CASE==2
                
                
                tomControls e2
                % initial guess
                x01 = { icollocate({...
                    xa == INITCOND(1)/2
                    na == INITCOND(2)/2
                    xs == INITCOND(3)/2
                    ns == INITCOND(4)/2
                    xh == INITCOND(5)/2
                    nh == INITCOND(6)/2
                    np == INITCOND(7)/2
                    })
                    collocate(e2==.2)
                    };
                
                cbox_E_free={
                    0.0 <= q(2)*collocate(e2)
                    };
                
                cbox_E_fix={
                    Efix2 == icollocate(e2)
                    };
                
                if tau>0
                    
                    x01 = { icollocate({...
                        xa == subs(xa,x)
                        na == subs(na,x)
                        xs == subs(xs,x)
                        ns == subs(ns,x)
                        xh == subs(xh,x)
                        nh == subs(nh,x)
                        np == subs(np,x)
                        })
                        collocate({
                        e2==subs(e2,x)
                        })
                        };
                end
                
            elseif CASE==3
                
                tomControls e3
                % initial guess
                x01 = { icollocate({...
                    xa == INITCOND(1)/2
                    na == INITCOND(2)/2
                    xs == INITCOND(3)/2
                    ns == INITCOND(4)/2
                    xh == INITCOND(5)/2
                    nh == INITCOND(6)/2
                    np == INITCOND(7)/2
                    })
                    collocate(e3==.2)};
                
                cbox_E_free={
                    0.0 <= q(3)*collocate(e3)};
                
                cbox_E_fix={Efix3 == icollocate(e3)};
                
                if tau>0
                    
                    x01 = { icollocate({...
                        xa == subs(xa,x)
                        na == subs(na,x)
                        xs == subs(xs,x)
                        ns == subs(ns,x)
                        xh == subs(xh,x)
                        nh == subs(nh,x)
                        np == subs(np,x)
                        })
                        collocate({e3==subs(e3,x)
                        })
                        };
                end
                
            elseif CASE==4
                tomControls e1 e2
                % initial guess
                x01 = { icollocate({...
                    xa == INITCOND(1)/2
                    na == INITCOND(2)/2
                    xs == INITCOND(3)/2
                    ns == INITCOND(4)/2
                    xh == INITCOND(5)/2
                    nh == INITCOND(6)/2
                    np == INITCOND(7)/2
                    })
                    collocate(e1==.1)
                    collocate(e2==.2)
                    };
                
                cbox_E_free={ 0.0 <= q(1)*collocate(e1)
                    0.0 <= q(2)*collocate(e2)
                    };
                
                cbox_E_fix={ Efix1 == icollocate(e1)
                    Efix2 == icollocate(e2)};
                
                if tau>0
                    
                    x01 = { icollocate({...
                        xa == subs(xa,x)
                        na == subs(na,x)
                        xs == subs(xs,x)
                        ns == subs(ns,x)
                        xh == subs(xh,x)
                        nh == subs(nh,x)
                        np == subs(np,x)
                        })
                        collocate({
                        e2==subs(e2,x)
                        e1==subs(e1,x)})
                        };
                end
                
                
                
            elseif CASE==5
                tomControls e1 e3
                % initial guess
                x01 = { icollocate({...
                    xa == INITCOND(1)/2
                    na == INITCOND(2)/2
                    xs == INITCOND(3)/2
                    ns == INITCOND(4)/2
                    xh == INITCOND(5)/2
                    nh == INITCOND(6)/2
                    np == INITCOND(7)/2
                    })
                    collocate(e1==.1)
                    collocate(e3==.2)};
                
                cbox_E_free={ 0.0 <= q(1)*collocate(e1)
                    0.0 <= q(3)*collocate(e3)};
                
                cbox_E_fix={ Efix1 == icollocate(e1)
                    Efix3 == icollocate(e3)};
                
                if tau>0
                    
                    x01 = { icollocate({...
                        xa == subs(xa,x)
                        na == subs(na,x)
                        xs == subs(xs,x)
                        ns == subs(ns,x)
                        xh == subs(xh,x)
                        nh == subs(nh,x)
                        np == subs(np,x)
                        })
                        collocate({e3==subs(e3,x)
                        
                        e1==subs(e1,x)})
                        };
                end
                
            elseif CASE==6
                tomControls e2 e3
                % initial guess
                x01 = { icollocate({...
                    xa == INITCOND(1)/2
                    na == INITCOND(2)/2
                    xs == INITCOND(3)/2
                    ns == INITCOND(4)/2
                    xh == INITCOND(5)/2
                    nh == INITCOND(6)/2
                    np == INITCOND(7)/2
                    })
                    collocate(e2==.2)
                    collocate(e3==.2)};
                
                cbox_E_free={
                    0.0 <= q(2)*collocate(e2)
                    0.0 <= q(3)*collocate(e3)};
                
                cbox_E_fix={
                    Efix2 == icollocate(e2)
                    Efix3 == icollocate(e3)};
                
                if tau>0  % does not for varflag =1
                    
                    x01 = { icollocate({...
                        xa == subs(xa,x)
                        na == subs(na,x)
                        xs == subs(xs,x)
                        ns == subs(ns,x)
                        xh == subs(xh,x)
                        nh == subs(nh,x)
                        np == subs(np,x)
                        })
                        collocate({e3==subs(e3,x)
                        e2==subs(e2,x)
                        })
                        };
                end
            else
                tomControls e1 e2 e3
                % initial guess
                x01 = { icollocate({...
                    xa == INITCOND(1)/2
                    na == INITCOND(2)/2
                    xs == INITCOND(3)/2
                    ns == INITCOND(4)/2
                    xh == INITCOND(5)/2
                    nh == INITCOND(6)/2
                    np == INITCOND(7)/2
                    })
                    collocate(e1==.1)
                    collocate(e2==.2)
                    collocate(e3==.2)};
                
                cbox_E_free={ 0.0 <= q(1)*collocate(e1)
                    0.0 <= q(2)*collocate(e2)
                    0.0 <= q(3)*collocate(e3)};
                
                cbox_E_fix={ Efix1 == icollocate(e1)
                    Efix2 == icollocate(e2)
                    Efix3 == icollocate(e3)};
                
                if tau>0  % does not work for varflag =2
                    
                    x01 = { icollocate({...
                        xa == subs(xa,x)
                        na == subs(na,x)
                        xs == subs(xs,x)
                        ns == subs(ns,x)
                        xh == subs(xh,x)
                        nh == subs(nh,x)
                        np == subs(np,x)
                        })
                        collocate({e3==subs(e3,x)
                        e2==subs(e2,x)
                        e1==subs(e1,x)})
                        };
                end
                
            end
            
            
            
            
            
            % Call to the parameter file that has the ecological parameters for the
            % different ecological models.
            
            
            CC_Equations_1 % calls a file with all equations
            
            
            if PELICAN==1
                cbox_pelican={0 <= icollocate(np)};
            else
                cbox_pelican={ConsP <= icollocate(np)};
            end
            
            
           
            % Calling the solver using the SNOPT algorithm
            if tau==0 % only used for initial guess
                
                % These were used to solve for the initial guess when the
                % lag was turned off. don't need them now that you are
                % using a saved intitial guess. if you deviate from the
                % base parameters you might need to generate a good initial
                % guess for your case.
                
                %[Solution, Result]=ezsolve(objective, {cbnd,ceq,cbox_stocks,cbox_E_fix,cbox_case,cbox_pelican,cbox_F_pos},x01,options); % 
                %[Solution, Result]=ezsolve(objective, {cbnd,ceq,cbox_stocks,cbox_E_free,cbox_case,cbox_pelican,cbox_F_pos},Solution,options);
                
            else % main set of analysis with tau=4
                
                [Solution, Result]=ezsolve(objective, {cbnd,ceq,cbox_stocks,cbox_pelican,cbox_F_pos},x01,options);
            end
            
            % Below is code used to address issues with the solver when and
            % if needed during the sensitivity analysis
            
            if tau>0&&abs(Result.Inform)>3
                CASE
                PELICAN
                VARFLAG
                Prob = sym2prob(objective,{cbnd,ceq,cbox_stocks,cbox_pelican,cbox_F_pos},Solution,options);
                Prob.SOL.moremem = 1e8;
                Prob.SOL.optPar(35) = 300000;
                Prob.SOL.optPar(36) = 300000;
                Prob.SOL.optPar(30) = 750000;
                Result = tomRun('snopt',Prob,1);
                Solution = getSolution(Result);
                %[Solution, Result]=ezsolve(objective, {cbnd,ceq,cbox_stocks,cbox_E_free,cbox_case,cbox_pelican,cbox_F_pos},Solution,options);
            end
            
            if tau>0&&abs(Result.Inform)>3
                CASE
                PELICAN
                VARFLAG
                disp('Starting over with the fixed path constraint imposed to get Initial Guess')
                [Solution, Result]=ezsolve(objective, {cbnd,ceq,cbox_stocks,cbox_E_fix,cbox_pelican,cbox_F_pos},x01,options);
                
                Prob = sym2prob(objective,{cbnd,ceq,cbox_stocks,cbox_pelican,cbox_F_pos},Solution,options);
                Prob.SOL.moremem = 1e8;
                Prob.SOL.optPar(35) = 50000;
                Prob.SOL.optPar(36) = 50000;
                Prob.SOL.optPar(30) = 500000;
                Result = tomRun('snopt',Prob,1);
                Solution = getSolution(Result);
                [Solution, Result]=ezsolve(objective, {cbnd,ceq,cbox_stocks,cbox_pelican,cbox_F_pos},Solution,options);
            end
        
            
            if tau>0&& abs(Result.Inform)>3
                CASE
                PELICAN
                VARFLAG
                disp('Going to Knitro solver')
                options.solver = 'knitro';
                [Solution, Result]=ezsolve(objective, {cbnd,ceq,cbox_stocks,cbox_pelican,cbox_F_pos},x01,options);
            end
            
            if abs(Result.Inform)>3
                CASE
                PELICAN
                VARFLAG
                disp('Last try')
                Prob = sym2prob(objective,{cbnd,ceq,cbox_stocks,cbox_pelican,cbox_F_pos},Solution,options);
                Prob.SOL.moremem = 1e8;
                Prob.SOL.optPar(35) = 50000;
                Prob.SOL.optPar(36) = 50000;
                Prob.SOL.optPar(30) = 500000;
                Result = tomRun('snopt',Prob,1);
                Solution = getSolution(Result);
            end
            
            
            
            if abs(Result.Inform)<4
                
                SolutionSavesALL{VARFLAG+1,PELICAN,CASE,i}=Solution;
                
                if PELICAN==1
                    SolutionSaves1{VARFLAG+1,PELICAN,CASE,i}=Solution;
                end
                
                % Storing the results with and without pelican
                Timenop{VARFLAG+1,PELICAN,CASE,i}=subs(icollocate(t),Solution);
                Xanop{VARFLAG+1,PELICAN,CASE,i}=subs(icollocate(xa),Solution);
                Xsnop{VARFLAG+1,PELICAN,CASE,i}=subs(icollocate(xs),Solution);
                Xhnop{VARFLAG+1,PELICAN,CASE,i}=subs(icollocate(xh),Solution);
                Nanop{VARFLAG+1,PELICAN,CASE,i}=subs(icollocate(na),Solution);
                Nsnop{VARFLAG+1,PELICAN,CASE,i}=subs(icollocate(ns),Solution);
                Nhnop{VARFLAG+1,PELICAN,CASE,i}=subs(icollocate(nh),Solution);
                Npnop{VARFLAG+1,PELICAN,CASE,i}=subs(icollocate(np),Solution);
                NPV{VARFLAG+1,PELICAN,CASE,i}=-Result.f_k;
                Xa=subs(icollocate(xa),Solution);
                Xs=subs(icollocate(xs),Solution);
                Xh=subs(icollocate(xh),Solution);
                Np=subs(icollocate(np),Solution);
                NPP=subs(icollocate(np),Solution);
                
                if PELICAN==1
                    Opt_pelicanT{Varflag,CASE}=NPP(end);
                end
                
                
                if CASE==1
                    
                    E1nop{VARFLAG+1,PELICAN,CASE,i}=subs(icollocate(e1),Solution);
                    E2nop{VARFLAG+1,PELICAN,CASE,i}=max(0,((P(2)*q(2)*Xs-c0(2))./c1(2)));
                    E3nop{VARFLAG+1,PELICAN,CASE,i}=max(0,((P(3)*q(3)*Xh-c0(3))./c1(3)));
                    
                elseif CASE==2
                    
                    E1nop{VARFLAG+1,PELICAN,CASE,i}=max(0,((P(1)*q(1)*Xa-c0(1))./c1(1)));
                    E2nop{VARFLAG+1,PELICAN,CASE,i}=subs(icollocate(e2),Solution);
                    E3nop{VARFLAG+1,PELICAN,CASE,i}=max(0,((P(3)*q(3)*Xh-c0(3))./c1(3)));
                    
                elseif CASE==3
                    
                    E1nop{VARFLAG+1,PELICAN,CASE,i}=max(0,((P(1)*q(1)*Xa-c0(1))./c1(1)));
                    E2nop{VARFLAG+1,PELICAN,CASE,i}=max(0,((P(2)*q(2)*Xs-c0(2))./c1(2)));
                    E3nop{VARFLAG+1,PELICAN,CASE,i}=subs(icollocate(e3),Solution);
                    
                elseif CASE==4
                    E1nop{VARFLAG+1,PELICAN,CASE,i}=subs(icollocate(e1),Solution);
                    E2nop{VARFLAG+1,PELICAN,CASE,i}=subs(icollocate(e2),Solution);
                    E3nop{VARFLAG+1,PELICAN,CASE,i}=max(0,((P(3)*q(3)*Xh-c0(3))./c1(3)));
                    
                elseif CASE==5
                    E1nop{VARFLAG+1,PELICAN,CASE,i}=subs(icollocate(e1),Solution);
                    E2nop{VARFLAG+1,PELICAN,CASE,i}=max(0,((P(2)*q(2)*Xs-c0(2))./c1(2)));
                    E3nop{VARFLAG+1,PELICAN,CASE,i}=subs(icollocate(e3),Solution);
                    
                elseif CASE==6
                    E1nop{VARFLAG+1,PELICAN,CASE,i}=max(0,((P(1)*q(1)*Xa-c0(1))./c1(1)));
                    E2nop{VARFLAG+1,PELICAN,CASE,i}=subs(icollocate(e2),Solution);
                    E3nop{VARFLAG+1,PELICAN,CASE,i}=subs(icollocate(e3),Solution);
                    
                else
                    E1nop{VARFLAG+1,PELICAN,CASE,i}=subs(icollocate(e1),Solution);
                    E2nop{VARFLAG+1,PELICAN,CASE,i}=subs(icollocate(e2),Solution);
                    E3nop{VARFLAG+1,PELICAN,CASE,i}=subs(icollocate(e3),Solution);
                    
                end
                
            else
                disp('Did not solve!!!!!!!!!!!!')
                CASE
                PELICAN
                VARFLAG
                dbstop  % this stops the program from running any further and is only called when the algorithm did not converge to optimal solution.
                SolutionSaves1{VARFLAG+1,PELICAN,CASE,i}=[];
                
            end
            
            
            if CASE==7
                % only calculating the profit with the case tau=4
                Time=subs(icollocate(t),Solution);
                E1=subs(icollocate(e1),Solution);
                E2=subs(icollocate(e2),Solution);
                E3=subs(icollocate(e3),Solution);
                Xa=subs(icollocate(xa),Solution);
                Xs=subs(icollocate(xs),Solution);
                Xh=subs(icollocate(xh),Solution);
                Np=subs(icollocate(np),Solution);
            end
            
            
        end
    end
    
end


%% CASE 3 only

for PELICAN=1:5  % There is no reason to run all of these cases as Halibut alone can not increase P(T)
      
    for Varflag=1:4
        
        VARFLAG=Varflag-1;
        
        for CASE=3  % Just case 3 from the analysis Halibut only
            
            toms t
            
            % Lag is on
            tau=4;
            options=struct;
            
            if VARFLAG<3
                Phas=tomPhase('Phas',t,t0,T,NN(8));
                options.solver='snopt'; %,[],'Cheb'); %,'fem1s')% 'fem1') %,[],'Cheb');
            else
                Phas=tomPhase('Phas',t,t0,T,NN(8),[],'Cheb');
                options.solver='snopt';
            end
            
            
            
            setPhase(Phas)
            
            if PELICAN>1
                Efix1=.01; Efix2=.01; Efix3=.01;
            else
                Efix1=.15; Efix2=.15; Efix3=.15;
            end
            
            
            
            tomStates xs ns xa na xh nh np
            
            % Halibut only required its own initial guess for each variant
            if PELICAN==1
                if VARFLAG==0
                   
                   load (strcat(PATH,'InitialGuessVarflag_0'))
                   x=SolutionVarFlag{CASE};
                    
                elseif VARFLAG==1
                    
                    load (strcat(PATH,'InitialGuessVarflag_1')) 
                    x = SolutionSaves{VARFLAG+1,PELICAN,CASE,i};
                    
                elseif VARFLAG==2
                    
                    load (strcat(PATH,'InitialGuessVarflag_2'))
                    x = SolutionSaves{VARFLAG+1,PELICAN,CASE,i};
                    
                else
                    
                    load (strcat(PATH,'InitialGuessVarflag_3'))
                    x = SolutionSaves{VARFLAG+1,PELICAN,CASE,i};
                    
                end
            end
            
            if PELICAN>1
                x = SolutionSavesALL{VARFLAG+1,PELICAN-1,CASE,i}; % uses past runs as initial guess
            end
            
            [xazero, xszero, a,b,Mpj,tjuv, theta, vbk,kappa,d,wr,Y, alphaa, alphas, Cmax, Mp, amp_a, amp_s, per_a, per_s, start_a, start_s,...
                recbara, wra, Ma, winfa, kappaa, recbars, wrs, Ms, winfs, kappas,...
                z, phi, Kbp, theta_1, theta_2, theta_3, Gamma_a, Gamma_s, Gamma_other, Mbp, ...
                theta_1s, theta_2s, theta_3s, P,q,c0,c1]=CC_Param_File(flag, VARFLAG, capT,Param,7);
            
            
            tomControls e3
            
            
            cbox_E_free={
                0.0 <= q(3)*collocate(e3)};
            
            cbox_E_fix={Efix3 == icollocate(e3)};
            
         
            x01 = { icollocate({...
                xa == subs(xa,x)
                na == subs(na,x)
                xs == subs(xs,x)
                ns == subs(ns,x)
                xh == subs(xh,x)
                nh == subs(nh,x)
                np == subs(np,x)
                })
                collocate({e3==subs(e3,x)
                })
                };
            
            % Call to the parameter file that has the ecological parameters for the
            % different ecological models.          
            CC_Equations_1 % calls a file with all equations
            
            
            if PELICAN==1||CASE==3
                cbox_pelican=[];
            else
                cbox_pelican={ConsP <= icollocate(np)};
            end
           
            [Solution, Result]=ezsolve(objective, {cbnd,ceq,cbox_stocks,cbox_E_free,cbox_case,cbox_pelican,cbox_F_pos},x01,options);
            
            if tau>0&&abs(Result.Inform)>3
                CASE
                PELICAN
                VARFLAG
                Prob = sym2prob(objective,{cbnd,ceq,cbox_stocks,cbox_E_free,cbox_case,cbox_pelican,cbox_F_pos},Solution,options);
                Prob.SOL.moremem = 1e8;
                Prob.SOL.optPar(35) = 300000;
                Prob.SOL.optPar(36) = 300000;
                Prob.SOL.optPar(30) = 750000;
                Result = tomRun('snopt',Prob,1);
                Solution = getSolution(Result);
                %[Solution, Result]=ezsolve(objective, {cbnd,ceq,cbox_stocks,cbox_E_free,cbox_case,cbox_pelican,cbox_F_pos},Solution,options);
            end
            
            if tau>0&&abs(Result.Inform)>3
                CASE
                PELICAN
                VARFLAG
                disp('Starting over with the fixed path constraint imposed to get Initial Guess')
                [Solution, Result]=ezsolve(objective, {cbnd,ceq,cbox_stocks,cbox_E_fix,cbox_case,cbox_pelican,cbox_F_pos},x01,options);
                
                Prob = sym2prob(objective,{cbnd,ceq,cbox_stocks,cbox_E_free,cbox_case,cbox_pelican,cbox_F_pos},Solution,options);
                Prob.SOL.moremem = 1e8;
                Prob.SOL.optPar(35) = 50000;
                Prob.SOL.optPar(36) = 50000;
                Prob.SOL.optPar(30) = 500000;
                Result = tomRun('snopt',Prob,1);
                Solution = getSolution(Result);
                [Solution, Result]=ezsolve(objective, {cbnd,ceq,cbox_stocks,cbox_E_free,cbox_case,cbox_pelican,cbox_F_pos},Solution,options);
            end
            
            if tau>0&& abs(Result.Inform)>3
                CASE
                PELICAN
                VARFLAG
                disp('Going to Knitro solver')
                options.solver = 'knitro';
                [Solution, Result]=ezsolve(objective, {cbnd,ceq,cbox_stocks,cbox_E_free,cbox_case,cbox_pelican,cbox_F_pos},x01,options);
            end
            
            if abs(Result.Inform)>3
                CASE
                PELICAN
                VARFLAG
                disp('Last try')
                Prob = sym2prob(objective,{cbnd,ceq,cbox_stocks,cbox_E_free,cbox_case,cbox_pelican,cbox_F_pos},Solution,options);
                Prob.SOL.moremem = 1e8;
                Prob.SOL.optPar(35) = 50000;
                Prob.SOL.optPar(36) = 50000;
                Prob.SOL.optPar(30) = 500000;
                Result = tomRun('snopt',Prob,1);
                Solution = getSolution(Result);
            end
            
            if abs(Result.Inform)<4
                
                SolutionSavesALL{VARFLAG+1,PELICAN,CASE,i}=Solution;
                
                if PELICAN==1
                    SolutionSaves1{VARFLAG+1,PELICAN,CASE,i}=Solution;
                end
                
                % Storing the results with and without pelican
                Timenop{VARFLAG+1,PELICAN,CASE,i}=subs(icollocate(t),Solution);
                Xanop{VARFLAG+1,PELICAN,CASE,i}=subs(icollocate(xa),Solution);
                Xsnop{VARFLAG+1,PELICAN,CASE,i}=subs(icollocate(xs),Solution);
                Xhnop{VARFLAG+1,PELICAN,CASE,i}=subs(icollocate(xh),Solution);
                Nanop{VARFLAG+1,PELICAN,CASE,i}=subs(icollocate(na),Solution);
                Nsnop{VARFLAG+1,PELICAN,CASE,i}=subs(icollocate(ns),Solution);
                Nhnop{VARFLAG+1,PELICAN,CASE,i}=subs(icollocate(nh),Solution);
                Npnop{VARFLAG+1,PELICAN,CASE,i}=subs(icollocate(np),Solution);
                NPV{VARFLAG+1,PELICAN,CASE,i}=-Result.f_k;
                Xa=subs(icollocate(xa),Solution);
                Xs=subs(icollocate(xs),Solution);
                Xh=subs(icollocate(xh),Solution);
                Np=subs(icollocate(np),Solution);
                
                NPP=subs(icollocate(np),Solution);
                
                if PELICAN==1
                    Opt_pelicanT{Varflag,CASE}=NPP(end);
                end
                
                
                
                E1nop{VARFLAG+1,PELICAN,CASE,i}=max(0,((P(1)*q(1)*Xa-c0(1))./c1(1)));
                E2nop{VARFLAG+1,PELICAN,CASE,i}=max(0,((P(2)*q(2)*Xs-c0(2))./c1(2)));
                E3nop{VARFLAG+1,PELICAN,CASE,i}=subs(icollocate(e3),Solution);
                
                
                
            else
                disp('Did not solve!!!!!!!!!!!!')
                CASE
                PELICAN
                VARFLAG
                dbstop
                SolutionSaves1{VARFLAG+1,PELICAN,CASE,i}=[];
                
            end
            
            
            
            
        end
    end
    if PELICAN==1
        if VARFLAG==0
            save (strcat(PATH,'InitialGuessVarflag_0_tau4_NPV'), 'NPV')
            save (strcat(PATH,'InitialGuessVarflag_0_tau4'), 'SolutionSaves1')
        elseif VARFLAG==1
            save (strcat(PATH,'InitialGuessVarflag_1_tau4_NPV'), 'NPV')
            save (strcat(PATH,'InitialGuessVarflag_1_tau4'), 'SolutionSaves1')
        elseif VARFLAG==2
            save (strcat(PATH,'InitialGuessVarflag_2_tau4_NPV'), 'NPV')
            save (strcat(PATH,'InitialGuessVarflag_2_tau4'), 'SolutionSaves1')
        else
            save (strcat(PATH,'InitialGuessVarflag_3_tau4_NPV'), 'NPV')
            save (strcat(PATH,'InitialGuessVarflag_3_tau4'), 'SolutionSaves1')
        end
        

        save (strcat(PATH,'optimal_pelicanT'), 'Opt_pelicanT')
           
    end
end


DATA=[Timenop Xanop Nanop Xsnop Nsnop Xhnop Nhnop Npnop E1nop E2nop E3nop];
save (strcat(PATH,'SolutionSavesAll'), 'SolutionSavesALL')
save (strcat(PATH,'Solution_Tim'),'DATA')
save (strcat(PATH,'Solution_NPV'),'NPV')


%% NPV calculations
time=Time(Time<30); % Using the first 30 periods and not full 40 to avoid 
% the influence of the transversality condition on the NPV calculations.
% See dynamics to note that systems are settled down

L=length(time); trun=1:L;
for VarFlag=1:4
    VARFLAG=VarFlag-1;
   
    for CASE=1:7
     
     E1= E1nop{VARFLAG+1,PELICAN,CASE,i}(trun);
     E2=E2nop{VARFLAG+1,PELICAN,CASE,i}(trun);
     E3=E3nop{VARFLAG+1,PELICAN,CASE,i}(trun);
     Xa=Xanop{VARFLAG+1,PELICAN,CASE,i}(trun);
     Xs=Xsnop{VARFLAG+1,PELICAN,CASE,i}(trun);
     Xh=Xhnop{VARFLAG+1,PELICAN,CASE,i}(trun);   
        Ra=P(1)*q(1)*E1.*(Xa.^BETA(1))-c0(1)*(E1)-c1(1)*(E1.^2); 
        Rs=P(2)*q(2)*E2.*(Xs.^BETA(2))-c0(2)*(E2)-c1(2)*(E2.^2);
        Rh=P(3)*q(3)*E3.*(Xh.^BETA(3))-c0(3)*(E3)-c1(3)*(E3.^2);
        
        f = fittype('pchipinterp');
                     g1a=fit(time,exp(-dis.*time).*Ra,f);
                     int1a=integrate(g1a,0:1:time(end),0);
                     g2s=fit(time,exp(-dis.*time).*Rs,f);
                     int2s=integrate(g2s,0:1:time(end),0);
                     g3h=fit(time,exp(-dis.*time).*Rh,f);
                     int3h=integrate(g3h,0:1:time(end),0);
                     
                     NpvA(VarFlag,CASE)=int1a(end);
                     NpvS(VarFlag,CASE)=int2s(end);
                     NpvH(VarFlag,CASE)=int3h(end)*10; % scaled up for graphing purposes, does not change the relative results
        
    end
end

NPVDATA=[NpvA,NpvS,NpvH];
save (strcat(PATH,'NetPresent'),'NPVDATA')

%% This can be run separate to save time as all of the data is saved above and then loaded back in the file
CA_Current_Analysis_2
