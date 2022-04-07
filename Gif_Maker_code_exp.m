%COVID Gif Maker
clear;close all; clc
%from COVIDexp_MainSwitch.m

    %INPUT parameters:
        %t is  time
        %beta is the transmission rate
        %sigma is the rate at which exposed individuals become infected
        %omega is the rate at which asymptomatic individuals become
        %symptomatic or recovered 
        %epsilon is the rate at which symptomatic people become quarantined
        %gamma is the rate at which quarantined individuals become
        %hospitalized or recovered 
        %1/aR is the duration at which hospitalized individuals who
        %recover remain hospitalized
        %1/aD is the duration at which hospitalized individuals who die
        %remain hospitalized
        %x is the rate at which hospitalized individuals die which is a
        %function of 
        %X(x hat) is the lethality rate of hospitalized individuals with
        %access to an ICU bed 
        %p is the fraction of symptomatic individudals out of exposed
        %individuals 
        %q is the fraction of hospitalized individuals out of quarantined
        %individuals 
        %v is the fraction of recoved individuals out of asymptomatic
        %indivuals 
        %a is the fraction of susceptible people that enter confinement
        %b is the fraction of confined individuals that become susceptible
        %tlock is the lockdown time in days
        %tlift is the lockdown lift time in days 
        
        
    %OUTPUT:
        %N is the population array
        %N(1) is the population of susceptibles = S
        %N(2) is the population in the incubation period of disease progression = E
        %N(3) is the population of infectious individuals who do not show symptoms of COVID-19 = Ia 
        %N(4) is the population of infection individuals who show symptoms = Is
        %N(5) is the population of infectious individuals, hospitalized with symptoms = H
        %N(6) is the population of symptomatic infectious individuals who are isolated = Q
        %N(7) is the population of individuals who survived COVID-19 = R
        %N(8) is the population of individuals who did not survive = D
        %N(9) is the population of individuals who have not been infected with
        %COVID-19 and have isolated 
        %Ntot=total population = N(1)+N(2)+N(3)+N(4)+N(5)+N(6)+N(7)-N(8)+N(9)
        %Phi is the rate at which susceptible individuals become confined 
        %Psi is the rate at which confined individuals become susceptible
    
    %Parameter values
    lambda = [0:0.2:1];
    beta = (1.5)./(2+exp(15.*(lambda-0.5)));  
    sigma = 1/6;  
    omega = 1/14;
    epsilon = .5;   
    gamma = 1/5;  
    aR= 1/12;
    aD= 1/14;
    x= 0.04;
    p= 0.6; 
    q=.9;
    v=.8; 
    a=.5;
    b=1;
    tlock=20;
    tlift=80;
    

    N10=9000099;      %Initial population of S 
    %10,000,099 = average state population
    N20=100;      %Initial population of E
    %100 = paper inital exposed value
    %other populations N0=0 for paper values 
    N30=0;      %Initial population of Ia
    N40=0;      %Initial population of Is
    N50=0;      %Initial population of H 
    N60=0;      %Initial population of Q
    N70=0;      %Initial population of R
    N80=0;      %Initial population of D
    N90=0;      %Initial population of C
    tend=3*365;    %Simulation length in days
   

%creating the gif (from matlab help)
beta_var=[beta(1),beta(2),beta(3),beta(4),beta(5),beta(6)] %range for parameter sweep
nImages = length(beta_var);
fig = figure;
for idx = 1:nImages
    %runs simulation from 0 to tlock    
    N0=[N10 N20 N30 N40 N50 N60 N70 N80 N90];
    [t1,N1] = ode45(@(t,N) COVIDexp_ODE(t,N,beta_var(idx),sigma,omega,epsilon,gamma,aR,aD,x,p,q,v,a,b,tlock,tlift),[0:1:tlock],N0);
    Ntot1=N1(:,1)+N1(:,2)+N1(:,3)+N1(:,4)+N1(:,5)+N1(:,6)+N1(:,7)-N1(:,8)+N1(:,9);

    %runs simulation from tlock to tlift 
    %initial condtions are end conditions of last run & confined population
    %increased according to a 
    N0=[N1(end,1)*(1-a) N1(end,2) N1(end,3) N1(end,4) N1(end,5) N1(end,6) N1(end,7) N1(end,8) (N1(end,9)+ a*N1(end,1))];
    [t2,N2] = ode45(@(t,N) COVIDexp_ODE(t,N,beta_var(idx),sigma,omega,epsilon,gamma,aR,aD,x,p,q,v,a,b,tlock,tlift),[tlock:1:tlift],N0);
    Ntot2=N2(:,1)+N2(:,2)+N2(:,3)+N2(:,4)+N2(:,5)+N2(:,6)+N2(:,7)-N2(:,8)+N2(:,9);

    %runs simulation from tlift to tend
    %inital conditions are end conditions of last run & confined population
    %decreased according to b 
    N0=[(N2(end,1)+ b*N2(end,9)) N2(end,2) N2(end,3) N2(end,4) N2(end,5) N2(end,6) N2(end,7) N2(end,8) N2(end,9)*(1-b)];
    [t3,N3] = ode45(@(t,N) COVIDexp_ODE(t,N,beta_var(idx),sigma,omega,epsilon,gamma,aR,aD,x,p,q,v,a,b,tlock,tlift),[tlift:1:tend],N0);
    Ntot3=N3(:,1)+N3(:,2)+N3(:,3)+N3(:,4)+N3(:,5)+N3(:,6)+N3(:,7)-N3(:,8)+N3(:,9);

    %puts all three runs together 
    Ntot= [Ntot1; Ntot2; Ntot3];
    t= [t1;t2;t3];
    N= [ N1; N2; N3];
    
    plot(t,N(:,2),'linewidth',2); %this plots exposed population 
    xlabel('Time (days)');
    ylabel('Exposed Population');
    axis([0 500 0 3.5*10^6]); % need to be adjusted with every parameter, but must be set for gif to keep constant axes 
    title('Beta Increasing Parameter Sweep');
    legend('beta = 0','beta = .2','beta = .4','beta = .6','beta = .8','beta = 1');
    hold on
    drawnow 
    % Capture the plot as an image 
    frame = getframe(fig); 
    im{idx} = frame2im(frame);
end 

close;

filename = 'testAnimated.gif'; %names the gif file, saves in same folder
for idx = 1:nImages
    [A,map] = rgb2ind(im{idx},256); 
    % Write to the GIF File 
    if  idx == 1
        imwrite(A,map,filename,'gif', 'Loopcount',inf,'DelayTime',1); 
    else 
        imwrite(A,map,filename,'gif','WriteMode','append', 'DelayTime', 1); 
    end 
    
end

web(filename);%displays gif in web so you can check that it worked easily 