%initialize time
step = .1;%time step value
t=0:step:100; %time(t) from 0 to 100 miliseconds with a time step of .01

%given information
I(1:499) = .10;
I(500:numel(t)) = .15;

gbarNa=120; %socium channel conductance (mSeimens/cm^2)
gbarK= 36; %potasium conductance (mSeimens/cm^2)
gbarL= 0.3; %leakage conductance (mSeimens/cm^2)
Ek= -12; %potasium potential (mV)
ENa= 115; %sodium potential (mV)
EL= 10.6; %leakage potential (mV)
Vm= 0; %resting membrane potential(mV)
Cm=1; %membrane capacitance (microfarad/cm^2)

%Gating variable equations%
alpha_m = 0.1*((25-Vm)/exp((25-Vm)/10)-1);
beta_m = 4*exp(-Vm/18);

alpha_n = 0.1*((10-Vm)/exp((10-Vm)/10)-1);
beta_n = .125*exp(-Vm/80);

alpha_h = .07*exp(-Vm/20);
beta_h = 1/(exp((30-Vm)/10)+1);
%Initializing values%
m(1)= alpha_m/(alpha_m + beta_m);
n(1)= alpha_n/(alpha_n + beta_n);
h(1)= alpha_h/(alpha_h + beta_h);

%Derivative equations need solving. Derivative equations will need to be
%evaludated at each time step from 0 to 100miliseconds. 

for time=1:numel(t)-1
    %first diff EQ to solve: dm/dt dn/dt & dh/dt for each time point
    alpha_m(time) = 0.1*((25-Vm(time))/exp((25-Vm(time)/10)-1));
    beta_m(time) = 4*exp(-Vm(time)/18);

    alpha_n(time) = 0.1*((10-Vm(time))/exp((10-Vm(time)/10)-1));
    beta_n(time) = .125*exp(-Vm(time)/80);

    alpha_h(time) = .07*exp(-Vm(time)/20);
    beta_h(time) = 1/(exp((30-Vm(time)/10)+1));



%Current Equations%
I_K = (n(time)^4) * gbarK * (Vm(time)-Ek);
I_Na = (m(time)^3) * gbarNa * h(time) * (Vm(time)-ENa);
I_L = gbarL * (Vm(time)-EL);
I_ion = I(time)- I_K - I_Na - I_L;

%Euler Approximation of Derivatives%
Vm(time+1) = Vm(time) + step * I_ion/Cm; %membrane voltage at a specified time
m(time+1) = m(time) + step*(alpha_m(time) * (1-m(time)) - beta_m(time) * m(time));
n(time+1) = n(time) + step*(alpha_n(time) * (1-n(time)) - beta_n(time) * n(time));
h(time+1) = h(time) + step*(alpha_h(time) * (1-h(time)) - beta_h(time) * h(time));
end
Vm = Vm - 70;


%membrane voltage (Vm) plot

plot(t,Vm)
hold on
legend({'Membrane Voltage'})
xlabel('time (miliseconds)')
ylabel('Voltage (mili-volts)')
title('Membrane Volgate over time')

%plotting conductance (g)

figure

hold on
p2 = plot(t,gbarNa*(m.^3).*h, 'r')
xlabel('time (miliseconds)')
ylabel('Conductance')
title('conductance of K+ and Na+ over time')
