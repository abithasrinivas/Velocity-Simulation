%Velocity Simulation to observe velocity profile 
%
% Input functions required : blochsim 
%                            newvelsim
%                            gensech180
%
%Author: Sai Abitha Srinivas 
clear all;
close all;

accel = 0 ; % in cm/ms^2
vel = 0  ; % in cm/ms
pos0 = 0 ;  % cm

homogeneity = 'perfect'; % 'perfect'  % 'bad_B0_B1'

switch homogeneity
    case 'bad_B0'
        off_resonance =    100e-3 / 42576 ; % (KHz --> Tesla) ... gammabar is in kHz/T
        
    case 'bad_B0_B1'
        off_resonance =   100e-3 / 42576 ; % (KHz --> Tesla) ... gammabar is in kHz/T
        
    case 'perfect'
        off_resonance = 0%  100e-3 / 42576 ; % (KHz --> Tesla) ... gammabar is in kHz/T
end

segmentpulse = 'bhard20';%change as required
refpulse = 'bhard180' ; %change as required
n = 8; %change as per segmentpulse FA
dt = 1e-3;

[B1 Gz] =newvelsim(segmentpulse,refpulse,n,dt)
Mzfinal =[]; 
isControl = 0; %change as reqd 

% control case: ASL has two cases 
if isControl
    Gz = abs(Gz);
end
 
Bx = real(B1);
By = imag(B1);

NSTEPS = length(B1);

% total duration of the simulation interval (in ms)
duration = NSTEPS*dt;  % ms.
vel_range =[-100:100]*1e-3; 
 
for vel = vel_range  % cm / msec
    
    t = linspace(0,duration, NSTEPS)'; % mseconds.
    zpos = pos0 + vel*t + 0.5*accel*(t.^2);
    Bz = zpos.*Gz;
    Bz = Bz + off_resonance;  
    T1 = 2290;  %ms
    T2 = 68;   %ms        
    beff = [Bx By  Bz];
    Mi = [0 0 1]';
    M = blochsim(Mi, beff, T1, T2, dt, NSTEPS);     
    Mzfinal=[Mzfinal; M(end,3)];

end

figure(3)
t=[0:length(Gz)-1]*dt;
subplot(321)
area(t, Gz/max(Gz));
grid on
subplot(323)
plot(t, abs(B1)/max(abs(B1)),'r');
grid on
subplot(325) 
plot(t, angle(B1),'g');
grid on
hold on 

subplot(122)
hold on
plot(vel_range' * 1e3, Mzfinal) 
hold on 
plot ( 0, min (Mzfinal), '*');
Profiledepth = min(Mzfinal);
[d0 dist0] = min(abs(Mzfinal(80:100)-0));
[d1 dist1] = min(abs(Mzfinal(100:120)-0));
plot ([dist0 dist1], 'g','LineWidth', 2);
axis ([ -100 100 -1 1])
grid on

