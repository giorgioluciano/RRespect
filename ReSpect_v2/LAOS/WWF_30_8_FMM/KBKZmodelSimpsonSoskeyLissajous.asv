function [time2, strain, strainRt , tau12t2, ElStr, VisStr] = KBKZmodelSimpsonSoskeyLissajous(para)

% KBKZ model from Song 2020 

% gi and tau i from reSpect code of Shanbhag

gi = [177464	54980	35272.4	22493.7	16039.9	13572.9	9004.74 ...
    6926.17	527.529	2.7341  ];

taui = [0.00020788	0.00288704	0.0234895	0.163525	1.03726 ...
    5.24526	15.6349	35.5925	167.211	4810.48];


warning('off');
PointsPerPeriod = 125;
NoPeriods = 400;
CutPeriods = 300;

etai = taui.*gi;

omega = para(1);
Amp = para(2);
aSos = para(3);
bSos = para(4);
mint = 0;
maxt = 2*pi/omega*NoPeriods;
tspan = [mint maxt];
tstep=2*pi/omega/PointsPerPeriod;

time = mint:tstep:maxt;

dim = size(time, 2);

tau12 = zeros(1, dim);

PreSum = @(tp, t) Amp.*(sin(omega*t )-sin(omega*tp)) .*...
        1.0 ./(1.0+aSos .*(abs( Amp.*(sin(omega*t )- ...
        sin(omega*tp)) )).^bSos );
    
Integrand = @(tp, t, PreSumn) PreSumn.*sum(gi./taui.*exp(-(t-tp)./taui) );

SimpsonOdd = @(integrand, h) (integrand(1) + ... 
    4*sum(integrand(2:2:end-1) ) + 2*sum(integrand(3:2:end-2) ) ...
    + integrand(end))*h/3;
SimpsonThreeEight = @(integrand, h) 3*h/8*(integrand(1) + ...
    3*integrand(2)+3*integrand(3)+integrand(4)  );
for i = 1:dim
   integrand2 = zeros(1,i);
   time2 = mint:tstep:time(i);
   PreSumn = PreSum(time2, time(i));
    for j = 1:i
     integrand2(j) = Integrand(time(j), time(i), PreSumn(j) );
    end
    
    
    if i == 1
        tau12(i) = 0;
    elseif i == 2
        tau12(i) = (integrand2(1)+integrand2(2))/2*tstep;
    elseif mod(i,2) == 1 
        tau12(i) = SimpsonOdd(integrand2, tstep);
    else tau12(i) = SimpsonOdd(integrand2(1:end-3), tstep) + ...
            SimpsonThreeEight(integrand2(end-3:end), tstep);
    end
  
  % integral(@(tp) Integrand(tp, time(i) ), 0, time(i));

end

etas = 0.01;
tau12s = etas*Amp*omega*cos(omega*time);
tau12t = tau12 + tau12s;

time2 = time(PointsPerPeriod*CutPeriods+1 :end);
tau12t2 = tau12t(PointsPerPeriod*CutPeriods + 1: end) ;

strain = Amp*sin(omega.*time2 );
strainRt = Amp*omega*cos(omega.*time2 );

Fs = 1000;            % Sampling frequency  
% Fs = size(time,2)/time(end);
Fs = size(time2,2)/(2*pi*(NoPeriods-CutPeriods) ); 
% 100 is no. of periods from the Shear Stress Calc algorithm
T = 1/Fs;             % Sampling period       
L = 1500;             % Length of signal
L = size(time2,2);
t = (0:L-1)*T;        % Time vector
S = 0.7*sin(2*pi*50*t+50) + sin(2*pi*120*t);
S = tau12t2;
X = S ;
% plot(1000*t(1:L),X(1:L))
% title('Signal Corrupted with Zero-Mean Random Noise')
% xlabel('t (milliseconds)')
% ylabel('X(t)')

% saveas(gcf, 'sinewaveexample.pdf');
Y = fft(X);

P2 = abs(Y/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);

aSos = 2*pi*Fs*(0:(L/2))/L;
% plot(f,P1) 
% title('Single-Sided Amplitude Spectrum of X(t)')
% xlabel('f (Hz)');
% ylabel('|P1(f)|');
% axis([0 6 0 100] );
% axis 'auto y';
harmInd = 1:2:61;
harmonics = interp1(aSos,P1, harmInd);

harmonics = harmonics/Amp;
% disp('abs harmonics')
% disp(harmonics);

P2r = real(Y/L);
P1r = P2r(1:L/2+1);
P1r(2:end-1) = 2*P1r(2:end-1);

harmonics = interp1(aSos,P1r, harmInd);
harmonicsR = harmonics/Amp;
% disp('real modulus harmonics')
% disp(harmonicsR/Amp);

P2i = imag(Y/L);
P1i = P2i(1:L/2+1);
P1i(2:end-1) = 2*P1i(2:end-1);

harmonics = interp1(aSos,P1i, harmInd);
harmonicsI = -harmonics/Amp;
% disp('imaginary Modulus harmonics')
% disp(harmonicsI/Amp);

% saveas(gcf, 'fourierexample.pdf');

delta = atan(harmonicsR(1)/harmonicsI(1));
time2 = time(PointsPerPeriod*(NoPeriods-1)+1 :end);
tau12t2 = tau12t(PointsPerPeriod*(NoPeriods) + 1: end) ;

strain = Amp*sin(omega.*time2 );
strainRt = Amp*omega*cos(omega.*time2 );

ElStr = 0.0;
VisStr = 0.0;
dimHarm = size(harmonicsI, 2);
for i = 1:dimHarm
    ElStr = ElStr + Amp*(harmonicsI(i)*sin( (2*(i-1)+1)*omega.*time2));
    VisStr = VisStr + Amp*(harmonicsR(i)*cos( (2*(i-1)+1)*omega.*time2));
    
end

end