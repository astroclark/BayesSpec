% plot_samps.m
%
% James Clark, jclark.astro@gmail.com
%

function plot_samps(Data)

freq=Data.samples.freq;
tau=Data.samples.tau;
amp=Data.samples.amp;
t0=Data.samples.t0;


figure('Position', [100, 50, 1000, 800]);

% --- Frequency
[f,xi]=ksdensity(freq,'support','positive');
subplot(4,3,1)
hist(freq,50);
hold on
xlabel('frequency')
% plot([1500,1500],get(gca,'Ylim'),'r')

subplot(4,3,2)
plot(xi,f)
hold on
xlabel('frequency')
% plot([1500,1500],get(gca,'Ylim'),'r')

subplot(4,3,3)
plot(freq)
hold on
xlabel('sample #')
ylabel('frequency')
% plot(get(gca,'Xlim'),[1500,1500],'r')
% 

% --- Damping Time
[f,xi]=ksdensity(tau,'support','positive');
subplot(4,3,4)
hist(tau,50);
hold on
xlabel('tau')
% plot([0.25,0.25],get(gca,'Ylim'),'r')

subplot(4,3,5)
plot(xi,f)
hold on
xlabel('tau')
% plot([0.25,0.25],get(gca,'Ylim'),'r')

subplot(4,3,6)
plot(tau)
hold on
xlabel('sample #')
ylabel('tau')
% plot(get(gca,'Ylim'),[0.25,0.25],'r')


% --- Start Time
[f,xi]=ksdensity(amp,'support','positive');
subplot(4,3,7)
hist(amp,50);
hold on
xlabel('amp')
% plot([0.25,0.25],get(gca,'Ylim'),'r')

subplot(4,3,8)
plot(xi,f)
hold on
xlabel('amp')
% plot([0.25,0.25],get(gca,'Ylim'),'r')

subplot(4,3,9)
plot(amp)
hold on
xlabel('sample #')
ylabel('amp')
% plot(get(gca,'Ylim'),[0.25,0.25],'r')


% --- Amplitude
[f,xi]=ksdensity(t0,'support','positive');
subplot(4,3,10)
hist(t0,50);
hold on
xlabel('t0')
% plot([3,3],get(gca,'Ylim'),'r')

subplot(4,3,11)
plot(xi,f)
hold on
xlabel('t0')
% plot([3,3],get(gca,'Ylim'),'r')

subplot(4,3,12)
plot(t0)
hold on
xlabel('sample #')
ylabel('t0')
% plot(get(gca,'Ylim'),[3,3],'r')
