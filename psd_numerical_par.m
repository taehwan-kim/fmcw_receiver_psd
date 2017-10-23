clc; clear; close all;

set(0,'DefaultFigureWindowStyle','normal');
addpath ./lib;

%% constants

tau_c       = 0.16e-6;      % coherence time 
tau         = 5*tau_c;    % round-trip delay
wp          = 2*pi*60e6;     % tuning bandwidth

Deltaw      = 2/tau_c;      % frequency noise density at DC

%% pll related constants

Deltaw_pll = Deltaw * 0.0001;
wp1         = 2*pi*1e5;
wz1         = 2*pi*1e3;

%% laser phase noise PSD

% original laser with white frequency noise assumption 
s_phin = @(w) Deltaw ./ w.^2;

% colored frequency noise from tuning bandwidth limitation 
tf_lp = @(w) 1 ./ (1+(w/wp).^2) ;
s_phin_lp = @(w) s_phin(w) .* tf_lp(w);

% noise-shaped case for pll
tf_pll = @(w) 0.0001 * (1+(w/wz1).^2) ./ (1+(w/wp).^2) ./ (1+(w/wp1).^2);
s_phin_pll = @(w) s_phin(w) .* tf_pll(w);

% photocurrent phase noise density from the laser phase noise PSD
s_dphin = @(w, u) 4*sin(w .* u/2).^2 .* s_phin(w);
s_dphin_lp = @(w, u) 4*sin(w .* u/2).^2 .* s_phin_lp(w);
s_dphin_pll = @(w, u) 4*sin(w .* u/2).^2 .* s_phin_pll(w);

% photocurrent noise 
s_theta = @(w, u, tau) 2*s_dphin(w, tau) + 2*s_dphin(w, u) - s_dphin(w, u+tau) - s_dphin(w, u-tau);
s_theta_lp = @(w, u, tau) 2*s_dphin_lp(w, tau) + 2*s_dphin_lp(w, u) - s_dphin_lp(w, u+tau) - s_dphin_lp(w, u-tau);
s_theta_pll = @(w, u, tau) 2*s_dphin_pll(w, tau) + 2*s_dphin_pll(w, u) - s_dphin_pll(w, u+tau) - s_dphin_pll(w, u-tau);

% tmax = 1e-3;
% wmin = 1e3;
wmax = 1000e6;
samples = 3000;
sampling_period = 2e-9;
tmax = samples * sampling_period;
frequency_bin = 1/(2*tmax);

%s = (-1*tmax:sampling_period:tmax);
s = (0:sampling_period:tmax);
f = (-1*samples*frequency_bin+frequency_bin/2:frequency_bin:samples*frequency_bin-frequency_bin/2);
% linspace(-1*tmax,tmax,samples+1);
ws=10.^(1:1:20);

sigma_dphin = zeros(1,length(s));
sigma_theta = zeros(1,length(s));
sigma_dphin_lp = zeros(1,length(s));
sigma_theta_lp = zeros(1,length(s));
sigma_dphin_pll = zeros(1,length(s));
sigma_theta_pll = zeros(1,length(s));

parfor i=1:length(s)

    sigma_dphin(i) = (1/(pi)) * integral(@(w)s_dphin(w, s(i)), 0, Inf);
    sigma_theta(i) = (1/(pi)) * integral(@(w)s_theta(w, s(i), tau), 0, Inf);
    
    sigma_dphin_lp(i) = (1/(pi)) * integral(@(w)s_dphin_lp(w, s(i)), 0, Inf);
    sigma_theta_lp(i) = (1/(pi)) * integral(@(w)s_theta_lp(w, s(i), tau), 0, Inf);
    
    sigma_dphin_pll(i) = (1/(pi)) * integral(@(w)s_dphin_pll(w, s(i)), 0, wmax,'ArrayValued',true);
    sigma_theta_pll(i) = (1/(pi)) * integral(@(w)s_theta_pll(w, s(i), tau), 0, wmax,'ArrayValued',true);

end

sigma_dphin = [fliplr(sigma_dphin(2:end)), sigma_dphin];
sigma_theta = [fliplr(sigma_theta(2:end)), sigma_theta];

sigma_dphin_lp = [fliplr(sigma_dphin_lp(2:end)), sigma_dphin_lp];
sigma_theta_lp = [fliplr(sigma_theta_lp(2:end)), sigma_theta_lp];

sigma_dphin_pll = [fliplr(sigma_dphin_pll(2:end)), sigma_dphin_pll];
sigma_theta_pll = [fliplr(sigma_theta_pll(2:end)), sigma_theta_pll];


r_i0 = exp(-0.5*sigma_theta);
s_i0 = 10*log10(abs(fft(r_i0(1:end-1))));
s_i0 = [s_i0(samples+2:end), s_i0(1:samples+1)];

r_i0_lp = exp(-0.5*sigma_theta_lp);
s_i0_lp = 10*log10(abs(fft(r_i0_lp(1:end-1))));
s_i0_lp = [s_i0_lp(samples+2:end), s_i0_lp(1:samples+1)];

r_i0_pll = exp(-0.5*sigma_theta_pll);
s_i0_pll = 10*log10(abs(fft(r_i0_pll(1:end-1))));
s_i0_pll = [s_i0_pll(samples+2:end), s_i0_pll(1:samples+1)];

% h = figure();
% plot(f/1e6, s_i0);
% hold on;
% plot(f/1e6, s_i0_lp);
% hold on;
% plot(f/1e6, s_i0_pll);
% hold on;
% plot(f(3000)/1e6, s_i0(3000));
% hold on;
% plot(f(3000)/1e6, s_i0_lp(3000));
% hold on;
% plot(f(3000)/1e6, s_i0_pll(3000));
% plot_graph(h, 3, 'temp.pdf');
% 
% t1 = text(f(3000)/1e6,s_i0(3000),' ? temp');
% set(t1,'FontSize',15);
% set(t1,'FontName','Helvetica Neue');
% % set(t1,'Interpreter','latex');
% % set(t1,'String', '? temp');
% 
% t2 = text(f(3000)/1e6,s_i0_lp(3000),'temp ? ');
% set(t2,'FontSize',15);
% set(t2,'HorizontalAlignment','right');
% set(t2,'FontName','Helvetica Neue');

