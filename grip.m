%% ================================================
%%  SIR MODEL ŠIRENJA GRIPA - Seminarski rad
%%  Osnove matematičkog modeliranja
%% ================================================
 
clear; clc; close all;
 
%% --- 1. PODACI ---
N = 15000;
t_pod = (0:29)';
 
I_pod = [150; 205; 280; 378; 508; 674; 887; 1150; 1454; 1797; ...
         2150; 2474; 2726; 2861; 2855; 2724; 2492; 2197; 1882; 1573; ...
         1288; 1043; 832; 663; 519; 406; 318; 246; 193; 148];
 
R_pod = [0; 51; 120; 214; 341; 512; 740; 1040; 1428; 1920; ...
         2528; 3255; 4091; 5012; 5978; 6943; 7864; 8706; 9449; 10084; ...
         10615; 11051; 11403; 11685; 11908; 12084; 12222; 12330; 12414; 12479];
 
S_pod = N - I_pod - R_pod;
 
%% --- 2. PROCENA PARAMETARA ---
% Procena r iz dR/dt = r * I
r_niz = (R_pod(2:end) - R_pod(1:end-1)) ./ I_pod(1:end-1);
r = mean(r_niz);
fprintf('Procenjeni r = %.4f dan^-1  (prosecno trajanje bolesti = %.2f dana)\n', r, 1/r);
 
% Procena a iz dI/dt = a*S*I - r*I
dI = I_pod(2:end) - I_pod(1:end-1);
a_niz = (dI + r*I_pod(1:end-1)) ./ (S_pod(1:end-1).*I_pod(1:end-1)/N);
a = mean(a_niz(1:13));   % prosek preko faze rasta (dani 0-12)
fprintf('Procenjeni a = %.4f dan^-1\n', a);
 
S0 = S_pod(1)/N;
R0_val = a * S0 / r;
fprintf('Osnovni reprodukcioni broj R0 = %.3f\n\n', R0_val);
 
%% --- 3. SIR model ---
sir_ode = @(t, y, aa, rr) [-aa*y(1)*y(2); ...
                             aa*y(1)*y(2) - rr*y(2); ...
                             rr*y(2)];
 
y0 = [S_pod(1); I_pod(1); 0] / N;   % pocetni uslovi (udjeli)
tspan = [0 29];
opts = odeset('RelTol',1e-8, 'AbsTol',1e-10);
 
%% --- 4. Simulacija svih scenarija ---
 
% Scenarij 0: Osnovni
[t0, Y0] = ode45(@(t,y) sir_ode(t,y,a,r), tspan, y0, opts);
 
% Scenarij 1: Vakcinacija 20%
% 20% populacije prelazi iz S u R pre pocetka sezone (trajni imunitet)
y1    = y0;
y1(1) = y1(1) - 0.20;
y1(3) = y1(3) + 0.20;
R0_init_vak = 0.20;   % inicijalno vakcinisani (nisu bili bolesni)
[t1, Y1] = ode45(@(t,y) sir_ode(t,y,a,r), tspan, y1, opts);
 
% Scenarij 2: Karantin (20% zarazenih izolovano -> a_ef = 0.8*a)
a2 = 0.8 * a;
[t2, Y2] = ode45(@(t,y) sir_ode(t,y,a2,r), tspan, y0, opts);
 
% Scenarij 3: Lekovi (r povecano za 20%)
r3 = 1.20 * r;
[t3, Y3] = ode45(@(t,y) sir_ode(t,y,a,r3), tspan, y0, opts);
 
%% --- 5. Grafik 1: Osnovni scenario + podaci ---
figure('Position',[100 100 820 520]);
plot(t0, Y0(:,1)*N, 'b-', 'LineWidth', 2.2); hold on;
plot(t0, Y0(:,2)*N, 'r-', 'LineWidth', 2.2);
plot(t0, Y0(:,3)*N, 'g-', 'LineWidth', 2.2);
plot(t_pod, S_pod, 'bs', 'MarkerFaceColor', 'b', 'MarkerSize', 5);
plot(t_pod, I_pod, 'rs', 'MarkerFaceColor', 'r', 'MarkerSize', 5);
plot(t_pod, R_pod, 'gs', 'MarkerFaceColor', 'g', 'MarkerSize', 5);
xlabel('Vreme (dani)', 'FontSize', 13);
ylabel('Broj osoba', 'FontSize', 13);
title(sprintf('SIR model - Osnovni scenario (a=%.3f, r=%.3f, R_0=%.2f)', a, r, R0_val), ...
      'FontSize', 14);
legend('S(t) - model', 'I(t) - model', 'R(t) - model', ...
       'S - podaci', 'I - podaci', 'R - podaci', 'Location', 'east');
grid on;
saveas(gcf, 'fig_basic.eps', 'epsc');
 
%% --- 6. Grafik 2: Poredjenje svih scenarija (broj zarazenih) ---
figure('Position',[100 100 860 520]);
plot(t0, Y0(:,2)*N, 'k-',  'LineWidth', 2.5); hold on;
plot(t1, Y1(:,2)*N, 'b--', 'LineWidth', 2.2);
plot(t2, Y2(:,2)*N, 'r:',  'LineWidth', 2.4);
plot(t3, Y3(:,2)*N, 'm-.', 'LineWidth', 2.2);
xlabel('Vreme (dani)', 'FontSize', 13);
ylabel('Broj zarazenih osoba', 'FontSize', 13);
title('Poredjenje scenarija kontrole sirenja gripa', 'FontSize', 14);
legend(sprintf('Osnovni scenario (R_0 = %.2f)', R0_val), ...
       '1. Vakcinacija 20%', ...
       '2. Karantin 20%', ...
       '3. Lekovi (r +20%)', ...
       'Location', 'northeast', 'FontSize', 11);
grid on;
saveas(gcf, 'fig_scenariji.eps', 'epsc');
 
%% --- 7. Tabelarni prikaz rezultata ---
fprintf('\n%-12s %6s %10s %8s %14s\n', 'Scenario', 'R0', 'max I', 'dan max', 'Ukupno zar.');
scen     = {'Osnovni', 'Vakcinacija', 'Karantin', 'Lekovi'};
Yarr     = {Y0, Y1, Y2, Y3};
tarr     = {t0, t1, t2, t3};
R0s      = [a*S0/r,  a*(S0-0.2)/r,  a2*S0/r,  a*S0/r3];
% R0_init: udeo koji je na pocetku u R a nije prebolelo bolest (samo vakcinacija)
R0_inits = [0,       R0_init_vak,   0,         0];
 
for k = 1:4
    I = Yarr{k}(:,2)*N;
    [maxI, idx] = max(I);
    tmax = tarr{k}(idx);
    % Ukupno zarazenih = svi koji su prosli kroz I tokom simulacije
    % = (1 - S_inf - inicijalno_vakcinisani) * N
    total_infected = (1 - Yarr{k}(end,1) - R0_inits(k)) * N;
    fprintf('%-12s %6.2f %10.0f %8.1f %14.0f\n', ...
            scen{k}, R0s(k), maxI, tmax, total_infected);
end
