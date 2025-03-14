% Figure 4 in the Results section - Simulate the impact of lipid-to-aSyn ratio on the mechanistic model of pure and lipidic aggregation

function simulate_asyn_agg_model_Figure4

clear all
close all
clc

%% General settings

% Parameters - lipidic scenario calibration

n = 2;
kn = 0; 
kmaxI = 3.09e-03;
Kp = 5.25;
kmaxII_p = 7.89e-03;
kmaxII_l = 0;
kmaxII_star = 0;
Ks = 2.13;
kc1 = 5.91e-01; 
kdis = 0.0024;
kc2 = 0.15*kc1;

kl2 = 7.92;
y = 457.5;
kon = 9.5*kl2; 
KD = 3.8e-01;
kbind_off = 1.63e-03;
kbind_on = kbind_off/KD;
gamma = 30;
Ms = 6000;

kl1_lip = 0;                  % lipidic aggregation (no pure aggregation)
kl1_lip_pure = 7.92;          % pure and lipidic aggregation

% Initial Conditions (Concentration \muM)
A_0         = 0;              % type-A oligomer concentration
B_0         = 0;              % type-B oligomer concentration
F_star_0    = 0;              % freshly generated fibril number concentration
Mf_star_0   = 0;              % the amount of monomers in freshly generated fibril
Fp_0        = 0;              % pure fibril number concentration
Mp_0        = 0;              % the amount of monomers in pure fibrils
Fl_0        = 0;              % lipidic fibril number concentration
Mfl_0       = 0;              % the amount of monomers in lipidic fibrils
Vf_b_0      = 0;              % concentration of monomer-coated SUVs bounded to fibrils
Vb_0        = 0;              % monomer-coated vesicle concentration

% The initial lipid-to-aSyn ratio - varying DMPS lipid concentration
DMPS_tot = [0 1 5 10 20 25 50 100 200 300 325 350 375 400 450 500 linspace(600,5000,20)];   % lipid concentrations
ncases = length(DMPS_tot);
asyn_tot = 50 * ones(1,ncases);                                                             % initial monomer concentration
init_ratio = DMPS_tot./asyn_tot;                                                            % initial lipid-to-aSyn ratios

% Color palette
load('colormap.mat','colors_lip','colors_lip_pure');

% ODE options
ODEoptions = odeset('RelTol',1e-5,'AbsTol',1e-6);

% Final time point for simulations (Time h)
tf = 120; 


% Initializing the aggregation metrics 

% Steady states
Mfp_st = zeros(1,ncases);
Mfl_st = zeros(1,ncases);
Mfp_tilde_st = zeros(1,ncases);
Mfl_tilde_st = zeros(1,ncases);
Mf_star_st = zeros(1,ncases);
Mftot_st = zeros(1,ncases);
Motot_st = zeros(1,ncases);
M_st = zeros(1,ncases);
Mb_st = zeros(1,ncases);
Mfb_st = zeros(1,ncases);

% Half times
Mfp_t12 = zeros(1,ncases);
Mfl_t12 = zeros(1,ncases);
Mfp_tilde_t12 = zeros(1,ncases);
Mfl_tilde_t12 = zeros(1,ncases);
Mftot_t12 = zeros(1,ncases);
Mf_star_t12 = zeros(1,ncases);

% Composition indexes
Icross = zeros(1,ncases);
Iself = zeros(1,ncases);
Ip = zeros(1,ncases);
Il = zeros(1,ncases);
Istar = zeros(1,ncases);
Ftot2_AUC = zeros(1,ncases);
Fp2_AUC = zeros(1,ncases);
Fp2_tilde_AUC = zeros(1,ncases);
Fl2_AUC = zeros(1,ncases);
Fl2_tilde_AUC = zeros(1,ncases);
Fstar2_AUC = zeros(1,ncases);

p_star = zeros(1,ncases);
pp = zeros(1,ncases);
pl = zeros(1,ncases);

% Olig abundance and peak time
A_peak_abundance = zeros(1,ncases);
A_peak_time = zeros(1,ncases);

B_peak_abundance = zeros(1,ncases);
B_peak_time = zeros(1,ncases);

% Olig propensity to aggregate - cumulative probability - integral of kc2B(t)/(kc2B(t) + kdA(t)) within the given time range
olig_prop = zeros(1,ncases);

% F* propensity to undergo one aggregation pathway or the other
F_star_prop = zeros(1,ncases);

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Lipidic aggregation only %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%

kl1 = kl1_lip;

% Define the parameter vector v
rates_lip = [n kn kmaxI Kp kmaxII_p Ks kc1 kdis kc2 kl1 kl2 y kon KD kbind_off kbind_on gamma Ms kmaxII_l kmaxII_star];
v = rates_lip;
colors = colors_lip;

% Parameter labels
v_labels = {'n' 'kn' 'kmaxI' 'Kp' 'kmaxII_p' 'Ks' 'kc1' 'kdis' 'kc2' 'kl1' 'kl2' 'y' 'kon' 'KD' 'kbind_off' 'kbind_on' 'gamma' 'Ms' 'kmaxII_l' 'kmaxII_star'};

% Display the table of parameter values
par_tab = array2table(v,'VariableNames',v_labels);
disp(par_tab)

% Aggregation pathways involved
if (kl1 == 0 && any(DMPS_tot))
    disp('Lipidic aggregation')
elseif (kl1 ~= 0 && any(DMPS_tot))
    disp('Pure protein and lipidic aggregation')  
elseif kl1 ~=0
    disp('Pure protein aggregation') 
elseif kl1 == 0
    disp('No aggregation')
end

% Deterministic simulations

for i = 1:ncases

    % Indexes for colors
    j = i;

    % Initial conditions for mtot and DMPStot (in the form of Vf_0)
    M_0 = asyn_tot(i);
    
    DMPSL_0 = DMPS_tot(i)/gamma;
    Vf_0 = gamma/Ms * DMPSL_0;

    % Initial condition vector
    x_init = [M_0, A_0, B_0, F_star_0, Mf_star_0, Fp_0, Mp_0, Fl_0, Mfl_0, Vf_b_0, Vf_0, Vb_0];

    [t,x] = ode23s(@(t,x) asyn_agg_model(t,x,v),[0 tf],x_init, ODEoptions);
    M = x(:,1);          % M     = the amount of free monomers in solution
    A = x(:,2);          % A     = type-A oligomer number concentration
    B = x(:,3);          % B     = type-B oligomer number concentration
    F_star = x(:,4);     % Fl    = freshly generated fibril number concentration
    Mf_star = x(:,5);    % Mfl   = the amount of monomers in freshly generated fibril
    Fp = x(:,6);         % Fp    = pure fibrils fibril number concentration
    Mfp = x(:,7);        % Mfp   = the amount of monomers in pure fibrils
    Fl = x(:,8);         % Fl    = lipidic fibril number concentration
    Mfl = x(:,9);        % Mfl   = the amount of monomers in lipidic fibrils
    Vf_b = x(:,10);      % Vf_b  = the amount of monomer-coated SUVs bounded to fibrils
    Vf = x(:,11);        % Vf    = the amount of free vesicles in solution
    Vb = x(:,12);        % Vb    = the amount of monomer-coated lipid vesicles in solution

    Mf_b = Ms/gamma * Vf_b;       % Mf_b       = the amount of monomers in SUVs bounded to fibrils
    DMPS_gamma = Ms/gamma * Vf;   % DMPS_gamma = Ms/gamma * Vf (DMPS_gamma = free vesicles in solution such that DMPS_free = gamma * DMPS_gamma (NB: DMPS_bound = Ms * Vb = gamma * Mb))
    M_b = Ms/gamma * Vb;          % Mb         = the amount of lipid-bounded monomers
    MA = n * A;
    MB = n * B;
    Mfp_tilde = Mfp - Mf_star;
    Mfl_tilde = Mfl - Mf_star;
    Mftot = Mfl + Mfp - Mf_star;  % Mftot = total fibril mass concentration
    Fp_tilde = Fp - F_star;
    Fl_tilde = Fl - F_star;
    Ftot = Fp + Fl - F_star;      % Ftot  = total fibril number concentration   
    mtot = M + M_b + MA + MB + Mf_b + Mfp + Mfl - Mf_star;

    % Plateau levels
    Mf_star_st(i) = Mf_star(end);
    Mfp_st(i) = Mfp(end);
    Mfl_st(i) = Mfl(end);
    Mfp_tilde_st(i) = Mfp_tilde(end);
    Mfl_tilde_st(i) = Mfl_tilde(end);
    Mftot_st(i) = Mftot(end);
    Motot_st(i) = MA(end) + MB(end);
    M_st(i) = M(end);
    Mb_st(i) = M_b(end);
    Mfb_st(i) = Mf_b(end);

    % Half times
    [val,indx_p] = min(abs(Mfp-Mfp(end)/2));
    Mfp_t12(i) = t(indx_p);

    [val,indx_l] = min(abs(Mfl-Mfl(end)/2));
    Mfl_t12(i) = t(indx_l);

    [val,indx_p_tilde] = min(abs(Mfp_tilde-Mfp_tilde(end)/2));
    Mfp_tilde_t12(i) = t(indx_p_tilde);

    [val,indx_l_tilde] = min(abs(Mfl_tilde-Mfl_tilde(end)/2));
    Mfl_tilde_t12(i) = t(indx_l_tilde);

    [val,indx_tot] = min(abs(Mftot-Mftot(end)/2));
    Mftot_t12(i) = t(indx_tot);

    [val,indx_star] = min(abs(Mf_star-Mf_star(end)/2));
    Mf_star_t12(i) = t(indx_star);

    % Ftot, Fp, Fl (^2) AUC
    time_int = t';
    Ftot_aux = (Ftot.^2)';
    Ftot2_AUC(i) = trapz(time_int,Ftot_aux,2);
    Fp_aux = (Fp.^2)';
    Fp2_AUC(i) = trapz(time_int,Fp_aux,2);
    Fp_tilde_aux = (Fp_tilde.^2)';
    Fp2_tilde_AUC(i) = trapz(time_int,Fp_tilde_aux,2);
    Fl_aux = (Fl.^2)';
    Fl2_AUC(i) = trapz(time_int,Fl_aux,2);
    Fl_tilde_aux = (Fl_tilde.^2)';
    Fl2_tilde_AUC(i) = trapz(time_int,Fl_tilde_aux,2);
    F_star_aux = (F_star.^2)';
    Fstar2_AUC(i) = trapz(time_int,F_star_aux,2);

    % LB composition indexes
    Ip(i) = Fp2_tilde_AUC(i)./Ftot2_AUC(i);
    Il(i) = Fl2_tilde_AUC(i)./Ftot2_AUC(i);
    Istar(i) = Fstar2_AUC(i)./Ftot2_AUC(i);
    Iself(i) = Ip(i) + Il(i) + Istar(i);
    Icross(i) = 1 - Iself(i);

    % Portion of pure, lipidic, newly generated fibrils (end point levels)
    p_star(i) = F_star(end)/Ftot(end);
    pp(i) = Fp_tilde(end)/Ftot(end);
    pl(i) = Fl_tilde(end)/Ftot(end);

    % Olig abundance and peak time
    [val,indx] = max(A);
    A_peak_time(i) = t(indx);
    A_peak_abundance(i) = val;

    [val,indx] = max(B);
    B_peak_time(i) = t(indx);
    B_peak_abundance(i) = val;

    % Olig propensity to aggregation
    time_int = t(2:end)';
    A_contrib = kdis * A(2:end);
    B_contrib = kc2 * B(2:end);
    olig_prop_int = B_contrib./(B_contrib + A_contrib);
    olig_prop_int = olig_prop_int';
    olig_prop(i) = trapz(time_int, olig_prop_int, 2);

    % F_star propensity
    time_int = t(2:end)';
    Vb_contrib = kon * Vb(2:end);
    M_contrib = kl1 * M(2:end);
    F_star_prop_int = Vb_contrib./(Vb_contrib + M_contrib);
    F_star_prop_int = F_star_prop_int';
    F_star_prop(i) = trapz(time_int, F_star_prop_int, 2)/tf;

    % Number concentrations
    figure(1)

    % Type-A oligomer number concentration
    subplot(6,2,1)
    plot(t, A,'Color', colors(j,:), 'LineWidth',2);
    hold on
    
    % Type-B oligomer number concentration
    subplot(6,2,3)
    plot(t, B,'Color', colors(j,:), 'LineWidth',2);
    hold on

    % Newly generated fibril number concentration
    subplot(6,2,5)
    plot(t, F_star,'Color', colors(j,:), 'LineWidth',2);
    hold on

    % Pure fibril number concentration
    subplot(6,2,7)
    plot(t, Fp_tilde,'Color', colors(j,:),'LineWidth',2);
    hold on

    % Lipidic fibril number concentration
    subplot(6,2,9)
    plot(t, Fl_tilde,'Color', colors(j,:),'LineWidth',2);
    hold on

    % Total fibril number concentratio
    subplot(6,2,11)
    plot(t, Ftot,'Color', colors(j,:),'LineWidth',2);
    hold on

    % Fibril mass concentrations
    figure(2)
    
    subplot(4,2,1)
    plot(t, Mf_star/mtot(1),'Color', colors(j,:),'LineWidth',2);
    hold on

    subplot(4,2,3)
    plot(t, Mfp_tilde/mtot(1),'Color', colors(j,:),'LineWidth',2);
    hold on

    subplot(4,2,5)
    plot(t, Mfl_tilde/mtot(1),'Color', colors(j,:),'LineWidth',2);
    hold on

    subplot(4,2,7)
    plot(t, Mftot/mtot(1),'Color', colors(j,:),'LineWidth',2);
    hold on

    % Aggregation metrics
    figure(3)
    
    % Abundance
    subplot(2,5,1)
    p1 = plot(init_ratio(i),A_peak_abundance(i),'o','MarkerSize',10,'MarkerFaceColor', colors(j,:),'Color',colors(j,:));
    hold on

    subplot(2,5,2)
    p2 = plot(init_ratio(i),B_peak_abundance(i),'o','MarkerSize',10,'MarkerFaceColor', colors(j,:),'Color',colors(j,:));
    hold on

    subplot(2,5,3)
    p3 = plot(init_ratio(i),Mfl_tilde_st(i)/asyn_tot(i),'o','MarkerSize',10,'MarkerFaceColor', colors(j,:),'Color',colors(j,:));
    hold on

    subplot(2,5,4)
    p4 = plot(init_ratio(i),Mftot_st(i)/asyn_tot(i),'o','MarkerSize',10,'MarkerFaceColor', colors(j,:),'Color',colors(j,:));
    hold on

    subplot(2,5,5)
    if i ==1
        plot(init_ratio(i),p_star(i),'square','MarkerSize',10,'MarkerEdgeColor',[0.5 0.5 0.5]);
        hold on
        plot(init_ratio(i),pp(i),'^','MarkerSize',10,'MarkerEdgeColor',[0.5 0.5 0.5]);
        plot(init_ratio(i),pl(i),'o','MarkerSize',10,'MarkerEdgeColor',[0.5 0.5 0.5]);
    end

    p5a = plot(init_ratio(i),p_star(i),'square','MarkerSize',10,'MarkerFaceColor', colors(j,:),'Color',colors(j,:));
    p5b = plot(init_ratio(i),pp(i),'^','MarkerSize',10,'MarkerFaceColor', colors(j,:),'Color',colors(j,:));
    p5c = plot(init_ratio(i),pl(i),'o','MarkerSize',10,'MarkerFaceColor', colors(j,:),'Color',colors(j,:));

    % Half-times
    subplot(2,5,6)
    p6 = plot(init_ratio(i),A_peak_time(i),'o','MarkerSize',10,'MarkerFaceColor', colors(j,:),'Color',colors(j,:));
    hold on

    subplot(2,5,7)
    p7 = plot(init_ratio(i),B_peak_time(i),'o','MarkerSize',10,'MarkerFaceColor', colors(j,:),'Color',colors(j,:));
    hold on

    subplot(2,5,8)
    p8 = plot(init_ratio(i),Mfl_tilde_t12(i),'o','MarkerSize',10,'MarkerFaceColor', colors(j,:),'Color',colors(j,:));
    hold on

    subplot(2,5,9)
    p9 = plot(init_ratio(i),Mftot_t12(i),'o','MarkerSize',10,'MarkerFaceColor', colors(j,:),'Color',colors(j,:));
    hold on
   
    subplot(2,5,10)
    p10 = plot(init_ratio(i),F_star_prop(i),'o','MarkerSize',10,'MarkerFaceColor', colors(j,:),'Color',colors(j,:));
    hold on

    p1.HandleVisibility = 'off';
    p2.HandleVisibility = 'off';
    p3.HandleVisibility = 'off';
    p4.HandleVisibility = 'off';
    p5a.HandleVisibility = 'off';
    p5b.HandleVisibility = 'off';
    p5c.HandleVisibility = 'off';
    p6.HandleVisibility = 'off';
    p7.HandleVisibility = 'off';
    p8.HandleVisibility = 'off';
    p9.HandleVisibility = 'off';
    p10.HandleVisibility = 'off';

    if i ~= 1
       p7.HandleVisibility = 'off';
    else
       p7.HandleVisibility = 'on'; 
    end

    % Other species in the system (monomers and lipidic vesicles)
    figure(4)
    
    subplot(4,2,1)
    plot(t, Vb,'Color', colors(j,:),'LineWidth',2);
    hold on

    subplot(4,2,3)
    plot(t, Vf,'Color', colors(j,:),'LineWidth',2);
    hold on

    subplot(4,2,5)
    plot(t, M,'Color', colors(j,:),'LineWidth',2);
    hold on

    subplot(4,2,7)
    plot(t, M_b,'Color', colors(j,:),'LineWidth',2);
    hold on

end

% Monomer concentration of aSyn species
figure(5)
subplot(1,2,1)
m_index = [M_st; Mb_st; Motot_st;  Mfb_st; Mf_star_st; Mfp_tilde_st; Mfl_tilde_st]./mtot(1);
m_index = m_index';
x_indices = 1:length(init_ratio);
b = bar(x_indices, m_index, 'stacked', 'BarWidth',1);

% Apply colors
color_palette = [5 59 108; 200 222 236; 11 90 82; 97 23 112; 90 157 204; 46 157 148; 133 91 165;]./255;  
for k = 1:size(m_index, 2)
    b(k).FaceColor = color_palette(k, :);
    b(k).EdgeColor = [1 1 1];
end


% LB composition indexes
figure(6)
subplot(1,2,1)
composition_index = [Istar; Ip; Il; Icross]';
x_indices = 1:length(init_ratio);
b = bar(x_indices, composition_index, 'stacked', 'BarWidth',1);

% Apply colors
color_palette = [90 157 204; 46 157 148; 133 91 165; 145 206 137]./255;
for k = 1:size(composition_index, 2)
    b(k).FaceColor = color_palette(k, :);
    b(k).EdgeColor = [1 1 1];
end

%% Figure settings

% Figure 1 settings
f1 = figure(1);
var_labels = {'$$A$$', '$$B$$', '$$F^*$$', '$$\tilde{F_p}$$','$$\tilde{F_l}$$','$$F_{TOT}$$'};
for i = 1:6
    j = 2*i-1;
    subplot(6,2,j)
    ylh = ylabel(var_labels{i},'FontWeight','bold','rotation',0,'HorizontalAlignment','right','VerticalAlignment','top','Interpreter', 'LaTeX');
    ylh.Units = 'pixels';
    ylh.Position = [-46.499999999999986,41.75003779423582,0];
    set(ylh,'rotation',0)
    if strcmp(var_labels{i}, '$$\tilde{F_p}$$')
        ylim([0 max(max(Fl_tilde), max(F_star))])
    end
    if j == 11
       xlabel('Time (hours)','Interpreter', 'LaTeX');
    end

     set(gca,'Fontsize',15)
end

ax1 = subplot(6,2,11);
colormap(ax1,colors);
c = colorbar;
title(c,'$\frac{DMPS_{TOT}}{m_{TOT}}$','Interpreter','latex')
caxis([min(init_ratio) max(init_ratio)])

% Figure 2 settings
f2 = figure(2);
var_labels = {'$$\frac{M_{F,*}}{m_{TOT}}$$', '$$\frac{\tilde{M}_{F,p}}{m_{TOT}}$$','$$\frac{\tilde{M}_{F,l}}{m_{TOT}}$$','$$\frac{M_{F,TOT}}{m_{TOT}}$$'};
for i = 1:4
    j = 2*i-1;
    subplot(4,2,j)
    ylh = ylabel(var_labels{i},'FontWeight','bold','rotation',0,'HorizontalAlignment','right','VerticalAlignment','middle','Interpreter', 'LaTeX');
    ylh.Units = 'pixels';
    if strcmp(var_labels{i}, '$$\frac{M_{F,TOT}}{m_{TOT}}$$')
        ylh.Position = [-37.49999999999997,41.75003779423582,0];
    else
        ylh.Position = [-48.49999999999997,41.75003779423582,0]; 
    end
    set(ylh,'rotation',0)
    if strcmp(var_labels{i}, '$$\frac{\tilde{M}_{F,p}}{m_{TOT}}$$')
       ylim([0 1])
    end
    if j == 7
       xlabel('Time (hours)','Interpreter', 'LaTeX');
    end
    set(gca,'Fontsize',15)
end

% Figure 3 settings
f3 = figure(3);
subplot(2,5,1)
title('$$A$$ peak level','Interpreter', 'LaTeX');
set(gca,'Fontsize',15)
xlabel('$$R$$','Interpreter', 'LaTeX');
ylabel('Concentration ($$\mu M$$)','Interpreter', 'LaTeX')
hold on

subplot(2,5,2)
title('$$B$$ peak level','Interpreter', 'LaTeX');
set(gca,'Fontsize',15)
xlabel('$$R$$','Interpreter', 'LaTeX');
ylabel('Concentration ($$\mu M$$)','Interpreter', 'LaTeX')
hold on

subplot(2,5,3)
title('$$\tilde{M}_{F,l}$$ relative abundance','Interpreter', 'LaTeX');
set(gca,'Fontsize',15)
xlabel('$$R$$','Interpreter', 'LaTeX');
ylabel('Normalized to $$m_{TOT}$$ (a.u.)','Interpreter', 'LaTeX');

subplot(2,5,4)
title('$$M_{F,TOT}$$ relative abundance','Interpreter', 'LaTeX');
set(gca,'Fontsize',15)
xlabel('$$R$$','Interpreter', 'LaTeX');
ylabel('Normalized to $$m_{TOT}$$ (a.u.)','Interpreter', 'LaTeX');

subplot(2,5,5)
title('Fibril relative abundances','Interpreter', 'LaTeX');
set(gca,'Fontsize',15)
xlabel('$$R$$','Interpreter', 'LaTeX');
ylabel('Normalized to $$F_{TOT}$$ (a.u.)','Interpreter', 'LaTeX');
hold on
legend('$$F^*$$','$$\tilde{F_p}$$', '$$\tilde{F_l}$$','Interpreter','Latex','Fontsize',12)

subplot(2,5,6)
title('$$A$$ peak time','Interpreter', 'LaTeX');
set(gca,'Fontsize',15)
xlabel('$$R$$','Interpreter', 'LaTeX');
ylabel('$$t_{1/2}$$ (h)','Interpreter', 'LaTeX');
set(gca, 'YScale', 'log')
hold on

subplot(2,5,7)
title('$B$ peak time','Interpreter', 'LaTeX');
set(gca,'Fontsize',15)
xlabel('$$R$$','Interpreter', 'LaTeX');
ylabel('$$t_{1/2}$$ (h)','Interpreter', 'LaTeX');
set(gca, 'YScale', 'log');
hold on

subplot(2,5,8)
title('$$\tilde{M}_{F,l}$$ half-time','Interpreter', 'LaTeX');
set(gca,'Fontsize',15)
xlabel('$$R$$','Interpreter', 'LaTeX');
ylabel('$$t_{1/2}$$ (h)','Interpreter', 'LaTeX');
set(gca, 'YScale', 'log')
hold on

subplot(2,5,9)
title('$$M_{F,TOT}$$ half-time','Interpreter', 'LaTeX');
set(gca,'Fontsize',15)
xlabel('$$R$$','Interpreter', 'LaTeX');
ylabel('$$t_{1/2}$$ (h)','Interpreter', 'LaTeX');
set(gca, 'YScale', 'log')
hold on

subplot(2,5,10)
title('$$F^*$$ propensity','Interpreter', 'LaTeX');
set(gca,'Fontsize',15)
xlabel('$$R$$','Interpreter', 'LaTeX');
ylabel('(a.u.)','Interpreter', 'LaTeX');
hold on

% Figure 4 settings
f4 = figure(4);
subplot(4,2,1)
ylh = ylabel('$$V_b$$','rotation',0,'HorizontalAlignment','right','VerticalAlignment','top','Interpreter','latex');
ylh.Units = 'pixels';
ylh.Position = [-46.499999999999986,41.75003779423582,0];
set(ylh,'rotation',0)
set(gca,'Fontsize',15)
hold on

subplot(4,2,3)
ylh = ylabel('$$V_f$$','rotation',0,'HorizontalAlignment','right','VerticalAlignment','top','Interpreter','latex');
ylh.Units = 'pixels';
ylh.Position = [-46.499999999999986,41.75003779423582,0];
set(ylh,'rotation',0)
set(gca,'Fontsize',15)
hold on

subplot(4,2,5)
ylh = ylabel('$$m$$','rotation',0,'HorizontalAlignment','right','VerticalAlignment','top','Interpreter','latex');
ylh.Units = 'pixels';
ylh.Position = [-46.499999999999986,41.75003779423582,0];
set(ylh,'rotation',0)
set(gca,'Fontsize',15)
hold on

subplot(4,2,7)
ylh = ylabel('$$m_b$$','rotation',0,'HorizontalAlignment','right','VerticalAlignment','top','Interpreter','latex');
ylh.Units = 'pixels';
ylh.Position = [-46.499999999999986,41.75003779423582,0];
set(ylh,'rotation',0)
set(gca,'Fontsize',15)
xlabel('Time (hours)','Interpreter', 'LaTeX');
hold on

% Figure 5 settings
f5 = figure(5);
set(gca,'Fontsize',15)
hold on
legend({'$$m$$', '$$m_b$$', '$$M_{O,TOT}$$', '$$m_{F,b}$$', '$$M_{F,*}$$', '$$\tilde{M}_{F,p}$$', '$$\tilde{M}_{F,l}$$'}, 'Location', 'northwest', 'Fontsize',16,'Interpreter', 'LaTex');
xlabel('$$R$$','Interpreter', 'LaTex','FontSize', 18);
formatted_xticks = arrayfun(@(x) sprintf('%.1f', x), init_ratio, 'UniformOutput', false);
xticks(x_indices);
xticklabels(formatted_xticks);
xtickangle(45);
ax = gca;
ax.XAxis.FontSize = 10;
ylim([0 1])
text(mean(xlim), 1.05 * max(ylim), 'Lipidic aggregation', 'FontSize', 14, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');

% Figure 6 settings
f6 = figure(6);
set(gca,'Fontsize',15)
hold on
legend({'$$I^*$$', '$$I^p$$', '$$I^l$$', '$$I^{cross}$$'}, 'Location', 'northwest', 'Fontsize',16,'Interpreter', 'LaTex');
xlabel('$$R$$','Interpreter', 'LaTex');
formatted_xticks = arrayfun(@(x) sprintf('%.1f', x), init_ratio, 'UniformOutput', false);
xticks(x_indices);
xticklabels(formatted_xticks);
xtickangle(45); 
ax = gca;
ax.XAxis.FontSize = 10;
text(mean(xlim), 1.05 * max(ylim), 'Lipidic aggregation', 'FontSize', 14, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Pure and lipidic aggregation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%

kl1 = kl1_lip_pure;

% Define the parameter vector
rates_pure_lip = [n kn kmaxI Kp kmaxII_p Ks kc1 kdis kc2 kl1 kl2 y kon KD kbind_off kbind_on gamma Ms kmaxII_l kmaxII_star];
v = rates_pure_lip;
colors = colors_lip_pure;

% Parameter labels
v_labels = {'n' 'kn' 'kmaxI' 'Kp' 'kmaxII_p' 'Ks' 'kc1' 'kdis' 'kc2' 'kl1' 'kl2' 'y' 'kon' 'KD' 'kbind_off' 'kbind_on' 'gamma' 'Ms' 'kmaxII_l' 'kmaxII_star'};

% Display the table of parameter values
par_tab = array2table(v,'VariableNames',v_labels);
disp(par_tab)

% Aggregation pathways involved
if (kl1 == 0 && any(DMPS_tot))
    disp('Lipidic aggregation')
elseif (kl1 ~= 0 && any(DMPS_tot))
    disp('Pure protein and lipidic aggregation')  
elseif kl1 ~=0
    disp('Pure protein aggregation') 
elseif kl1 == 0
    disp('No aggregation')
end

% Table of parameter values
par_tab = array2table(v,'VariableNames',v_labels);
disp(par_tab)

% Deterministic simulations

for i = 1:ncases
    
    % Indexes for colors
    j = i;

    % Initial conditions for mtot and DMPStot (in the form of Vf_0)
    M_0 = asyn_tot(i);
    
    DMPSL_0 = DMPS_tot(i)/gamma;
    Vf_0 = gamma/Ms * DMPSL_0;

    % Initial condition vector
    x_init = [M_0, A_0, B_0, F_star_0, Mf_star_0, Fp_0, Mp_0, Fl_0, Mfl_0, Vf_b_0, Vf_0, Vb_0];

    [t,x] = ode23s(@(t,x) asyn_agg_model(t,x,v),[0 tf],x_init, ODEoptions);
    M = x(:,1);          % M     = the amount of free monomers in solution
    A = x(:,2);          % A     = type-A oligomer number concentration
    B = x(:,3);          % B     = type-B oligomer number concentration
    F_star = x(:,4);     % Fl    = freshly generated fibril number concentration
    Mf_star = x(:,5);    % Mfl   = the amount of monomers in freshly generated fibril
    Fp = x(:,6);         % Fp    = pure fibrils fibril number concentration
    Mfp = x(:,7);        % Mfp   = the amount of monomers in pure fibrils
    Fl = x(:,8);         % Fl    = lipidic fibril number concentration
    Mfl = x(:,9);        % Mfl   = the amount of monomers in lipidic fibrils
    Vf_b = x(:,10);      % Vf_b  = the amount of monomer-coated SUVs bounded to fibrils
    Vf = x(:,11);        % Vf    = the amount of free vesicles in solution
    Vb = x(:,12);        % Vb    = the amount of monomer-coated lipid vesicles in solution

    Mf_b = Ms/gamma * Vf_b;       % Mf_b       = the amount of monomers in SUVs bounded to fibrils
    DMPS_gamma = Ms/gamma * Vf;   % DMPS_gamma = Ms/gamma * Vf (DMPS_gamma = free vesicles in solution such that DMPS_free = gamma * DMPS_gamma (NB: DMPS_bound = Ms * Vb = gamma * Mb))
    M_b = Ms/gamma * Vb;          % Mb         = the amount of lipid-bounded monomers
    MA = n * A;
    MB = n * B;
    Mfp_tilde = Mfp - Mf_star;
    Mfl_tilde = Mfl - Mf_star;
    Mftot = Mfl + Mfp - Mf_star;  % Mftot = total fibril mass concentration
    Fp_tilde = Fp - F_star;
    Fl_tilde = Fl - F_star;
    Ftot = Fp + Fl - F_star;      % Ftot  = total fibril number concentration   
    mtot = M + M_b + MA + MB + Mf_b + Mfp + Mfl - Mf_star;

    % Plateau levels
    Mf_star_st(i) = Mf_star(end);
    Mfp_st(i) = Mfp(end);
    Mfl_st(i) = Mfl(end);
    Mfp_tilde_st(i) = Mfp_tilde(end);
    Mfl_tilde_st(i) = Mfl_tilde(end);
    Mftot_st(i) = Mftot(end);
    Motot_st(i) = MA(end) + MB(end);
    M_st(i) = M(end);
    Mb_st(i) = M_b(end);
    Mfb_st(i) = Mf_b(end);

    % Half times
    [val,indx_p] = min(abs(Mfp-Mfp(end)/2));
    Mfp_t12(i) = t(indx_p);

    [val,indx_l] = min(abs(Mfl-Mfl(end)/2));
    Mfl_t12(i) = t(indx_l);

    [val,indx_p_tilde] = min(abs(Mfp_tilde-Mfp_tilde(end)/2));
    Mfp_tilde_t12(i) = t(indx_p_tilde);

    [val,indx_l_tilde] = min(abs(Mfl_tilde-Mfl_tilde(end)/2));
    Mfl_tilde_t12(i) = t(indx_l_tilde);

    [val,indx_tot] = min(abs(Mftot-Mftot(end)/2));
    Mftot_t12(i) = t(indx_tot);

    [val,indx_star] = min(abs(Mf_star-Mf_star(end)/2));
    Mf_star_t12(i) = t(indx_star);

    % Ftot, Fp, Fl (^2) AUC
    time_int = t';
    Ftot_aux = (Ftot.^2)';
    Ftot2_AUC(i) = trapz(time_int,Ftot_aux,2);
    Fp_aux = (Fp.^2)';
    Fp2_AUC(i) = trapz(time_int,Fp_aux,2);
    Fp_tilde_aux = (Fp_tilde.^2)';
    Fp2_tilde_AUC(i) = trapz(time_int,Fp_tilde_aux,2);
    Fl_aux = (Fl.^2)';
    Fl2_AUC(i) = trapz(time_int,Fl_aux,2);
    Fl_tilde_aux = (Fl_tilde.^2)';
    Fl2_tilde_AUC(i) = trapz(time_int,Fl_tilde_aux,2);
    F_star_aux = (F_star.^2)';
    Fstar2_AUC(i) = trapz(time_int,F_star_aux,2);

    % Composition indexes
    Ip(i) = Fp2_tilde_AUC(i)./Ftot2_AUC(i);
    Il(i) = Fl2_tilde_AUC(i)./Ftot2_AUC(i);
    Istar(i) = Fstar2_AUC(i)./Ftot2_AUC(i);
    Iself(i) = Ip(i) + Il(i) + Istar(i);
    Icross(i) = 1 - Iself(i);

    % Portion of pure, lipidic, newly generated fibrils (end point levels)
    p_star(i) = F_star(end)/Ftot(end);
    pp(i) = Fp_tilde(end)/Ftot(end);
    pl(i) = Fl_tilde(end)/Ftot(end);

    % Olig abundance and peak time
    [val,indx] = max(A);
    A_peak_time(i) = t(indx);
    A_peak_abundance(i) = val;

    [val,indx] = max(B);
    B_peak_time(i) = t(indx);
    B_peak_abundance(i) = val;

    % Olig propensity to aggregation
    F_star_AUC(i) = trapz(time_int,F_star_aux,2);
    time_int = t(2:end)';
    A_contrib = kdis * A(2:end);
    B_contrib = kc2 * B(2:end);
    olig_prop_int = B_contrib./(B_contrib + A_contrib);
    olig_prop_int = olig_prop_int';
    olig_prop(i) = trapz(time_int, olig_prop_int, 2);

    % F_star propensity
    time_int = t(2:end)';
    Vb_contrib = kon * Vb(2:end);
    M_contrib = kl1 * M(2:end);
    F_star_prop_int = Vb_contrib./(Vb_contrib + M_contrib);
    F_star_prop_int = F_star_prop_int';
    F_star_prop(i) = trapz(time_int, F_star_prop_int, 2)/tf;

    % Number concentrations
    figure(1)

    % Type-A oligomer number concentration
    subplot(6,2,2)
    plot(t, A,'Color', colors(j,:), 'LineWidth',2);
    hold on
    
    % Newly generated fibril number concentration
    subplot(6,2,4)
    plot(t, B,'Color', colors(j,:), 'LineWidth',2);
    hold on

    % Type-B oligomer number concentration
    subplot(6,2,6)
    plot(t, F_star,'Color', colors(j,:), 'LineWidth',2);
    hold on

    % Pure fibril number concentration
    subplot(6,2,8)
    plot(t, Fp_tilde,'Color', colors(j,:),'LineWidth',2);
    hold on

    % Lipidic fibril number concentration
    subplot(6,2,10)
    plot(t, Fl_tilde,'Color', colors(j,:),'LineWidth',2);
    hold on

    % Total fibril number concentratio
    subplot(6,2,12)
    plot(t, Ftot,'Color', colors(j,:),'LineWidth',2);
    hold on
   
    % Fibril mass concentrations
    figure(2)
    
    subplot(4,2,2)
    plot(t, Mf_star/mtot(1),'Color', colors(j,:),'LineWidth',2);
    hold on

    subplot(4,2,4)
    plot(t, Mfp_tilde/mtot(1),'Color', colors(j,:),'LineWidth',2);
    hold on

    subplot(4,2,6)
    plot(t, Mfl_tilde/mtot(1),'Color', colors(j,:),'LineWidth',2);
    hold on

    subplot(4,2,8)
    plot(t, Mftot/mtot(1),'Color', colors(j,:),'LineWidth',2);
    hold on


    % Aggregation metrics
    figure(3)
    
    % Abundance
    subplot(2,5,1)
    p1 = plot(init_ratio(i),A_peak_abundance(i),'o','MarkerSize',10,'MarkerFaceColor', colors(j,:),'Color',colors(j,:));
    hold on

    subplot(2,5,2)
    p2 = plot(init_ratio(i),B_peak_abundance(i),'o','MarkerSize',10,'MarkerFaceColor', colors(j,:),'Color',colors(j,:));
    hold on

    subplot(2,5,3)
    p3 = plot(init_ratio(i),Mfl_tilde_st(i)/asyn_tot(i),'o','MarkerSize',10,'MarkerFaceColor', colors(j,:),'Color',colors(j,:));
    hold on

    subplot(2,5,4)
    p4 = plot(init_ratio(i),Mftot_st(i)/asyn_tot(i),'o','MarkerSize',10,'MarkerFaceColor', colors(j,:),'Color',colors(j,:));
    hold on

    subplot(2,5,5)
    if i ==1
        plot(init_ratio(i),p_star(i),'square','MarkerSize',10,'MarkerEdgeColor',[0.5 0.5 0.5]);
        hold on
        plot(init_ratio(i),pp(i),'^','MarkerSize',10,'MarkerEdgeColor',[0.5 0.5 0.5]);
        plot(init_ratio(i),pl(i),'o','MarkerSize',10,'MarkerEdgeColor',[0.5 0.5 0.5]);
    end

    p5a = plot(init_ratio(i),p_star(i),'square','MarkerSize',10,'MarkerFaceColor', colors(j,:),'Color',colors(j,:));
    p5b = plot(init_ratio(i),pp(i),'^','MarkerSize',10,'MarkerFaceColor', colors(j,:),'Color',colors(j,:));
    p5c = plot(init_ratio(i),pl(i),'o','MarkerSize',10,'MarkerFaceColor', colors(j,:),'Color',colors(j,:));

    % Half-times
    subplot(2,5,6)
    p6 = plot(init_ratio(i),A_peak_time(i),'o','MarkerSize',10,'MarkerFaceColor', colors(j,:),'Color',colors(j,:));
    hold on

    subplot(2,5,7)
    p7 = plot(init_ratio(i),B_peak_time(i),'o','MarkerSize',10,'MarkerFaceColor', colors(j,:),'Color',colors(j,:));
    hold on

    subplot(2,5,8)
    p8 = plot(init_ratio(i),Mfl_tilde_t12(i),'o','MarkerSize',10,'MarkerFaceColor', colors(j,:),'Color',colors(j,:));
    hold on

    subplot(2,5,9)
    p9 = plot(init_ratio(i),Mftot_t12(i),'o','MarkerSize',10,'MarkerFaceColor', colors(j,:),'Color',colors(j,:));
    hold on
   
    subplot(2,5,10)
    p10 = plot(init_ratio(i),F_star_prop(i),'o','MarkerSize',10,'MarkerFaceColor', colors(j,:),'Color',colors(j,:));
    hold on

    p1.HandleVisibility = 'off';
    p2.HandleVisibility = 'off';
    p3.HandleVisibility = 'off';
    p4.HandleVisibility = 'off';
    p5a.HandleVisibility = 'off';
    p5b.HandleVisibility = 'off';
    p5c.HandleVisibility = 'off';
    p6.HandleVisibility = 'off';
    p7.HandleVisibility = 'off';
    p8.HandleVisibility = 'off';
    p9.HandleVisibility = 'off';
    p10.HandleVisibility = 'off';

    if i ~= 1
       p7.HandleVisibility = 'off';
    else
       p7.HandleVisibility = 'on'; 
    end

    % Other species in the system (monomers and lipidic vesicles)
    figure(4)
    
    subplot(4,2,2)
    plot(t, Vb,'Color', colors(j,:),'LineWidth',2);
    hold on

    subplot(4,2,4)
    plot(t, Vf,'Color', colors(j,:),'LineWidth',2);
    hold on

    subplot(4,2,6)
    plot(t, M,'Color', colors(j,:),'LineWidth',2);
    hold on

    subplot(4,2,8)
    plot(t, M_b,'Color', colors(j,:),'LineWidth',2);
    hold on

end

% Monomer concentration of aSyn species
figure(5)
sgtitle({'Relative mass of aSyn species'}, 'FontSize', 20, 'FontWeight', 'bold');
subplot(1,2,2)
title('Lipidic and pure aggregation')
m_index = [M_st; Mb_st; Motot_st;  Mfb_st; Mf_star_st; Mfp_tilde_st; Mfl_tilde_st]./mtot(1);
m_index = m_index';
x_indices = 1:length(init_ratio);
b = bar(x_indices, m_index, 'stacked', 'BarWidth',1);

% Apply colors
color_palette = [5 59 108; 200 222 236; 11 90 82; 97 23 112; 90 157 204; 46 157 148; 133 91 165;]./255;  
for k = 1:size(m_index, 2)
    b(k).FaceColor = color_palette(k, :);
    b(k).EdgeColor = [1 1 1];
end

% LB composition indexes
figure(6)
sgtitle({'Lewy body composition indexes'}, 'FontSize', 20, 'FontWeight', 'bold');
subplot(1,2,2)
composition_index = [Istar; Ip; Il; Icross]';
x_indices = 1:length(init_ratio);
b = bar(x_indices, composition_index, 'stacked', 'BarWidth',1);

% Apply colors
color_palette = [90 157 204; 46 157 148; 133 91 165; 145 206 137]./255;
for k = 1:size(composition_index, 2)
    b(k).FaceColor = color_palette(k, :);
    b(k).EdgeColor = [1 1 1];
end

%% Figure settings for second columns

% Figure 1 settings
f1 = figure(1);
for i = 1:6
    j = 2*i;
    subplot(6,2,j)
    if j == 12
       xlabel('Time (hours)','Interpreter', 'LaTeX');
    end
    set(gca,'Fontsize',15)
end 
ax2 = subplot(6,2,12);
colormap(ax2,colors);
c = colorbar;
title(c,'$\frac{DMPS_{TOT}}{m_{TOT}}$','Interpreter','latex')
caxis([min(init_ratio) max(init_ratio)])

% Figure 2 settings
f2 = figure(2);
for i = 1:4
    j = 2*i;
    subplot(4,2,j)
    if j == 8
       xlabel('Time (hours)','Interpreter', 'LaTeX');
    end
    set(gca,'Fontsize',15)
end

% Figure 3 settings
f3 = figure(3);
subplot(2,5,5)
title('Fibril relative abundances', 'Interpreter','LaTex');
set(gca,'Fontsize',15)
xlabel('$$R$$','Interpreter','LaTex');
ylabel('Normalized to $$F_{TOT}$$ (a.u.)','Interpreter','LaTex');
hold on
legend('$$F^*$$','$$\tilde{F_p}$$', '$$\tilde{F_l}$$','Interpreter','Latex','Fontsize',12)
subplot(2,5,7)
legend('Lipidic aggregation','Lipidic & pure aggregation');

% Figure 4 settings
f4 = figure(4);
for i = 1:8
    subplot(4,2,i)
    if i == 8 
        xlabel('Time (hours)','Interpreter', 'LaTeX');
        hold on
    end
    set(gca,'Fontsize',15)
end

% Figure 5 settings
f5 = figure(5);
set(gca,'Fontsize',15)
hold on
legend({'$$m$$', '$$m_b$$', '$$M_{O,TOT}$$', '$$m_{F,b}$$', '$$M_{F,*}$$', '$$\tilde{M}_{F,p}$$', '$$\tilde{M}_{F,l}$$'}, 'Location', 'northwest', 'Fontsize',16,'Interpreter', 'LaTex');
xlabel('$$R$$','Interpreter', 'LaTex','FontSize',18);
formatted_xticks = arrayfun(@(x) sprintf('%.1f', x), init_ratio, 'UniformOutput', false);
xticks(x_indices);
xticklabels(formatted_xticks);
xtickangle(45);
ax = gca;
ax.XAxis.FontSize = 10;
ylim([0 1])
text(mean(xlim), 1.05 * max(ylim), 'Lipidic and pure aggregation', 'FontSize', 14, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');

% Figure 6 settings
f6 = figure(6);
set(gca,'Fontsize',15)
hold on
legend({'$$I^*$$', '$$I^p$$', '$$I^l$$', '$$I^{cross}$$'}, 'Location', 'northwest', 'Fontsize',16,'Interpreter', 'LaTex');
xlabel('$$R$$','Interpreter', 'LaTex','FontSize',18);
formatted_xticks = arrayfun(@(x) sprintf('%.1f', x), init_ratio, 'UniformOutput', false);
xticks(x_indices);
xticklabels(formatted_xticks);
xtickangle(45); 
ax = gca;
ax.XAxis.FontSize = 10;
text(mean(xlim), 1.05 * max(ylim), 'Lipidic and pure aggregation', 'FontSize', 14, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');

%% Final figure settings

% Figure 1 settings
figure(1)
for i = 1:12
    subplot(6,2,i)
    xlim([0 tf])
end

% Figure 2 settings
figure(2)
for i = 1:8
subplot(4,2,i)
xlim([0 tf])
    if i == 1 || i== 2
        ylim([0 0.045])
    end
    if i == 3 || i== 4
        ylim([0 1])
    end
    if i == 5 || i== 6
        ylim([0 0.65])
    end
    if i == 7 || i== 8
        ylim([0 1])
    end

end

% Figure 3 settings
figure(3)
for i = 1:10
    subplot(2,5,i)
    xlim([0 max(init_ratio)])
end

% Figure 4 settings
figure(4)
for i = 1:8
    subplot(4,2,i)
    xlim([0 tf])
    if i == 1 || i== 2
        ylim([0 0.3])
    end
    if i == 3 || i== 4
        ylim([0 1])
    end
    if i == 5 || i== 6
        ylim([0 55])
    end
    if i == 7 || i== 8
        ylim([0 50])
    end
end

%% Set the figure size
set(f1,'Position', [440,70,773,727]);
set(f2,'Position', [440,70,773,727]);
set(f3,'Units','centimeters','Position',[0,7.337777777777777,50.8,20.77861111111111]);
set(f4,'Position', [440,70,773,727]);
set(f5,'Position', [48,367,1380,430]);
set(f6,'Position', [48,367,1380,430]);

% %% Save the figures
% % print(f1, 'Figure4_1_PDpaper.png', '-dpng', '-r300');
% % print(f2, 'Figure4_2_PDpaper.png', '-dpng', '-r300');
% % print(f3, 'Figure4_3_PDpaper.png', '-dpng', '-r300');
% % print(f4, 'Figure4_4_PDpaper.png', '-dpng', '-r300');
% print(f5, 'Figure4_5_PDpaper.png', '-dpng', '-r300');
% print(f6, 'Figure4_6_PDpaper.png', '-dpng', '-r300');

end