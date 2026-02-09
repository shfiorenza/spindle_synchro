clear variables; 


fixed_overlap = true;
dynamic_equil = true;
overlap_length = 6.55;       % um
eta = 33.3;                  % pN*min/um
N_kfl = 1;                   % number of k fibers
N_kfr = 1;                   % number of k fibers
N_bfl = 1;                   % number of bridge fibers
N_bfr = 1;                   % number of bridge fibers
c_m_B = 26;                  % 1/um
c_m_K = 26;                  % 1/um
f_0 = 5;                     % pN    
v_0 = 3;                     % um/min
D_kc = 0.25;                 % um
R0_kc = 0.45;                % um
kappa_0 = 180;               % pN/um  
kappa_m_stretch = 750;       % pN/um^2
kappa_m_compress = 5000;     % pN/um^2
v_0__plus_k_base = 0.678;    % um/min
v_0__plus_b_base = 9.95;     % um/min   
v_0__minus_base = 1.16;      % um/min
l_0__plus_k = 2.3;           % um
l_0__plus_b = 9.4;           % um 
l_0__minus = 20.0;           % um
F_0__plus = 310;             % pN
F_0__minus = 995;            % pN
A_pk = 3.34;                 % unitless
A_pb = 1.0;                  % unitless
A_m = 0.706;                 % unitless

dt = 1e-2;                  % min   
t_run = 600;                 % min
n_datapoints = 6000;


eps_steady = 1e-6;
L_prev = 6.5;
l_kr_prev = 0.5*L_prev - D_kc - R0_kc/2.0;
l_kl_prev = 0.5*L_prev - D_kc - R0_kc/2.0;


L_vs_t = zeros(n_datapoints, 1);
lkl_vs_t = zeros(n_datapoints, 1);
lkr_vs_t = zeros(n_datapoints, 1);
lbl_vs_t = zeros(n_datapoints, 1);
lbr_vs_t = zeros(n_datapoints, 1);
iKC_vs_t = zeros(n_datapoints, 1);
O_BB_vs_t = zeros(n_datapoints, 1);
O_KB_L_vs_t = zeros(n_datapoints, 1);
O_KB_R_vs_t = zeros(n_datapoints, 1);
v_SPBL_vs_t = zeros(n_datapoints, 1);
v_SPBR_vs_t = zeros(n_datapoints, 1);
v_kfl_vs_t = zeros(n_datapoints, 1);
v_kfr_vs_t = zeros(n_datapoints, 1);
v_bfl_vs_t = zeros(n_datapoints, 1);
v_bfr_vs_t = zeros(n_datapoints, 1);
v_mbl_vs_t = zeros(n_datapoints, 1);
v_mkl_vs_t = zeros(n_datapoints, 1);
v_pbl_vs_t = zeros(n_datapoints, 1);
v_pkl_vs_t = zeros(n_datapoints, 1);
wt_Fpkl_vs_t = zeros(n_datapoints, 1);
wt_Fpbl_vs_t = zeros(n_datapoints, 1);
wt_Fmkl_vs_t = zeros(n_datapoints, 1);
wt_Fmbl_vs_t = zeros(n_datapoints, 1);
wt_FmklActual_vs_t = zeros(n_datapoints, 1);
wt_FmblActual_vs_t = zeros(n_datapoints, 1);
wt_Lpkl_vs_t = zeros(n_datapoints, 1);
wt_Lpbl_vs_t = zeros(n_datapoints, 1);
wt_Lmkl_vs_t = zeros(n_datapoints, 1);
wt_Lmbl_vs_t = zeros(n_datapoints, 1);
x_spbl_vs_t = zeros(n_datapoints, 1);
x_spbr_vs_t = zeros(n_datapoints, 1);
x_kcl_vs_t = zeros(n_datapoints, 1);
x_kcr_vs_t = zeros(n_datapoints, 1); 
x_bl_vs_t = zeros(n_datapoints, 1);
x_br_vs_t = zeros(n_datapoints, 1);
f_kcr_vs_t = zeros(n_datapoints, 1);
f_kcl_vs_t = zeros(n_datapoints, 1);
f_bbr_vs_t = zeros(n_datapoints, 1);
f_bbl_vs_t = zeros(n_datapoints, 1);
f_kbr_vs_t = zeros(n_datapoints, 1);
f_kbl_vs_t = zeros(n_datapoints, 1);
f_bkr_vs_t = zeros(n_datapoints, 1);
f_bkl_vs_t = zeros(n_datapoints, 1);
f_kl_vs_t = zeros(n_datapoints, 1);
f_bl_vs_t = zeros(n_datapoints, 1); 
f_kr_vs_t = zeros(n_datapoints, 1);
f_br_vs_t = zeros(n_datapoints, 1);

% for scanning parameter space
%{
v0pk = v_0__plus_k_base;
v0m = v_0__minus_base;
l0pk = l_0__plus_k;
l0m = l_0__minus;
Apk = A_pk;
Am = A_m; 
%for factor = 0.05:0.01:1.95;
for factor_one = 0.05:0.01:1.95 %0.05:0.01:1.95 %0.05:0.01:1.95 %0.05:0.01:1.95%[-1:0.1:1]
    for factor_two = 0.05:0.01:1.95 %0.05:0.01:1.95%[-1:0.1:1]
    dir = "spindle_scan/PaperData/full/Apk_vs_Am/";
    %dir = "spindle_scan/PaperData/full/v0m/";
    name = sprintf("run_%.2f_%.2f", factor_one, factor_two)
    %name = sprintf("run_%.2f", factor)
    %v_0__plus_k_base = v0pk*factor;%_one;
    %v_0__minus_base =  v0m*factor;%_two;
    %l_0__plus_k = l0pk * factor_one; 
    %l_0__minus =  l0m * factor_two;
    A_pk = Apk * factor_one; 
    A_m =  Am * factor_two;
%}

% for synchro plot (used in to generate columns in Fig. 2)
%{
%factor_p = [2.0, 1.0, 0.5];
factor_p = [1.25, 1.0, 0.8];
fig_synchro = figure('Position',[100 48 480 640], 'Units', 'normalized');
s1 = subplot(2,1,1);
s2 = subplot(2,1,2);
hold on
for i_c = 1:3
    v_0__plus_k_base = v0pk * factor_p(i_c);
    %v_0__minus_base = v0m * factor_m(i_c);
%}

% initial conditions
early_exit = false;
if L_prev > 6.5 && dynamic_equil
    x_spb_l = -L_prev/2.0; 
    x_spb_r = L_prev/2.0;
    L = L_prev; 
    l_k_r = l_kr_prev; 
    l_k_l = l_kl_prev; 
else
    x_spb_l = -3.25;
    x_spb_r = 3.25;
    L = x_spb_r - x_spb_l;
    l_k_r = 0.5*L - D_kc - R0_kc/2.0;
    l_k_l = 0.5*L - D_kc - R0_kc/2.0;
end
l_b_r = 5*L/8; 
l_b_l = 5*L/8; 
if fixed_overlap
    l_b_r = (L+overlap_length) / 2.0;
    l_b_l = (L+overlap_length) / 2.0;
end
v_SPB__R = 0.0;
v_SPB__L = 0.0;
v_bl__minus = 0.0;
v_br__minus = 0.0;
v_kl__minus = 0.0;
v_kr__minus = 0.0;
% numerical integration routine 
n_steps = int32(t_run / dt);
t_pickup = t_run / n_datapoints;
i_data = 0;
n_datapoints_actual = n_datapoints;
t_run_actual = n_datapoints_actual * t_pickup;
for i_step = 1 : 1 : n_steps
    % first, handle geometric parameters
    L = x_spb_r - x_spb_l;
    if L < 0  
        early_exit = true;
    end
    O_BB = max(0, min(l_b_l + l_b_r - L, L)); 
    O_KB__R = max(0, min(l_k_r + l_b_l - L, l_k_r));
    O_KB__L = max(0, min(l_k_l + l_b_r - L, l_k_l));
    dx = L - (l_k_l + l_k_r + 2*D_kc + R0_kc);
    kappa = kappa_0 + (dx > 0)*kappa_m_stretch*dx - (dx < 0)*kappa_m_compress*dx;

    wt_KC__L = (1 + kappa*dx/(N_kfl*F_0__plus));
    wt_KC__R = (1 + kappa*dx/(N_kfr*F_0__plus));

    v_0__plus_bl = v_0__plus_b_base * (A_pb - l_b_l / l_0__plus_b);
    v_0__plus_br = v_0__plus_b_base * (A_pb - l_b_r / l_0__plus_b);
    v_0__plus_kl = v_0__plus_k_base * (A_pk - l_k_l / l_0__plus_k);
    v_0__plus_kr = v_0__plus_k_base * (A_pk - l_k_r / l_0__plus_k);

    v_0__minus_bl = v_0__minus_base * (l_b_l / l_0__minus + A_m);
    v_0__minus_br = v_0__minus_base * (l_b_r / l_0__minus + A_m);
    v_0__minus_kl = v_0__minus_base * (l_k_l / l_0__minus + A_m); 
    v_0__minus_kr = v_0__minus_base * (l_k_r / l_0__minus + A_m);
    
    % convention is v = [v_bl, v_br, v_kl, v_kr, v_spb_L, v_spb_R]
    %                   [v(1), v(2), v(3), v(4), v(5),    v(6)   ]
    func = @(v)[v(1) - (N_bfl*N_bfr*O_BB*c_m_B*f_0*((v(6) - v(5) + v(2))/v_0 - 1) ...
                     - F_0__minus*N_bfl ...
                     + N_bfl*N_kfr*O_KB__R*c_m_K*f_0*((v(6) - v(5) + v(4))/v_0 - 1)) ...
                     /((F_0__minus*N_bfl)/v_0__minus_bl + (N_bfl*N_bfr*O_BB*c_m_B*f_0)/v_0 + (N_bfl*N_kfr*O_KB__R*c_m_K*f_0)/v_0); ...
                v(2) - (F_0__minus*N_bfr ...
                     + N_bfl*N_bfr*O_BB*c_m_B*f_0*((v(5) - v(6) + v(1))/v_0 + 1) ....
                     + N_bfr*N_kfl*O_KB__L*c_m_K*f_0*((v(5) - v(6) + v(3))/v_0 + 1)) ...
                     /((F_0__minus*N_bfr)/v_0__minus_br + (N_bfl*N_bfr*O_BB*c_m_B*f_0)/v_0 + (N_bfr*N_kfl*O_KB__L*c_m_K*f_0)/v_0); ...           
                v(3) - (dx*kappa - F_0__minus*N_kfl ...
                     + N_bfr*N_kfl*O_KB__L*c_m_K*f_0*((v(6) - v(5) + v(2))/v_0 - 1)) ...
                     /((F_0__minus*N_kfl)/v_0__minus_kl + (N_bfr*N_kfl*O_KB__L*c_m_K*f_0)/v_0); ...
                v(4) - (F_0__minus*N_kfr - dx*kappa ...
                     + N_bfl*N_kfr*O_KB__R*c_m_K*f_0*((v(5) - v(6) + v(1))/v_0 + 1)) ...
                     /((F_0__minus*N_kfr)/v_0__minus_kr + (N_bfl*N_kfr*O_KB__R*c_m_K*f_0)/v_0); ...     
                v(5) - (dx*kappa + N_bfl*N_bfr*O_BB*c_m_B*f_0*((v(6) - v(1) + v(2))/v_0 - 1) ...
                     + N_bfr*N_kfl*O_KB__L*c_m_K*f_0*((v(6) + v(2) - v(3))/v_0 - 1) ...
                     + N_bfl*N_kfr*O_KB__R*c_m_K*f_0*((v(6) - v(1) + v(4))/v_0 - 1)) ...
                     /(eta + (N_bfl*N_bfr*O_BB*c_m_B*f_0)/v_0 + (N_bfr*N_kfl*O_KB__L*c_m_K*f_0)/v_0 + (N_bfl*N_kfr*O_KB__R*c_m_K*f_0)/v_0); ...
                v(6) - (N_bfl*N_bfr*O_BB*c_m_B*f_0*((v(5) + v(1) - v(2))/v_0 + 1) ...
                     - dx*kappa + N_bfr*N_kfl*O_KB__L*c_m_K*f_0*((v(5) - v(2) + v(3))/v_0 + 1) ...
                     + N_bfl*N_kfr*O_KB__R*c_m_K*f_0*((v(5) + v(1) - v(4))/v_0 + 1)) ...
                     /(eta + (N_bfl*N_bfr*O_BB*c_m_B*f_0)/v_0 + (N_bfr*N_kfl*O_KB__L*c_m_K*f_0)/v_0 + (N_bfl*N_kfr*O_KB__R*c_m_K*f_0)/v_0)];        
    % use prev values as initial guess
    v0 = [v_bl__minus, v_br__minus, v_kl__minus, v_kr__minus, v_SPB__L, v_SPB__R];
    options = optimset('Display','off', 'Algorithm', 'levenberg-marquardt');
    soln = fsolve(func, v0, options);

    v_bl__minus = soln(1);
    v_br__minus = soln(2);
    v_kl__minus = soln(3);
    v_kr__minus = soln(4);
    v_SPB__L = soln(5);
    v_SPB__R = soln(6);

    % update SPB position
    x_spb_l = x_spb_l + v_SPB__L * dt;
    x_spb_r = x_spb_r + v_SPB__R * dt;
    
    % calculate microtubule growth 
    dl_bl = v_0__plus_bl + v_bl__minus;
    dl_br = v_0__plus_br - v_br__minus;
    if fixed_overlap
        dl_bl = -v_SPB__L;
        dl_br = v_SPB__R;
    end
    dl_kl = v_0__plus_kl*wt_KC__L + v_kl__minus;
    dl_kr = v_0__plus_kr*wt_KC__R - v_kr__minus;
    
    % update microtubule lengths
    l_b_r = max(0, l_b_r + dl_br*dt);
    l_b_l = max(0, l_b_l + dl_bl*dt);
    l_k_r = max(0, l_k_r + dl_kr*dt);
    l_k_l = max(0, l_k_l + dl_kl*dt);

    % check equilibration status
    L_curr = x_spb_r - x_spb_l; 
    delta_L = L_curr - L_prev;
    if L_curr > 60 && delta_L > 0.0 && i_step > 100 && dynamic_equil 
        early_exit = true;
    end
    if abs(delta_L) < eps_steady && i_step > 100 && dynamic_equil 
        early_exit = true;
    end
    if dynamic_equil
        L_prev = L_curr;
        l_kr_prev = l_k_r;
        l_kl_prev = l_k_l;
    end
    % store data
    if mod(double(i_step)*dt, t_pickup) == 0
        i_data = i_data + 1;
        L_vs_t(i_data) = L;
        iKC_vs_t(i_data) = dx + R0_kc + D_kc;
        lbl_vs_t(i_data) = l_b_l;
        lbr_vs_t(i_data) = l_b_r;
        lkl_vs_t(i_data) = l_k_l;
        lkr_vs_t(i_data) = l_k_r;
        O_BB_vs_t(i_data) = O_BB;
        O_KB_L_vs_t(i_data) = O_KB__L;
        O_KB_R_vs_t(i_data) = O_KB__R;
        v_SPBR_vs_t(i_data) = v_SPB__R;
        v_SPBL_vs_t(i_data) = v_SPB__L;
        v_kfl_vs_t(i_data) = v_kl__minus + v_SPB__L;
        v_kfr_vs_t(i_data) = v_kr__minus + v_SPB__R;
        v_bfl_vs_t(i_data) = v_bl__minus + v_SPB__L;
        v_bfr_vs_t(i_data) = v_br__minus + v_SPB__R;
        v_mbl_vs_t(i_data) = v_bl__minus;
        v_mkl_vs_t(i_data) = v_kl__minus;
        v_pbl_vs_t(i_data) = v_0__plus_bl;
        if fixed_overlap
            v_pbl_vs_t(i_data) = -v_SPB__L - v_bl__minus;
        end
        v_pkl_vs_t(i_data)= v_0__plus_kl*wt_KC__L;
        x_spbr_vs_t(i_data) = x_spb_r;
        x_spbl_vs_t(i_data) = x_spb_l;
        x_kcl_vs_t(i_data) = x_spb_l + l_k_l + D_kc/2.0;
        x_kcr_vs_t(i_data) = x_spb_r - l_k_r - D_kc/2.0;
        x_bl_vs_t(i_data) = x_spb_l + l_b_l;
        x_br_vs_t(i_data) = x_spb_r - l_b_r;
        % calculate forces
        F_KCL = kappa*dx;
        F_BBL = -N_bfl*N_bfr*O_BB*c_m_B*f_0*((v_SPB__L - v_SPB__R + v_bl__minus - v_br__minus)/v_0 + 1);
        F_KBL = -N_bfr*N_kfl*O_KB__L*c_m_K*f_0*((v_SPB__L - v_SPB__R - v_br__minus + v_kl__minus)/v_0 + 1);
        F_BKL = -N_bfl*N_kfr*O_KB__R*c_m_K*f_0*((v_SPB__L - v_SPB__R + v_bl__minus - v_kr__minus)/v_0 + 1);
        F_KCR = -kappa*dx;
        F_BBR = N_bfl*N_bfr*O_BB*c_m_B*f_0*((v_SPB__L - v_SPB__R + v_bl__minus - v_br__minus)/v_0 + 1);
        F_KBR = N_bfl*N_kfr*O_KB__R*c_m_K*f_0*((v_SPB__L - v_SPB__R + v_bl__minus - v_kr__minus)/v_0 + 1);
        F_BKR = N_bfr*N_kfl*O_KB__L*c_m_K*f_0*((v_SPB__L - v_SPB__R - v_br__minus + v_kl__minus)/v_0 + 1);
        f_kcl_vs_t(i_data) = F_KCL;
        f_bbl_vs_t(i_data) = F_BBL;
        f_kbl_vs_t(i_data) = F_KBL;
        f_bkl_vs_t(i_data) = F_BKL;
        f_kcr_vs_t(i_data) = F_KCR;
        f_bbr_vs_t(i_data) = F_BBR;
        f_kbr_vs_t(i_data) = F_KBR;
        f_bkr_vs_t(i_data) = F_BKR;
        f_bl_vs_t(i_data) = F_BBL + F_BKL;
        f_kl_vs_t(i_data) = F_KBL + F_KCL;
        f_br_vs_t(i_data) = F_BBR + F_BKR;
        f_kr_vs_t(i_data) = F_KBR + F_KCR;
        wt_Fpkl_vs_t(i_data) = wt_KC__L;
        wt_Fpbl_vs_t(i_data) = 1.0;
        wt_Fmkl_vs_t(i_data) = (1 - (F_KBL + F_KCL)/F_0__minus);
        wt_Fmbl_vs_t(i_data) = (1 - (F_BBL + F_BKL)/F_0__minus);
        wt_FmklActual_vs_t(i_data) = -v_kl__minus / v_0__minus_kl;
        wt_FmblActual_vs_t(i_data) = -v_bl__minus / v_0__minus_bl;
        wt_Lpkl_vs_t(i_data) = v_0__plus_kl/v_0__plus_k_base;
        wt_Lpbl_vs_t(i_data) = v_0__plus_bl/v_0__plus_b_base;
        if fixed_overlap
            wt_Lpbl_vs_t(i_data) = 1.0;
        end
        wt_Lmkl_vs_t(i_data) = v_0__minus_kl/v_0__minus_base;
        wt_Lmbl_vs_t(i_data) = v_0__minus_bl/v_0__minus_base;
        if early_exit
            n_datapoints_actual = i_data;
            t_run_actual = n_datapoints_actual * t_pickup;
            fprintf("Equil. time is %g\n", i_data*t_pickup);
            break;
        end
    end
end
%}

% Create a pseudo-kymograph to show system evolution over time
kymo = figure('Position', [50 50 405 540]);
plot(x_spbl_vs_t(1:n_datapoints_actual),linspace(0, t_run_actual, n_datapoints_actual), 'Linewidth', 3);
hold on
plot(x_spbr_vs_t(1:n_datapoints_actual),linspace(0, t_run_actual, n_datapoints_actual), 'Linewidth', 3); 
plot(x_kcl_vs_t(1:n_datapoints_actual),linspace(0, t_run_actual, n_datapoints_actual), 'Linewidth', 2, 'LineStyle',':');
plot(x_kcr_vs_t(1:n_datapoints_actual),linspace(0, t_run_actual, n_datapoints_actual), 'Linewidth', 2, 'LineStyle',':');
plot(x_bl_vs_t(1:n_datapoints_actual),linspace(0, t_run_actual, n_datapoints_actual), 'Linewidth', 2, 'LineStyle','--');
plot(x_br_vs_t(1:n_datapoints_actual),linspace(0, t_run_actual, n_datapoints_actual), 'Linewidth', 2,  'LineStyle','--'); 
set(gca, 'YDir','reverse');
ylabel("Time (minutes)", 'FontSize', 14);
xlabel("Position (um)", 'FontSize', 14);
xlim([-10 10]);
ylim([0 t_run_actual]);

% Plot relevant simulation stats as they evolve over time
comp = figure('Position',[100 48 1120 550], 'Units', 'normalized');
subplot(2, 3, 1)
plot(linspace(0, t_run_actual, n_datapoints_actual), L_vs_t(1:n_datapoints_actual), 'LineWidth', 3);
hold on; 
plot(linspace(0, t_run_actual, n_datapoints_actual), iKC_vs_t(1:n_datapoints_actual), 'Linewidth', 2);
plot(linspace(0, t_run_actual, n_datapoints_actual), O_BB_vs_t(1:n_datapoints_actual), 'LineWidth', 2);
plot(linspace(0, t_run_actual, n_datapoints_actual), O_KB_R_vs_t(1:n_datapoints_actual), 'LineWidth', 2);
plot(linspace(0, t_run_actual, n_datapoints_actual), O_KB_L_vs_t(1:n_datapoints_actual), 'LineWidth', 2);
%ylim([0 20]);
xlim([0 t_run_actual]);
xlabel("Time (minutes)", FontSize=14);
ylabel("Distance (microns)", FontSize=14);
legendLabel = {"L", "iKC", "O_C", "O_R", "O_L"};
legend(legendLabel, 'location' , 'northeastoutside', 'FontSize', 11);
hsp1 = get(gca, 'Position');                 
set(gca, 'Position', [hsp1(1:2) 0.15 hsp1(4)]) 
subplot(2, 3, 2)
plot(linspace(0, t_run_actual, n_datapoints_actual), f_kl_vs_t(1:n_datapoints_actual), 'LineWidth', 2);
hold on
plot(linspace(0, t_run_actual, n_datapoints_actual), f_bl_vs_t(1:n_datapoints_actual), 'LineWidth', 2);
plot(linspace(0, t_run_actual, n_datapoints_actual), f_kcl_vs_t(1:n_datapoints_actual), ':', 'LineWidth', 2);
plot(linspace(0, t_run_actual, n_datapoints_actual), f_bbl_vs_t(1:n_datapoints_actual), ':', 'Linewidth', 2);
plot(linspace(0, t_run_actual, n_datapoints_actual), f_kbl_vs_t(1:n_datapoints_actual), '-.', 'Linewidth', 2);
plot(linspace(0, t_run_actual, n_datapoints_actual), f_bkl_vs_t(1:n_datapoints_actual), '--', 'Linewidth', 2);
%ylim([-175 175]);
xlim([0 t_run_actual]);
yline(0, ':', 'LineWidth', 1.5)
xlabel("Time (minutes)", FontSize=14);
ylabel("Force (pN)", FontSize=14);
legendLabel = {"F_{kl}", "F_{bl}", "F_{kcl}", "F_{bbl}", "F_{kbl}", "F_{bkl}"};%, "F_{xll}"};
legend(legendLabel, 'location' , 'northeastoutside', 'FontSize', 11);
hsp1 = get(gca, 'Position');                 
set(gca, 'Position', [hsp1(1:2) 0.15 hsp1(4)]) 
subplot(2, 3, 3)
plot(linspace(0, t_run_actual, n_datapoints_actual),lkl_vs_t(1:n_datapoints_actual), 'LineWidth', 2);
hold on
plot(linspace(0, t_run_actual, n_datapoints_actual),lbl_vs_t(1:n_datapoints_actual), 'LineWidth', 2);
yline(l_0__plus_k*A_pk, '--', 'LineWidth', 1.5, 'Color', [0 0.4470 0.7410]);
yline(l_0__plus_b*A_pb, '--', 'LineWidth', 1.5, 'Color', [0.8500 0.3250 0.0980]);
%ylim([0 15]);
xlim([0 t_run_actual]);
xlabel("Time (minutes)", FontSize=14);
ylabel("Length (microns)", FontSize=14);
legendLabel = {"l_{kl}", "l_{bl}", "l_{0k}", "l_{0b}"};
legend(legendLabel, 'location' , 'northeastoutside', 'FontSize', 11);
hsp1 = get(gca, 'Position');                 
set(gca, 'Position', [hsp1(1:2) 0.15 hsp1(4)]) 
subplot(2, 3, 4)
plot(linspace(0, t_run_actual, n_datapoints_actual),v_SPBL_vs_t(1:n_datapoints_actual), 'LineWidth', 2, 'Color', [0 0 0 ]);
hold on
plot(linspace(0, t_run_actual, n_datapoints_actual),v_kfl_vs_t(1:n_datapoints_actual), '--', 'Linewidth', 2, 'Color', [0 0.4470 0.7410]);
plot(linspace(0, t_run_actual, n_datapoints_actual),v_bfl_vs_t(1:n_datapoints_actual), '--', 'Linewidth', 2, 'Color', [0.8500 0.3250 0.0980]);
plot(linspace(0, t_run_actual, n_datapoints_actual),v_mkl_vs_t(1:n_datapoints_actual), 'Linewidth', 2, 'Color', [0 0.4470 0.7410]);
plot(linspace(0, t_run_actual, n_datapoints_actual),v_mbl_vs_t(1:n_datapoints_actual), 'Linewidth', 2, 'Color', [0.8500 0.3250 0.0980]);
plot(linspace(0, t_run_actual, n_datapoints_actual),v_pkl_vs_t(1:n_datapoints_actual), ':', 'Linewidth', 2, 'Color', [0 0.4470 0.7410]);
plot(linspace(0, t_run_actual, n_datapoints_actual),v_pbl_vs_t(1:n_datapoints_actual), ':', 'Linewidth', 2, 'Color', [0.8500 0.3250 0.0980]);
%ylim([-3 3]);
xlim([0 t_run_actual]);
xlabel("Time (minutes)", 'FontSize', 14);
ylabel("Velocity (micron/min)", 'FontSize', 14);
legendLabel = {"v^{l}_{c}", "v_{kfl}", "v_{bfl}", "v^{-}_{kl}", "v^{-}_{bl}", "v^{+}_{kl}", "v^{+}_{bl}"};
legend(legendLabel, 'location' , 'northeastoutside', 'FontSize', 11);
hsp1 = get(gca, 'Position');                 
set(gca, 'Position', [hsp1(1:2) 0.15 hsp1(4)]) 
subplot(2, 3, 5)
plot(linspace(0, t_run_actual, n_datapoints_actual), wt_Fmkl_vs_t(1:n_datapoints_actual), 'LineWidth', 2, 'Color', [0 0.4470 0.7410]);
hold on; 
plot(linspace(0, t_run_actual, n_datapoints_actual), wt_Fpkl_vs_t(1:n_datapoints_actual), ':', 'Linewidth', 2, 'Color', [0 0.4470 0.7410]);
plot(linspace(0, t_run_actual, n_datapoints_actual), wt_Fmbl_vs_t(1:n_datapoints_actual), 'LineWidth', 2, 'Color', [0.8500 0.3250 0.0980]);
plot(linspace(0, t_run_actual, n_datapoints_actual), wt_Fpbl_vs_t(1:n_datapoints_actual), ':', 'LineWidth', 2, 'Color', [0.8500 0.3250 0.0980]);
%plot(linspace(0, t_run, n_datapoints), wt_FmklActual_vs_t, ':', 'LineWidth', 2, 'Color', [0 0 0]);
%plot(linspace(0, t_run, n_datapoints), wt_FmblActual_vs_t, ':', 'LineWidth', 2, 'Color', [0 0 0]/2.0);
%yline(1, '--','LineWidth', 1, 'Color', [0.25 0.25 0.25]);
%ylim([0.8 1.2]);
xlim([0 t_run_actual]);
xlabel("Time (minutes)", FontSize=14);
ylabel("F-dep. factor (unitless)", FontSize=14);
legendLabel = {"wt F^{-}_{kl}", "wt F^{+}_{kl}", "wt F^{-}_{bl}", "wt F^{+}_{bl}"};
legend(legendLabel, 'location' , 'northeastoutside', 'FontSize', 11);
hsp1 = get(gca, 'Position');                 
set(gca, 'Position', [hsp1(1:2) 0.15 hsp1(4)]) 
subplot(2, 3, 6)
plot(linspace(0, t_run_actual, n_datapoints_actual), wt_Lmkl_vs_t(1:n_datapoints_actual), 'LineWidth', 2, 'Color', [0 0.4470 0.7410]);
hold on
plot(linspace(0, t_run_actual, n_datapoints_actual), wt_Lpkl_vs_t(1:n_datapoints_actual), ':', 'Linewidth', 2, 'Color', [0 0.4470 0.7410]);
plot(linspace(0, t_run_actual, n_datapoints_actual), wt_Lmbl_vs_t(1:n_datapoints_actual), 'LineWidth', 2, 'Color', [0.8500 0.3250 0.0980]);
plot(linspace(0, t_run_actual, n_datapoints_actual), wt_Lpbl_vs_t(1:n_datapoints_actual), ':', 'LineWidth', 2, 'Color', [0.8500 0.3250 0.0980]);
yline(1, '--','LineWidth', 1, 'Color', [0.25 0.25 0.25]);
%ylim([0.5 1.5]);
xlim([0 t_run_actual]);
xlabel("Time (minutes)", FontSize=14);
ylabel("L-dep. factor (unitless)", FontSize=14);
legendLabel = {"wt L^{-}_{kl}", "wt L^{+}_{kl}", "wt L^{-}_{bl}", "wt L^{+}_{bl}"};
legend(legendLabel, 'location' , 'northeastoutside', 'FontSize', 11);
hsp1 = get(gca, 'Position');                 
set(gca, 'Position', [hsp1(1:2) 0.15 hsp1(4)]) 
%}

%% synchro plot -- used to create columns in Fig. 2
%{
fontSize = 14;
%hold on
axes(s1);
hold on
plot(linspace(0, t_run_actual, n_datapoints_actual), L_vs_t(1:n_datapoints_actual), 'LineWidth', 3, 'Color', colors(i_c, :));
ylabel("Distance (um)", FontSize=fontSize);
ylim([0 15]);
xlim([0 30]);
if i_c == 3
    legendLabel = {'v_0^+ > v_0^-';'v_0^+ = v_0^-';'v_0^+ < v_0^-';};
    legend(legendLabel, 'location', 'northeast', 'AutoUpdate','off');
end
buffer = 0.15;
s1.Position(1) = buffer;
s1.Position(2) = 0.4 + buffer;
s1.Position(3) = 1.0 - 1.5*buffer;  
s1.Position(4) = 0.4;
ax = gca;
ax.FontSize = fontSize;
ax.FontName = "Arial";
axes(s2);
hold on
plot(linspace(0, t_run_actual, n_datapoints_actual),abs(v_pkl_vs_t(1:n_datapoints_actual)), '--', 'Linewidth', 2, 'Color', colors(i_c, :));%[0 0.4470 0.7410]);
plot(linspace(0, t_run_actual, n_datapoints_actual),abs(v_mkl_vs_t(1:n_datapoints_actual)), ':','Linewidth', 2, 'Color', colors(i_c, :));%[0 0.4470 0.7410]);
%plot(linspace(0, t_run_actual, n_datapoints_actual),abs(v_mbl_vs_t(1:n_datapoints_actual)), '--', 'Linewidth', 2, 'Color', 'g');%[0.8500 0.3250 0.0980]);
%plot(linspace(0, t_run_actual, n_datapoints_actual),abs(v_pbl_vs_t(1:n_datapoints_actual)), ':', 'Linewidth', 2, 'Color', 'g');%[0.8500 0.3250 0.0980]);
ylabel("Velocity (um/min)", 'FontSize', fontSize);
xlabel("Time (minutes)");
ylim([0.75 1.75]);
xlim([0 30]);
if i_c == 2
    legendLabel = {'v^+'; 'v^-'};
    legend(legendLabel, 'location', 'northeast', 'AutoUpdate','off');
end
s2.Position(1) = s1.Position(1);
s2.Position(2) = buffer/1.5;
s2.Position(3) = s1.Position(3);
s2.Position(4) = s1.Position(4);
ax = gca;
ax.FontSize = fontSize;
ax.FontName = "Arial";
hold off
drawnow;
if i_c == 3
   %exportgraphics(gcf, "plot_plusLDep_2.0_1.0_0.5.pdf")   
   exportgraphics(gcf, "plot_noLDep_1.25_1.0_0.8.pdf")
end
%}
%%
% output data
L = L_vs_t(n_datapoints_actual);
O_BB = O_BB_vs_t(n_datapoints_actual);
O_KB_L = O_KB_L_vs_t(n_datapoints_actual);
O_KB_R = O_KB_R_vs_t(n_datapoints_actual);
l_kl = lkl_vs_t(n_datapoints_actual);
l_kr = lkr_vs_t(n_datapoints_actual);
l_bl = lbl_vs_t(n_datapoints_actual);
l_br = lbr_vs_t(n_datapoints_actual);
KC_stretch = iKC_vs_t(n_datapoints_actual);
v_centro_L = v_SPBL_vs_t(n_datapoints_actual);
v_centro_R = v_SPBR_vs_t(n_datapoints_actual);
v_kfl = v_kfl_vs_t(n_datapoints_actual);
v_kfr = v_kfr_vs_t(n_datapoints_actual);
v_bfl = v_bfl_vs_t(n_datapoints_actual);
v_bfr = v_bfr_vs_t(n_datapoints_actual);
v_mkl = v_mkl_vs_t(n_datapoints_actual);
v_mbl = v_mbl_vs_t(n_datapoints_actual);
v_pkl = v_pkl_vs_t(n_datapoints_actual);
v_pbl = v_pbl_vs_t(n_datapoints_actual);
f_kcl = f_kcl_vs_t(n_datapoints_actual);
f_bbl = f_bbl_vs_t(n_datapoints_actual);
f_kbl = f_kbl_vs_t(n_datapoints_actual);
f_bkl = f_bkl_vs_t(n_datapoints_actual);
f_kcr = f_kcr_vs_t(n_datapoints_actual);
f_bbr = f_bbr_vs_t(n_datapoints_actual);
f_kbr = f_kbr_vs_t(n_datapoints_actual);
f_bkr = f_bkr_vs_t(n_datapoints_actual);
f_br = f_br_vs_t(n_datapoints_actual);
f_bl = f_bl_vs_t(n_datapoints_actual);
f_kr = f_kr_vs_t(n_datapoints_actual);
f_kl = f_kl_vs_t(n_datapoints_actual);
wt_f_pbl = wt_Fpbl_vs_t(n_datapoints_actual);
wt_f_pkl = wt_Fpkl_vs_t(n_datapoints_actual);
wt_f_mbl = wt_Fmbl_vs_t(n_datapoints_actual);
wt_f_mkl = wt_Fmkl_vs_t(n_datapoints_actual);
wt_f_mbl_act = wt_FmblActual_vs_t(n_datapoints_actual);
wt_f_mkl_act = wt_FmklActual_vs_t(n_datapoints_actual);
wt_l_pbl = wt_Lpbl_vs_t(n_datapoints_actual);
wt_l_pkl = wt_Lpkl_vs_t(n_datapoints_actual);
wt_l_mbl = wt_Lmbl_vs_t(n_datapoints_actual);
wt_l_mkl = wt_Lmkl_vs_t(n_datapoints_actual);
fprintf("%g | %g | %g | %g | %g  (%g)\n", L, O_BB, KC_stretch, v_kfr, v_bfr, f_bl)

% save data to file -- used for parameter scans
%{
t_window = 5;
n_window = ceil(n_datapoints * t_window/t_run);
v_avg = mean(v_SPBR_vs_t(1:n_window));
t_eq = n_datapoints_actual * t_pickup;
params = table(...
    fixed_overlap, dynamic_equil, overlap_length, eta, N_kfl, N_kfr, N_bfl, N_bfr, c_m_B, c_m_K, f_0, v_0, D_kc, R0_kc, ...
    kappa_0, kappa_m_stretch, kappa_m_compress, v_0__plus_k_base, v_0__plus_b_base, v_0__minus_base, l_0__plus_k, l_0__plus_b, l_0__minus, ...
    F_0__plus, F_0__minus, A_pk, A_pb, A_m, dt, t_run, n_datapoints, 'VariableNames', ...
    {'fixed_overlap', 'dynamic_equil', 'overlap_length', 'eta', 'N_kfl', 'N_kfr', 'N_bfl', 'N_bfr', 'c_m_B','c_m_K', 'f_0', 'v_0', 'D_kc', 'R0_kc', ...
    'kappa0', 'kappa_m_stretch', 'kappa_m_compress', 'v_0__plus_k_base', 'v_0__plus_b_base', 'v_0__minus_base', 'l_0__plus_k', 'l_0__plus_b', 'l_0__minus', ...
    'F_0__plus', 'F_0__minus', 'A_pk', 'A_pb', 'A_m', 'dt', 't_run', 'n_datapoints'});
writetable(rows2vars(params), dir + "params_" + name + ".txt", 'Delimiter', '\t');
output = table(...
    L, O_BB, O_KB_L, O_KB_R, l_kl, l_kr, l_bl, l_br, KC_stretch, ...
    v_centro_L, v_centro_R, v_kfl, v_kfr, v_bfl, v_bfr, v_mkl, v_mbl, v_pkl, v_pbl, ...
    f_kcl, f_bbl, f_kbl, f_bkl, f_kcr, f_bbr, f_kbr, f_bkr, f_br, f_bl, f_kr, f_kl, ...
    wt_f_pbl, wt_f_pkl, wt_f_mbl, wt_f_mkl, wt_f_mbl_act, wt_f_mkl_act, wt_l_pbl, wt_l_pkl, wt_l_mbl, wt_l_mkl,...
    v_avg, t_eq, 'VariableNames', ...
    {'L', 'O_BB', 'O_KB_L', 'O_KB_R', 'l_kl', 'l_kr', 'l_bl', 'l_br', 'KC_stretch', ...
    'v_c_L', 'v_c_R', 'v_kfl', 'v_kfr', 'v_bfl', 'v_bfr', 'v_mkl', 'v_mbl', 'v_pkl', 'v_pbl', ...
    'f_kcl', 'f_bbl', 'f_kbl', 'f_bkl', 'f_kcr', 'f_bbr', 'f_kbr', 'f_bkr', 'f_br', 'f_bl', 'f_kr', 'f_kl', ...
    'wt_f_pbl', 'wt_f_pkl', 'wt_f_mbl', 'wt_f_mkl', 'wt_f_mbl_act', 'wt_f_mkl_act', 'wt_l_pbl', 'wt_l_pkl', 'wt_l_mbl', 'wt_l_mkl', ...
    'v_avg','t_eq'});
writetable(rows2vars(output), dir + "stats_" + name + ".txt", 'Delimiter', '\t');
%saveas(kymo, dir + "kymo_" + name + ".png", 'png');
%saveas(comp, dir + "comp_" + name + ".png", 'png');
%close(kymo);
%close(comp);
    end
end
%}
