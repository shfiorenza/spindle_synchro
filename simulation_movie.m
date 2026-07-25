clear variables; 

perturb_type = 3; 
perturb_labels = ["Constant L, increased F"; "Const L, decreased F"; "Increased L, const F"; "Decreased L, const F"];
perturb_movieNames = ["mov_S1"; "mov_S2"; "mov_S3"; "mov_S4"];

make_video = true;

movie_duration = 20;

output_movie_name = perturb_movieNames(perturb_type);

t_perturb1 = 20;
t_perturb2 = 40;
% graphical parameters -- changes appearance of model schematic 
flux_inv = 1;
r_centro = 1.5;
r_kc = 0.25;
y_kc = r_centro - 2*r_kc;
mt_height = 0.5;
line_weight = 1.5;
flux_weight = 1.5;
mag_scale = 0.4;
mt_color = [0.9 0.9 0.9];

t_run = 60; %600                 % min
n_datapoints = 600; %6000;
t_pickup = t_run / n_datapoints;

if make_video
    v = VideoWriter(output_movie_name, 'MPEG-4');
    v.FrameRate = (t_run / t_pickup) / movie_duration;
    v.Quality = 100; %25;
    open(v);
end

params = [26, 180, 0.678, 9.95, 1.16, 2.3, 9.4, 20, 310, 995, 3.34, 1.0, 0.706];
i_param = 1;
%fixed constants
fixed_overlap = true;
dynamic_equil = false;
overlap_length = 6.55; 
eta = 33.3;                 % pN*min/um
N_kfl = 1;                  % number of k fibers
N_kfr = 1;                  % number of k fibers
N_bfl = 1;                  % number of bridge fibers
N_bfr = 1;                  % number of bridge fibers
c_m_B = params(i_param, 1); %30; %params(i_param, 1);   % 1/um
c_m_K = params(i_param, 1);
f_0 = 5;                    % pN    
v_0 = 3; %4; %3; %4;%3;                    % um/min
D_kc = 0.25;                % um
R0_kc = 0.45;                   % um
kappa_0 = params(i_param, 2);  % pN/um  
kappa_m_stretch = 750;
kappa_m_compress = 5000;
v_0__plus_k_base = params(i_param, 3);   % um/min
v_0__plus_b_base = params(i_param, 4);     % um/min   
v_0__minus_base = params(i_param, 5);     % um/min %1.4*
l_0__plus_k = params(i_param, 6);          % um
l_0__plus_b = params(i_param, 7);           % um 1.1
l_0__minus = params(i_param, 8);          % um 0.
F_0__plus = params(i_param, 9);           % pN
F_0__minus = params(i_param, 10);        % pN
A_pk = params(i_param, 11); 
A_pb = params(i_param, 12); 
A_m = params(i_param, 13);
dt = 1e-2;                  % min   

fig1 = figure;
set(fig1, 'Position', [50 50 1200 600]);

L_vs_t = NaN(n_datapoints, 1);
v_pkl_vs_t = NaN(n_datapoints, 1);
v_mkl_vs_t = NaN(n_datapoints, 1);

% initial conditions
early_exit = false;
eps_steady = 1e-6;
L_prev = 6.5;
l_kr_prev = 0.5*L_prev - D_kc - R0_kc/2.0;
l_kl_prev = 0.5*L_prev - D_kc - R0_kc/2.0;
if L_prev > 6.5 && dynamic_equil
    x_spb_l = -L_prev/2.0; %-3.25;
    x_spb_r = L_prev/2.0; %3.25;
    L = L_prev; %x_spb_r - x_spb_l;
    l_k_r = l_kr_prev; %0.5*L - D_kc - R0_kc/2.0;
    l_k_l = l_kl_prev; %0.5*L - D_kc - R0_kc/2.0;
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
delta_br = 0;
delta_bl = 0;
delta_kr = 0;
delta_kl = 0;
i_perturb1 = t_perturb1 / t_run * n_steps;
i_perturb2 = t_perturb2 / t_run * n_steps;

for i_step = 1 : 1 : n_steps
    % first, handle geometric parameters
    L = x_spb_r - x_spb_l;
    if L < 0  
        early_exit = true;
        %break;
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
    
    delta_br = mod(delta_br + v_br__minus*dt, flux_inv);
    delta_bl = mod(delta_bl - v_bl__minus*dt, flux_inv);
    delta_kr = mod(delta_kr + v_kr__minus*dt, flux_inv);
    delta_kl = mod(delta_kl - v_kl__minus*dt, flux_inv);
    
    if mod(double(i_step)*dt, t_pickup) == 0
        i_data = i_data + 1;
        % store data
        L_vs_t(i_data) = L;
        v_mkl_vs_t(i_data) = v_kl__minus;
        v_pkl_vs_t(i_data)= v_0__plus_kl*wt_KC__L;
        % plot figure
        clf;
        subplot(2, 2, [1,2])
        daspect([1,1,1])
        xlim([-10 10])
        ylim([-1.75 1.75])
        xticks([-10 -5 0 5 10]);
        yticks([]);

        x_kcl= x_spb_l + l_k_l + D_kc;
        x_kcr =x_spb_r - l_k_r - D_kc;


        if i_step == i_perturb1
            if perturb_type == 1
                A_pk = A_pk*1.3;
                A_pb = A_pb*1.3;
            elseif perturb_type == 2
                A_m = A_m*0.6;
            elseif perturb_type == 3
                 A_pk = A_pk*1.4;
                 A_pb = A_pb*1.4;
            elseif perturb_type == 4
                 A_m = A_m*1.2;
            else
                disp("Error. Invalid perturbation type.")
                return;
            end
        end
        if i_step == i_perturb2
            if perturb_type == 1
                A_m = A_m*1.3;
            elseif perturb_type == 2
                A_pk = A_pk*0.75;
                A_pb = A_pb*0.75;
            elseif perturb_type == 3
                A_m = A_m*0.8;
            elseif perturb_type == 4
                A_pk = A_pk*0.7;
                A_pb = A_pb*0.7;
            else
                disp("Error. Invalid perturbation type.")
                return;
            end
        end

        % centrosomes
        circle(x_spb_r, 0, r_centro, [1 1 1], line_weight); % green in fig: [21 96 130]/255
        circle(x_spb_l, 0, r_centro, [1 1 1], line_weight);

        % kinetochores
        circle(x_kcl, y_kc, r_kc, [1 1 1], line_weight); % blue in fig: [78 167 46]/255
        circle(x_kcr, y_kc, r_kc, [1 1 1], line_weight);
        % microtubules
        sc_p = 2.5;
        sc_m = 2.5;
        % right bridge
        rectangle('Position',[x_spb_r-l_b_r -mt_height/2 l_b_r mt_height], 'FaceColor', mt_color, 'LineWidth', line_weight);
        patch([x_spb_r-l_b_r x_spb_r-l_b_r x_spb_r-l_b_r+(sc_p-(A_pb - l_b_r / l_0__plus_b)) x_spb_r-l_b_r+(sc_p-(A_pb - l_b_r / l_0__plus_b))], [-mt_height/2 mt_height/2 mt_height/2 -mt_height/2], [0 0.447 0.7410], 'LineStyle', 'none', 'FaceAlpha', 'interp', 'FaceVertexAlphaData',[1 1 0.5 0.5]');
        patch([x_spb_r-sc_m*(l_b_r / l_0__minus + A_m) x_spb_r-sc_m*(l_b_r / l_0__minus + A_m) x_spb_r x_spb_r], [-mt_height/2 mt_height/2 mt_height/2 -mt_height/2], [0.85 0.325 0.0980], 'LineStyle', 'none', 'FaceAlpha', 'interp', 'FaceVertexAlphaData',[0.5 0.5 1 1]');
        % left bridge
        rectangle('Position',[x_spb_l -y_kc-mt_height/2 l_b_l mt_height], 'FaceColor', mt_color, 'LineWidth', line_weight);
        patch([x_spb_l+l_b_l x_spb_l+l_b_l x_spb_l+l_b_l-(sc_p-(A_pb - l_b_l / l_0__plus_b)) x_spb_l+l_b_l-(sc_p-(A_pb - l_b_l / l_0__plus_b))], [-y_kc-mt_height/2 -y_kc+mt_height/2 -y_kc+mt_height/2 -y_kc-mt_height/2], [0 0.447 0.7410], 'LineStyle', 'none', 'FaceAlpha', 'interp', 'FaceVertexAlphaData',[1 1 0.5 0.5]');
        patch([x_spb_l+sc_m*(l_b_l / l_0__minus + A_m) x_spb_l+sc_m*(l_b_l / l_0__minus + A_m) x_spb_l x_spb_l], [-y_kc-mt_height/2 -y_kc+mt_height/2 -y_kc+mt_height/2 -y_kc-mt_height/2], [0.85 0.325 0.0980], 'LineStyle', 'none', 'FaceAlpha', 'interp', 'FaceVertexAlphaData',[0.5 0.5 1 1]');
        % right k-fiber
        rectangle('Position',[x_spb_r-l_k_r y_kc-mt_height/2 l_k_r mt_height], 'FaceColor', mt_color, 'LineWidth', line_weight);
        patch([x_spb_r-l_k_r x_spb_r-l_k_r x_spb_r-l_k_r+(sc_p-(A_pk - l_k_r / l_0__plus_k)) x_spb_r-l_k_r+(sc_p-(A_pk - l_k_r / l_0__plus_k))], [y_kc-mt_height/2 y_kc+mt_height/2 y_kc+mt_height/2 y_kc-mt_height/2], [0 0.447 0.7410], 'LineStyle', 'none', 'FaceAlpha', 'interp', 'FaceVertexAlphaData',[1 1 0.5 0.5]');
        patch([x_spb_r-sc_m*(l_k_r / l_0__minus + A_m) x_spb_r-sc_m*(l_k_r / l_0__minus + A_m) x_spb_r x_spb_r], [y_kc-mt_height/2 y_kc+mt_height/2 y_kc+mt_height/2 y_kc-mt_height/2], [0.85 0.325 0.0980], 'LineStyle', 'none', 'FaceAlpha', 'interp', 'FaceVertexAlphaData',[0.5 0.5 1 1]');
        % left k-fiber
        rectangle('Position',[x_spb_l y_kc-mt_height/2 l_k_l mt_height], 'FaceColor', mt_color, 'LineWidth', line_weight);
        patch([x_spb_l+l_k_l x_spb_l+l_k_l x_spb_l+l_k_l-(sc_p-(A_pk - l_k_l / l_0__plus_k)) x_spb_l+l_k_l-(sc_p-(A_pk - l_k_l / l_0__plus_k))], [y_kc-mt_height/2 y_kc+mt_height/2 y_kc+mt_height/2 y_kc-mt_height/2], [0 0.447 0.7410], 'LineStyle', 'none', 'FaceAlpha', 'interp', 'FaceVertexAlphaData',[1 1 0.5 0.5]');
        patch([x_spb_l+sc_m*(l_k_l / l_0__minus + A_m) x_spb_l+sc_m*(l_k_l / l_0__minus + A_m) x_spb_l x_spb_l], [y_kc-mt_height/2 y_kc+mt_height/2 y_kc+mt_height/2 y_kc-mt_height/2], [0.85 0.325 0.0980], 'LineStyle', 'none', 'FaceAlpha', 'interp', 'FaceVertexAlphaData',[0.5 0.5 1 1]');

        % flux ticks for right bridge
        mag_br = (v_br__minus - mag_scale) / (v_0/2);
        for x_tick = x_spb_r-l_b_r:flux_inv:x_spb_r+2*flux_inv
            if x_tick+delta_br-mag_br < x_spb_r-l_b_r
                dxTick = (x_spb_r-l_b_r) - (x_tick+delta_br-mag_br);
                line([x_spb_r-l_b_r x_tick+delta_br], [(mt_height/2)/mag_br*dxTick mt_height/2], 'LineWidth', flux_weight, 'Color', 'k');
                line([x_spb_r-l_b_r x_tick+delta_br], [(-mt_height/2)/mag_br*dxTick -mt_height/2], 'LineWidth', flux_weight, 'Color', 'k');
            elseif x_tick+delta_br > x_spb_r
                if x_tick+delta_br-mag_br > x_spb_r
                    continue
                end
                dxTick = (x_spb_r) - (x_tick+delta_br);
                line([x_tick+delta_br-mag_br x_spb_r], [0 mt_height/2+(mt_height/2)/mag_br*dxTick], 'LineWidth', flux_weight, 'Color', 'k');
                line([x_tick+delta_br-mag_br x_spb_r], [0 -mt_height/2-(mt_height/2)/mag_br*dxTick], 'LineWidth', flux_weight, 'Color', 'k');
            else
                line([x_tick+delta_br-mag_br x_tick+delta_br], [0 mt_height/2], 'LineWidth', flux_weight, 'Color', 'k');
                line([x_tick+delta_br-mag_br x_tick+delta_br], [0 -mt_height/2], 'LineWidth', flux_weight, 'Color', 'k');
            end
        end
        % flux ticks for left bridge
        mag_bl = (v_bl__minus + mag_scale) / (v_0/2);
        for x_tick = x_spb_l+l_b_l:-flux_inv:x_spb_l-2*flux_inv
            if x_tick-delta_bl-mag_bl > x_spb_l+l_b_l
                dxTick = (x_tick-delta_bl-mag_bl) - (x_spb_l+l_b_l);
                line([x_spb_l+l_b_l x_tick-delta_bl], [-y_kc-(mt_height/2)/mag_bl*dxTick -y_kc+mt_height/2], 'LineWidth', flux_weight, 'Color', 'k');
                line([x_spb_l+l_b_l x_tick-delta_bl], [-y_kc+(mt_height/2)/mag_bl*dxTick -y_kc-mt_height/2], 'LineWidth', flux_weight, 'Color', 'k');
            elseif x_tick-delta_bl < x_spb_l
                if x_tick-delta_bl-mag_bl < x_spb_l
                    continue
                end
                dxTick = (x_spb_l) - (x_tick-delta_bl);
                line([x_tick-delta_bl-mag_bl x_spb_l], [-y_kc -y_kc-mt_height/2-(mt_height/2)/mag_bl*dxTick], 'LineWidth', flux_weight, 'Color', 'k');
                line([x_tick-delta_bl-mag_bl x_spb_l], [-y_kc -y_kc+mt_height/2+(mt_height/2)/mag_bl*dxTick], 'LineWidth', flux_weight, 'Color', 'k');
            else
                line([x_tick-delta_bl-mag_bl x_tick-delta_bl], [-y_kc -y_kc+mt_height/2], 'LineWidth', flux_weight, 'Color', 'k');
                line([x_tick-delta_bl-mag_bl x_tick-delta_bl], [-y_kc -y_kc-mt_height/2], 'LineWidth', flux_weight, 'Color', 'k');
            end
        end
        % flux ticks for right k-fiber
        mag_kr = (v_kr__minus - mag_scale) / (v_0/2);
        for x_tick = x_spb_r-l_k_r:flux_inv:x_spb_r+2*flux_inv
            if x_tick+delta_kr-mag_kr < x_spb_r-l_k_r
                dxTick = (x_spb_r-l_k_r) - (x_tick+delta_kr-mag_kr);
                line([x_spb_r-l_k_r x_tick+delta_kr], [y_kc+(mt_height/2)/mag_kr*dxTick y_kc+mt_height/2], 'LineWidth', flux_weight, 'Color', 'k');
                line([x_spb_r-l_k_r x_tick+delta_kr], [y_kc-(mt_height/2)/mag_kr*dxTick y_kc-mt_height/2], 'LineWidth', flux_weight, 'Color', 'k');
            elseif x_tick+delta_kr > x_spb_r
                if x_tick+delta_kr-mag_kr > x_spb_r
                    continue
                end
                dxTick = (x_spb_r) - (x_tick+delta_kr);
                line([x_tick+delta_kr-mag_kr x_spb_r], [y_kc y_kc+mt_height/2+(mt_height/2)/mag_kr*dxTick], 'LineWidth', flux_weight, 'Color', 'k');
                line([x_tick+delta_kr-mag_kr x_spb_r], [y_kc y_kc-mt_height/2-(mt_height/2)/mag_kr*dxTick], 'LineWidth', flux_weight, 'Color', 'k');
            else
                line([x_tick+delta_kr-mag_kr x_tick+delta_kr], [y_kc y_kc+mt_height/2], 'LineWidth', flux_weight, 'Color', 'k');
                line([x_tick+delta_kr-mag_kr x_tick+delta_kr], [y_kc y_kc-mt_height/2], 'LineWidth', flux_weight, 'Color', 'k');
            end
        end
        % flux ticks for left k-fiber
        mag_kl = (v_kl__minus + mag_scale) / (v_0/2);
        for x_tick = x_spb_l+l_k_l:-flux_inv:x_spb_l-2*flux_inv
            if x_tick-delta_kl-mag_kl > x_spb_l+l_k_l
                dxTick = (x_tick-delta_kl-mag_kl) - (x_spb_l+l_k_l);
                line([x_spb_l+l_k_l x_tick-delta_kl], [y_kc-(mt_height/2)/mag_kl*dxTick y_kc+mt_height/2], 'LineWidth', flux_weight, 'Color', 'k');
                line([x_spb_l+l_k_l x_tick-delta_kl], [y_kc+(mt_height/2)/mag_kl*dxTick y_kc-mt_height/2], 'LineWidth', flux_weight, 'Color', 'k');
            elseif x_tick-delta_kl < x_spb_l
                if x_tick-delta_kl-mag_kl < x_spb_l
                    continue
                end
                dxTick = (x_spb_l) - (x_tick-delta_kl);
                line([x_tick-delta_kl-mag_kl x_spb_l], [y_kc y_kc-mt_height/2-(mt_height/2)/mag_kl*dxTick], 'LineWidth', flux_weight, 'Color', 'k');
                line([x_tick-delta_kl-mag_kl x_spb_l], [y_kc y_kc+mt_height/2+(mt_height/2)/mag_kl*dxTick], 'LineWidth', flux_weight, 'Color', 'k');
            else
                line([x_tick-delta_kl-mag_kl x_tick-delta_kl], [y_kc y_kc+mt_height/2], 'LineWidth', flux_weight, 'Color', 'k');
                line([x_tick-delta_kl-mag_kl x_tick-delta_kl], [y_kc y_kc-mt_height/2], 'LineWidth', flux_weight, 'Color', 'k');
            end
        end

        str = sprintf('Time: %#.2f min', double(i_step) * dt);
        annotation('textbox', [0.125 0.65 .3 .3], 'String', str, 'FitBoxToText', 'on', 'FontSize', 18, 'LineStyle', 'none');
        xlabel("Position (μm)")
        set(gca, 'FontName', 'Arial')
        set(gca, 'FontSize', 16);
        set(gca,'box','off')

        subplot(2, 2, 3)
        plot(linspace(0, t_run, n_datapoints), L_vs_t, 'LineWidth', 2, 'Color', 'k')
        xlim([0 t_run])
        ylim([0 15])
        if perturb_type == 3
            ylim([0 20])
        end
        xlabel("Time (min)")
        ylabel("Spindle length (μm)");
        set(gca, 'FontName', 'Arial')
        set(gca, 'FontSize', 14);
        set(gca,'box','off')
        if i_step >= i_perturb1
            xline(t_perturb1, '--', 'LineWidth', 1.5);
            if perturb_type == 1
                str = ["Lose plus-end"; "proteins"];
            elseif perturb_type == 2
                str = ["Lose minus-end"; "proteins"];
            elseif perturb_type == 3
                str = ["Lose plus-end"; "proteins"];
            elseif perturb_type == 4
                str = ["Gain minus-end"; "proteins"];
            end
            annotation('textbox', [0.24 0 .3 .225], 'String', str, 'FitBoxToText', 'on', 'FontSize', 12, 'LineStyle', 'none');
        end
        if i_step >= i_perturb2
            xline(t_perturb2, '--', 'LineWidth', 1.5);
            if perturb_type == 1
                title("Rescue spindle length");
                str = ["Gain minus-end"; "proteins"];
            elseif perturb_type == 2
                title("Rescue spindle length");
                str = ["Gain plus-end"; "proteins"];
            elseif perturb_type == 3
                title("Increase spindle length")
                str = ["Lose minus-end"; "proteins"];
            elseif perturb_type == 4
                title("Decrease spindle length")
                str = ["Gain plus-end"; "proteins"];
            end
            annotation('textbox', [0.3505 0 .3 .225], 'String', str, 'FitBoxToText', 'on', 'FontSize', 12, 'LineStyle', 'none');
        end


        subplot(2, 2, 4)
        yyaxis left
        %plot(linspace(0, t_run, n_datapoints), v_kfr_vs_t, 'LineWidth', 2)
        plot(linspace(0, t_run, n_datapoints), v_pkl_vs_t, 'LineWidth', 2)
        hold on
        xlim([0 t_run])
        ylim([0 2]);
        ylabel("Plus-end poly. (μm/min)");
        %plot(linspace(0, t_run, n_datapoints), v_bfr_vs_t, 'LineWidth', 2);
        yyaxis right
        plot(linspace(0, t_run, n_datapoints), v_mkl_vs_t, 'LineWidth', 2)
        hold on
        set(gca, 'YDir','reverse')
        xlim([0 t_run])
        ylim([-2 0]);
        xlabel("Time (min)")
        ylabel("Minus-end depoly. (μm/min)");
        %legendLabel = {"Plus-end poly.", "Minus-end depoly."};
        %legend(legendLabel, 'location', 'southwest')
        set(gca, 'FontName', 'Arial')
        set(gca, 'FontSize', 14);
        set(gca,'box','off')
        if i_step >= i_perturb1
            xline(t_perturb1, '--', 'LineWidth', 1.5);
            y_pos = 0;
            if perturb_type == 1
                str = ["Lose plus-end"; "proteins"];
            elseif perturb_type == 2
                str = ["Lose minus-end"; "proteins"];
                y_pos = 0.2;   
            elseif perturb_type == 3
                str = ["Lose plus-end"; "proteins"];
            elseif perturb_type == 4
                str = ["Gain minus-end"; "proteins"];
            end
            annotation('textbox', [0.68 y_pos .3 .225], 'String', str, 'FitBoxToText', 'on', 'FontSize', 12, 'LineStyle', 'none');
        end
        if i_step >= i_perturb2
            xline(t_perturb2, '--', 'LineWidth', 1.5);
            y_pos = 0;
            if perturb_type == 1
                title("Increase poleward flux");
                str = ["Gain minus-end"; "proteins"];
            elseif perturb_type == 2
                title("Decrease poleward flux");
                str = ["Gain plus-end"; "proteins"];
                y_pos = 0.2;
            elseif perturb_type == 3
                title("Rescue poleward flux")
                str = ["Lose minus-end"; "proteins"];
            elseif perturb_type == 4
                title("Rescue poleward flux")
                str = ["Gain plus-end"; "proteins"];
            end
            annotation('textbox', [0.7905 y_pos .3 .225], 'String', str, 'FitBoxToText', 'on', 'FontSize', 12, 'LineStyle', 'none');
        end

        drawnow();
        if make_video
            writeVideo(v, getframe(fig1));
        end
        if early_exit
            n_datapoints_actual = i_data;
            t_run_actual = n_datapoints_actual * t_pickup;
            fprintf("Equil. time is %g\n", i_data*t_pickup);
            break;
        end
    end
end
if make_video
    close(v);
end

function h = circle(x, y, r, col, lw)
    d = r*2;
    px = x-r;
    py = y-r;
    h = rectangle('Position',[px py d d],'Curvature',[1,1], 'FaceColor', col, 'LineWidth', lw);
    %daspect([1,1,1])
end

