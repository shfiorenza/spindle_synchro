clear variables;

% different parameters scanned over
%% For Fig. 2
%{
var_dir = "v0m/"; // Also for Fig. 4
%var_label = "Baseline depoly. rate (um/min)";
var_label = "Depolymerization rate (um/min)";
var_seq = 0.1:0.1:5; 
val_nominal = 1.16; %1.16;
%xlims = [0.8*val_nominal 1.2*val_nominal]; %[0.0 1.9];
%}
%{
var_dir = "Apk/";
var_label = "K-fiber poly. l-dep const";
var_seq = [0:0.1:7.5];
val_nominal = 3.34;
%}
%{
var_dir = "Am/";
var_label = "Depoly. l-dep const";
var_seq = [0:0.1:7.5];
val_nominal = 0.706;
%}
%% for Fig. S3
%{
var_dir = "O/";
var_label = "Overlap (um)"; 
var_seq = 0:0.5:12;
val_nominal = 6.55;
xlims = [-0.5 12.5];
%}
%
var_dir = "cm/";
var_label = "Motor concentration (1/um)";
var_seq = [0:1:10,15:5:50,60:10:100,150:50:500];
xlims = [-25 525];
val_nominal = 26;
%}
%{
var_dir = "v0pk/";  
var_label = "Plus-end poly. rate (um/min)"; %"K-Fiber poly. (um/min)";
var_seq = 0.0:0.1:5;
val_nominal = 1.0;
xlims = [0.0 5.3];
%}

if ~exist('xlims', 'var')
    xlims = [min(var_seq) max(var_seq)];
end

folder = "spindle_scan/srebrenoParams/";
%folder = "spindle_scan/srebrenoParamsFineGrain_PaperData";
%folder = "spindle_scan/srebrenoParamsNewStats";

% For Fig. 2b
%dirs=[folder + "/noLDep_v0m1.5_v0pk1.5/v0pk/"; folder + "/noLDep_v0m1.5_v0pk1.5/v0m/"];
% For Fig. 2d
%dirs=[folder + "/noMinusLDep/Apk/"; folder + "/noMinusLDep/v0m/"];
% For Fig. 2f
%dirs=[folder + "/full/Apk/"; folder + "/full/Am/"];
% For Fig. 4
%dirs=["spindle_scan/srebrenoParamsAlt/full/" + var_dir];
% For Fig. S3
dirs=[folder + "/full/" + var_dir; folder + "/noMinusLDep/" + var_dir; folder + "/noPlusLDep/" + var_dir];

fontSize = 14;
line_styles = ["-"; "--"; ":"];
colors = [[0 0.4470 0.7410]; [0.8500 0.3250 0.0980]; [0.9290 0.6940 0.1250]; ...
    [0.4940 0.1840 0.5560]; [0.4660 0.6740 0.1880]; [0.3010 0.7450 0.9330]; ...
    [0.6350 0.0780 0.1840]];

n_datapoints = length(var_seq);
n_dirs = length(dirs);

L = zeros(n_dirs, n_datapoints);
lb = zeros(n_dirs,n_datapoints);
lk = zeros(n_dirs,n_datapoints);
v_kf = zeros(n_dirs ,n_datapoints);
v_bf = zeros(n_dirs, n_datapoints);
v_pk = zeros(n_dirs, n_datapoints);
v_pb = zeros(n_dirs, n_datapoints);
f_eq = zeros(n_dirs, n_datapoints);
f_bl = zeros(n_dirs, n_datapoints);
f_kl = zeros(n_dirs, n_datapoints);
f_bb = zeros(n_dirs, n_datapoints);
f_bk = zeros(n_dirs, n_datapoints);
f_kcl = zeros(n_dirs, n_datapoints);
for i_dir = 1:n_dirs
    for i_var=1:n_datapoints
        dir = dirs(i_dir);
        var = var_seq(i_var);
        name = sprintf("run_%.2f", var);
        output = readtable(dir + "stats_" + name + ".txt", 'Delimiter', '\t', NumHeaderLines=1);
        params = readtable(dir + "params_" + name + ".txt", 'Delimiter', '\t', NumHeaderLines=1);
        % variable name in output folder is different for older runs
        try
            L(i_dir, i_var) = output.Var2(find(strcmp('L', output.Var1)));
            f_eq(i_dir, i_var) = -1*output.Var2(find(strcmp('f_bl', output.Var1)));
            f_bl(i_dir, i_var) = -1*output.Var2(find(strcmp('f_bl', output.Var1)));
            f_kl(i_dir, i_var) = -1*output.Var2(find(strcmp('f_kl', output.Var1)));
            f_kcl(i_dir, i_var) = output.Var2(find(strcmp('f_kcl', output.Var1)));
            f_bb(i_dir, i_var) = output.Var2(find(strcmp('f_bbl', output.Var1)));
            f_bk(i_dir, i_var) = output.Var2(find(strcmp('f_bkl', output.Var1)));
        catch
            L(i_dir, i_var) = output.Var2(find(strcmp('separation', output.Var1)));
            f_eq(i_dir, i_var) = output.Var2(find(strcmp('f_eq', output.Var1)));
        end
        lb(i_dir,i_var) = output.Var2(find(strcmp('l_bl', output.Var1)));
        lk(i_dir,i_var) = output.Var2(find(strcmp('l_kl', output.Var1)));
        v_kf(i_dir,i_var) = -1*output.Var2(find(strcmp('v_kfl', output.Var1)));
        v_bf(i_dir,i_var) = -1*output.Var2(find(strcmp('v_bfl', output.Var1)));
        v_pk(i_dir,i_var) = output.Var2(find(strcmp('v_pkl', output.Var1)));
        v_pb(i_dir,i_var) = output.Var2(find(strcmp('v_pbl', output.Var1)));

    end
end
% get cut off lengths
Lmin = 6.5;
Lmax = 30;
i_start = 1 * ones(1, n_dirs);
i_end = length(var_seq) * ones(1, n_dirs);
i_last = length(var_seq);
for i_dir = 1 : n_dirs
    for i_var = 1 : length(var_seq)
        if L(i_dir, i_var) < Lmax && L(i_dir, i_var) > Lmin
            i_start(i_dir) = i_var;
            break;
        end
    end
    for i_var = length(var_seq) : -1 : 1
        if L(i_dir, i_var) < Lmax && L(i_dir, i_var) > Lmin
            i_end(i_dir) = i_var;
            break;
        end
    end
end
% ***length vs flux plot -- used for Fig. 2 ***
%{
fig_LvsF = figure('Position',[100 48 480 400], 'Units', 'normalized');
full = plot(v_kf(1, i_start(1):i_end(1)), L(1, i_start(1):i_end(1)), '.', 'MarkerSize', 32);%, LineWidth', 3);
hold on
for i_dir = 2 : n_dirs
    plot(v_kf(i_dir, i_start(i_dir):i_end(i_dir)), L(i_dir, i_start(i_dir):i_end(i_dir)), '.', 'MarkerSize', 32);%, LineWidth', 3);
end
ylabel("Spindle length (um)")
xlabel("Flux velocity (um/min)")
ylim([0 36]);
xlim([0.0 2]);
yticks([0 10 20 30]);
xticks([0 0.5 1.0 1.5 2.0]);
legendLabel = {'Plus-end perturbations'; 'Minus-end perturbations'};
legend(legendLabel);
box off
ax = gca;
ax.FontSize = fontSize;
ax.FontName = "Arial";
exportgraphics(fig_LvsF, "plot_LvsF_noLDep.pdf")
%}

% *** bridging vs k-fiber flux plots -- used to make Fig. 4 ***
%{
fig = figure('Position', [0 75 1440 360]);
i_dir = 1;
subplot(1, 3, 1);
plot(var_seq(i_start(i_dir):i_end(i_dir)), v_kf(i_dir, i_start(i_dir):i_end(i_dir)), 'LineWidth', 2, 'Color', colors(4, :), 'LineStyle', line_styles(i_dir, :));
hold on
plot(var_seq(i_start(i_dir):i_end(i_dir)), v_bf(i_dir, i_start(i_dir):i_end(i_dir)), 'LineWidth', 2, 'Color', colors(5, :), 'LineStyle', line_styles(i_dir, :));
plot(var_seq(i_end(i_dir):i_last), v_kf(i_dir, i_end(i_dir):i_last), 'LineWidth', 2, 'Color', [0.75 0.75 0.75], 'LineStyle', line_styles(i_dir, :));
plot(var_seq(i_end(i_dir):i_last), v_bf(i_dir, i_end(i_dir):i_last), 'LineWidth', 2, 'Color', [0.75 0.75 0.75], 'LineStyle', line_styles(i_dir, :));
plot(var_seq, var_seq, 'LineWidth', 1, 'Color', [0.5 0.5 0.55],'LineStyle','--')
xline(1.16, '--','LineWidth', 1, 'Color', [0.25 0.25 0.25]);
%yline(1.03, '--','LineWidth', 1, 'Color', colors(4, :));
%yline(1.5, '--','LineWidth', 1, 'Color', colors(5, :));
ylabel("Flux velocity (um/min)");
xlabel("Depoly. rate (um/min)");
yticks([0.5 1.0 1.5]);
xticks([0.8 1.2 1.6]);
ylim([0.4 1.65]);
xlim(xlims);
box off
ax = gca;
ax.FontSize = 14;
ax.FontName = 'Arial';
subplot(1, 3, 2)
plot(var_seq(i_start(i_dir):i_end(i_dir)), f_kl(i_dir, i_start(i_dir):i_end(i_dir)), 'LineWidth', 2, 'Color', colors(4, :), 'LineStyle', line_styles(i_dir, :))
hold on
plot(var_seq(i_start(i_dir):i_end(i_dir)), f_bl(i_dir, i_start(i_dir):i_end(i_dir)), 'LineWidth', 2, 'Color', colors(5, :), 'LineStyle', line_styles(i_dir, :))
plot(var_seq(i_end(i_dir):i_last), f_kl(i_dir, i_end(i_dir):i_last), 'LineWidth', 2, 'Color', [0.75 0.75 0.75], 'LineStyle', line_styles(i_dir, :))
plot(var_seq(i_end(i_dir):i_last), f_bl(i_dir, i_end(i_dir):i_last), 'LineWidth', 2, 'Color', [0.75 0.75 0.75], 'LineStyle', line_styles(i_dir, :))
%plot(var_seq(i_start(i_dir):i_end(i_dir)), f_kcl(i_dir, i_start(i_dir):i_end(i_dir)), 'LineWidth', 2, 'Color', colors(6, :), 'LineStyle', line_styles(i_dir, :))
%plot(var_seq(i_end(i_dir):i_last), f_kcl(i_dir, i_end(i_dir):i_last), 'LineWidth', 2, 'Color', [0.75 0.75 0.75], 'LineStyle', line_styles(i_dir, :))
%yline(1, '--','LineWidth', 1, 'Color', [0.25 0.25 0.25]);
xline(1.16, '--','LineWidth', 1, 'Color', [0.25 0.25 0.25]);
ylabel("Force on minus-end (pN)");
xlabel("Depoly. rate (um/min)");
yticks([-400 0 400]);
xticks([0.8 1.2 1.6]);
xlim(xlims);
%ylim([0.6 1.4]);
legendLabel = {"K-Fiber", "Bridge"};
legend(legendLabel, 'location', 'northeast')
box off
ax = gca;
ax.FontSize = 14;
ax.FontName = 'Arial';
subplot(1, 3, 3)
plot(var_seq(i_start(i_dir):i_end(i_dir)), lk(i_dir, i_start(i_dir):i_end(i_dir)), 'LineWidth', 2, 'Color', colors(4, :), 'LineStyle', line_styles(i_dir, :))
hold on
plot(var_seq(i_start(i_dir):i_end(i_dir)), lb(i_dir, i_start(i_dir):i_end(i_dir)), 'LineWidth', 2, 'Color', colors(5, :), 'LineStyle', line_styles(i_dir, :))
plot(var_seq(i_end(i_dir):i_last), lk(i_dir, i_end(i_dir):i_last), 'LineWidth', 2, 'Color', [0.75 0.75 0.75], 'LineStyle', line_styles(i_dir, :))
plot(var_seq(i_end(i_dir):i_last), lb(i_dir, i_end(i_dir):i_last), 'LineWidth', 2, 'Color', [0.75 0.75 0.75], 'LineStyle', line_styles(i_dir, :))
%yline(1, '--','LineWidth', 1, 'Color', [0.25 0.25 0.25]);
xline(1.16, '--','LineWidth', 1, 'Color', [0.25 0.25 0.25]);
ylabel("Microtubule length (um)");
xlabel("Depoly. rate (um/min)");
yticks([4 8 12]);
xticks([0.8 1.2 1.6]);
xlim(xlims);
ylim([2 13]);
box off
ax = gca;
ax.FontSize = 14;
ax.FontName = 'Arial';
%}

% *** isolated force plot  -- used to make Fig. S3 ***

fig = figure('Position', [50 75 540 540]);
plot(var_seq(i_start(1):i_end(1)), f_eq(1, i_start(1):i_end(1)), 'LineWidth', 3, 'Color', colors(3, :), 'LineStyle', line_styles(1, :));
hold on
for i_dir = 2 : n_dirs
    plot(var_seq(i_start(i_dir):i_end(i_dir)), f_eq(i_dir, i_start(i_dir):i_end(i_dir)), 'LineWidth', 3, 'Color', colors(3, :), 'LineStyle', line_styles(i_dir, :));
end
ylabel('Force (pN)','FontSize', fontSize);
xlabel(var_label,'FontSize', fontSize);
%ylim([-200 200]);
%ylim([-50 400]);
%ylim([-10 250]);
xlim(xlims);
ax = gca;
ax.FontSize = fontSize; 
ax.FontName = 'Arial';
box off
legendLabel = {'Full model', 'No minus-end l-dep', 'No plus-end l-dep'};
legend(legendLabel,'FontSize', fontSize, 'Location', 'southwest');%outside');
%}


