
clear variables;

param1 = 0.05:0.01:1.95;
param2 = 0.05:0.01:1.95;
%dir = "spindle_scan/srebrenoParamsFineGrain_PaperData/full/l0pk_vs_l0m/";
%dir = "spindle_scan/srebrenoParamsFineGrain_PaperData/full/Apk_vs_Am_unfixedO/";
%dir = "spindle_scan/srebrenoParamsFineGrain_PaperData/full/Apk_vs_Am/";
dir = "spindle_scan/srebrenoParamsFineGrain_PaperData/full/v0pk_vs_v0m/";
output_folder = "new_v"; %"new_figs_unfixedO";

i_dir = 3;
Lmin = 6.5; %2; %6.5;
Lmax = 30;
D_kc = 0.25;                    % um
R0_kc = 0.45;                   % um

n_i = length(param1);
n_j = length(param2);
Length = zeros(n_i, n_j);
Flux_k = zeros(n_i, n_j);
Flux_b = zeros(n_i, n_j);
Flux_norm = zeros(n_i, n_j);
Flux_overall = zeros(n_i, n_j);
Correlation = zeros(n_i, n_j);
Correlation_norm = zeros(n_i, n_j);
Correlation_overall = zeros(n_i, n_j);
KCstretch = zeros(n_i, n_j);
overlap = zeros(n_i, n_j); 
L0 = 11.8;
F0 = 1.03;
Fnorm0 = 0.7055;
Foverall0 = 1.14;
epsL = 0.0;
epsF = 0.0;
for i = 1 : 1 : n_i
    for j = 1 : 1 : n_j
        name = sprintf("run_%.2f_%.2f", param1(i), param2(j));
        output = readtable(dir + "stats_" + name + ".txt", 'Delimiter', '\t', NumHeaderLines=1);
        if i_dir == 1 || i_dir == 3
            Length(i, j) = output.Var2(find(strcmp('separation', output.Var1)));
        else
            Length(i, j) = output.Var2(find(strcmp('L', output.Var1)));
        end
        Flux_k(i, j) = -1*output.Var2(find(strcmp('v_kfl', output.Var1)));
        Flux_b(i, j) = -1*output.Var2(find(strcmp('v_bfl', output.Var1)));
        Flux_norm(i, j) = Flux_k(i, j) / Flux_b(i, j);
        Flux_overall(i, j) = (3*Flux_k(i, j) + Flux_b(i, j))/4;
        KCstretch(i, j) = output.Var2(find(strcmp('KC_stretch', output.Var1)));
        
        if i_dir == 1 || i_dir == 3
            overlap(i, j) = output.Var2(find(strcmp('overlap_C', output.Var1)));
        else
            overlap(i, j) = output.Var2(find(strcmp('O_BB', output.Var1)));
        end
        
        if i == n_i/2 + 1/2 && j == n_j/2 + 1/2
            continue;
        end
        deltaL = Length(i, j) - L0;
        deltaF = Flux_k(i, j) - F0;
        deltaFnorm = Flux_norm(i, j) - Fnorm0;
        deltaFoverall = Flux_overall(i, j) - Foverall0;
        Correlation(i, j) = get_correlation(deltaL, deltaF, epsL, epsF);
        Correlation_norm(i, j) = get_correlation(deltaL, deltaFnorm, epsL, epsF);
        Correlation_overall(i, j) = get_correlation(deltaL, deltaFoverall, epsL, epsF);
    end
end
%}

Length_bounded = Length;
Flux_k_bounded = Flux_k;
Flux_b_bounded = Flux_b;
Flux_norm_bounded = Flux_norm;
Flux_overall_bounded = Flux_overall;
Correlation_bounded = Correlation;
Correlation_norm_bounded = Correlation_norm;
Correlation_overall_bounded = Correlation_overall;
KCstretch_bounded = KCstretch - (R0_kc + D_kc);
overlap_bounded = overlap;

Length_bounded(Length_bounded < Lmin) = nan;
Length_bounded(Length_bounded > Lmax) = nan;
Flux_k_bounded(Length_bounded ~= Length_bounded) = nan;
Flux_b_bounded(Length_bounded ~= Length_bounded) = nan;
Flux_norm_bounded(Length_bounded ~= Length_bounded) = nan;
Flux_overall_bounded(Length_bounded ~= Length_bounded) = nan;
Correlation_bounded(Length_bounded ~= Length_bounded) = 0;
Correlation_norm_bounded(Length_bounded ~= Length_bounded) = 0;
Correlation_overall_bounded(Length_bounded ~= Length_bounded) = 0;
KCstretch_bounded(Length_bounded ~= Length_bounded) = nan;
overlap_bounded(Length_bounded ~= Length_bounded) = nan;

Length_cut = Length;
Flux_k_cut = Flux_k;
Flux_b_cut = Flux_b;
delta = 30;
for i = [1 : 96 - delta - 1, 96 + delta + 1 : 191]
    for j = [1 : 96 - delta - 1, 96 + delta + 1 : 191]
        Length_cut(i, j) = nan;
        Flux_k_cut(i, j) = nan;
        Flux_b_cut(i, j) = nan;
    end
end

% Phase contour plots
plot_phase(i_dir, true, Length_bounded, [6, floor(Lmin):3:Lmax], "Spindle length (um)", "length", output_folder);
plot_phase(i_dir, true, Flux_k_bounded, 0.0:0.2:2.0, "K-fiber flux (um/min)", "fluxK", output_folder);
plot_phase(i_dir, true, Flux_b_bounded, 0.0:0.2:2.0, "Bridging fiber flux (um/min)", "fluxB", output_folder);
plot_phase(i_dir, true, Flux_norm_bounded, 0.0:0.1:2.0, "Normalized flux (unitless)", "fluxNorm", output_folder);
plot_phase(i_dir, true, Flux_overall_bounded, 0.0:0.1:2.0, "Overall flux (um/min)", "fluxOverall", output_folder);
plot_phase(i_dir, true, KCstretch_bounded, -1.0:0.1:2.0, "Interkinetochore stretch (um)", "iKC", output_folder);
plot_phase(i_dir, true, overlap_bounded, -5:1:20, "Overlap length (um)", "O", output_folder);
plot_phase(i_dir, false, Correlation_bounded, [-1 -2/3 0 2/3 1], "Correlation domains", "corr", output_folder);
plot_phase(i_dir, false, Correlation_norm_bounded, [-1 -2/3 0 2/3 1], "Correlation domains (normalized flux)", "corrNorm", output_folder);
plot_phase(i_dir, false, Correlation_overall_bounded, [-1 -2/3 0 2/3 1], "Correlation domains (overall flux)", "corrOverall", output_folder);

%plot_phase(i_dir, false, Length_bounded, [6, floor(Lmin):1:Lmax], "Spindle length (um)", "length");
%plot_phase(i_dir, false, Flux_overal l_bounded, 0.0:0.1:2.0, "Overall flux (um/min)", "fluxOverall");

% Cross section cuts
ix_constL = find(abs(Length - L0) < 0.05); 
ix_constF = find(abs(Flux_k - F0) < 0.004); 
plot_section(i_dir, Lmax, Lmin, Length(:, 96), Flux_k(:, 96), Flux_b(:, 96), "plusOnly", output_folder);
plot_section(i_dir, Lmax, Lmin, Length(96, :), Flux_k(96, :), Flux_b(96, :), "minusOnly", output_folder);
plot_section(i_dir, Lmax, Lmin, Length(ix_constL), Flux_k(ix_constL), Flux_b(ix_constL), "constL", output_folder);
plot_section(i_dir, Lmax, Lmin, Length(ix_constF), Flux_k(ix_constF), Flux_b(ix_constF), "constF", output_folder);
plot_section(i_dir, Lmax, Lmin, diag(Length), diag(Flux_k), diag(Flux_b), "diag", output_folder);
plot_section(i_dir, Lmax, Lmin, diag(flip(Length)), diag(flip(Flux_k)), diag(flip(Flux_b)), "diagOppo", output_folder);
plot_section(i_dir, Lmax, Lmin, diag(Length_cut), diag(Flux_k_cut), diag(Flux_b_cut), "diag_zoom", output_folder);
plot_section(i_dir, Lmax, Lmin, diag(flip(Length_cut)), diag(flip(Flux_k_cut)), diag(flip(Flux_b_cut)), "diagOppo_zoom", output_folder);
ix = [66, 76, 86, 96, 106, 116, 126];
%ix = [66, 81, 96, 111, 126];
plot_section_discrete(i_dir, Lmax, Lmin, Length, Flux_k, Flux_b, ix, ix, "combine_slope1_Same_7pts", output_folder);
plot_section_discrete(i_dir, Lmax, Lmin, Length, Flux_k, Flux_b, flip(ix), ix, "combine_slope1_Oppo_7pts", output_folder);

function factor = get_correlation(deltaL, deltaF, epsL, epsF)
if deltaL > epsL && deltaF > epsF
    factor = 1;
elseif deltaL < -epsL && deltaF < -epsF
    factor = 2/3;
elseif deltaL < -epsL && deltaF > epsF
    factor = -2/3;
elseif deltaL > epsL && deltaF < -epsF
    factor = -1;
end
end

function plot_phase(i_dir, color, array, contours, plotname, filename, output_folder)
%array = smoothdata2(array, SmoothingFactor=0.001);
fontSize = 17;
v0p = 0.678;
v0m = 1.16;
switch i_dir
    case 1
        xlab = "Minus-end l-dep. strength (%)";
        ylab = "Plus-end l-dep. strength (%)";
        xticklab = compose("%i", [0 50 100 150 200]);
        yticklab = compose("%i", [0 50 100 150 200]);
        xlimits = [-4 162];
        ylimits = [-4 162];
        label = "l0";
    case 2
        %xlab = "Minus-end protein concentration (%)";
        %ylab = "Plus-end protein concentration (%)";
        xlab = "Change in minus-end protein number, \DeltaA^-( %)";
        ylab = "Change in plus-end protein number, -\DeltaA^+ (%)";
        xticklab = compose("%i", [0 50 100 150 200]-100);
        yticklab = compose("%i", [0 50 100 150 200]-100);
        %xlimits = [-4 175];
        %ylimits = [-4 170];
        xlimits = [-4 196];
        ylimits = [-4 196];
        %ylimits = [23 196];
        label = "A";
    case 3
        xlab = "Minus-end depolymerization rate, v^-_0 (um/min)";
        ylab = "Plus-end polymerization rate, v^+_0 (um/min)";
        xticklab = compose("%.1f",v0m * [0 0.5 1 1.5 2]);
        yticklab = compose("%.1f",v0p * [0 0.5 1 1.5 2]); 
        xlimits = [-4 150];
        ylimits = [-4 196];
        label = "v0";
end

fig = figure('Position',[0 50 720 540]);
hold on
switch i_dir
    case 1
        [M,c] = contourf(flip(flip(array), 2), contours, "ShowText","off", "FaceAlpha", 0.75, 'LineWidth', 1);
    case 2
        [M,c] = contourf(flip(array), contours, "ShowText","off", "FaceAlpha", 0.75, 'LineWidth', 1);
    case 3
        [M,c] = contourf(array, contours, "ShowText", "off", "FaceAlpha", 0.75, 'LineWidth', 1);
end
if ~color
    c.FaceColor = "none";
    c.LineWidth = 2;
else
    clabel(M, c, "FontSize", fontSize, "FontName", "Arial", "FontWeight", "normal");
end
title(plotname, "FontName", 'Arial');
xticks(([0 50 100 150 200]-4)) % subtract 4 since scan ranges from 5% to 195%, so x=1 is 5%
yticks(([0 50 100 150 200]-4)) % subtract 4 since scan ranges from 5% to 195%, so y=1 is 5%
xlabel(xlab);
ylabel(ylab);
xticklabels(xticklab);
yticklabels(yticklab);
xlim(xlimits);
ylim(ylimits);
box on
pbaspect([1 1 1]);
ax = gca;
ax.FontSize = fontSize;
ax.FontName = "Arial";
exportgraphics(fig, sprintf("%s/phase_%s_%s.pdf", output_folder,label, filename));
end

function plot_section(i_dir, Lmax, Lmin, L, F_k, F_b, filename, output_folder)
v0p = 0.678;
v0m = 1.16;
fontSize = 18;
invalid_color = [0.75 0.75 0.75];
switch i_dir
    case 1
        label = "l0";
        xticklab = compose("%i", [0 50 100 150 200]);
    case 2
        label = "A";
        xticklab = compose("%i", [0 50 100 150 200]-100);
    case 3
        label = "v0";
        if  (filename == "minusOnly")
            xticklab = compose("%.1f",v0m * [0 0.5 1 1.5 2]);
        else
            xticklab = compose("%.1f",v0p * [0 0.5 1 1.5 2]); 
        end
end
i_start = 1;
i_end = length(L);
for i_var = 1 : length(L)     
    if L(i_var) < Lmax && L(i_var) > Lmin
        i_start = i_var;
        break;
    end
end
for i_var = length(L) : -1 : 1
    if L(i_var) < Lmax && L(i_var) > Lmin
        i_end = i_var;
        break;
    end
end
fig = figure('Position', [50 50 400 400]);
yyaxis left
plot(i_start:i_end, smooth(L(i_start:i_end)), 'LineWidth', 3)
hold on
plot(1:i_start, L(1:i_start), 'LineStyle', '-', 'LineWidth', 3, 'Color', invalid_color)
plot(i_end:length(L), L(i_end:length(L)), 'LineStyle', '-',  'LineWidth', 3, 'Color', invalid_color)
ylim([0 Lmax + 2])
yticks([0, 10, 20, 30])
ylabel("Spindle length (um)", "FontSize", fontSize, "FontName", 'Arial')
yyaxis right
plot(i_start:i_end, smooth(F_k(i_start:i_end)), 'LineStyle', '-','LineWidth', 3, 'Color', [128 0 128]/255)
hold on
plot(i_start:i_end, smooth(F_b(i_start:i_end)), 'LineStyle', '-', 'LineWidth', 3, 'Color', [218 112 214]/255)
plot(1:i_start, F_k(1:i_start), 'LineStyle', '-', 'LineWidth', 3, 'Color', invalid_color)
plot(i_end:length(L), F_k(i_end:length(L)), 'LineStyle', '-',  'LineWidth', 3, 'Color', invalid_color)
plot(1:i_start, F_b(1:i_start), 'LineStyle', '-', 'LineWidth', 3, 'Color', invalid_color, 'Marker', 'none')
plot(i_end:length(L), F_b(i_end:length(L)), 'LineStyle', '-',  'LineWidth', 3, 'Color', invalid_color, 'Marker', 'none')
ylim([0.0 1.75])
yticks([0.0 0.5, 1.0, 1.5]);
ylabel("Flux velocity (um/min)", "FontSize", fontSize, "FontName", 'Arial');
%xlim([-10 length(L) + 10]);
xlim([-4 length(L) + 5]);
xticks([0, 50, 100, 150, 200]-4) % subtract 4 since scan ranges from 5% to 195%, so x=1 is 5%
%xlim([-4 196]);
%xlim([46 146]);
%xticks([0, length(L)])
%xticklabels([])
xticklabels(xticklab);
ax = gca;
ax.FontSize = fontSize;
ax.FontName = "Arial";
ax.YAxis(2).Color = [128 0 128]/255;
if i_dir == 1 %&& (filename == "plusOnly" || filename == "minusOnly")
    set(ax, 'xdir', 'reverse');
end
if i_dir == 2 && filename == "plusOnly"
    set(ax, 'xdir', 'reverse');
end

pbaspect([1 1 1]);
box off
exportgraphics(fig, sprintf("%s/plot_%s_LandF_%s.pdf", output_folder, label, filename))
end

function plot_section_discrete(i_dir, Lmax, Lmin, L, F_k, F_b, ix, iy, filename, output_folder)
fontSize = 18;
invalid_color = [0.75 0.75 0.75];
switch i_dir
    case 1
        label = "l0";
        xticklab = compose("%i", [0 50 100 150 200]);
    case 2
        label = "A";
        xticklab = compose("%i", [0 50 100 150 200]-100);
    case 3
        label = "v0";
        xticklab = compose("%.1f",v0m * [0 0.5 1 1.5 2]);
end

fig = figure('Position', [50 50 400 400]);
yyaxis left
plot(ix(1), L(iy(1), ix(1)), 'o', 'LineWidth', 3)
hold on
ylim([0 Lmax + 2])
yticks([0, 10, 20, 30])
ylabel("Spindle length (um)", "FontSize", fontSize, "FontName", 'Arial')
yyaxis right
plot(ix(1), F_k(iy(1), ix(1)), 'o', 'LineStyle', '-','LineWidth', 3, 'Color', [128 0 128]/255)
hold on
plot(ix(1), F_b(iy(1), ix(1)), 'o', 'LineStyle', '-', 'LineWidth', 3, 'Color', [218 112 214]/255)

for i = 2 : length(iy)
    disp(L(iy(i), ix(i)))
    if L(iy(i), ix(i)) > Lmax || L(iy(i), ix(i)) < Lmin
        yyaxis left
        plot(ix(i), L(iy(i), ix(i)), 'o', 'LineWidth', 3, 'Color', invalid_color)
        yyaxis right
        plot(ix(i), F_k(iy(i), ix(i)), 'o', 'LineStyle', '-','LineWidth', 3, 'Color', invalid_color)
        plot(ix(i), F_b(iy(i), ix(i)), 'o', 'LineStyle', '-', 'LineWidth', 3, 'Color', invalid_color)
    else
        yyaxis left
        plot(ix(i), L(iy(i), ix(i)), 'o', 'LineWidth', 3)
        yyaxis right
        plot(ix(i), F_k(iy(i), ix(i)), 'o', 'LineStyle', '-','LineWidth', 3, 'Color', [128 0 128]/255)
        plot(ix(i), F_b(iy(i), ix(i)), 'o', 'LineStyle', '-', 'LineWidth', 3, 'Color', [218 112 214]/255)
    end
end

ylim([0.0 1.75])
yticks([0.0 0.5, 1.0, 1.5]);
ylabel("Flux velocity (um/min)", "FontSize", fontSize, "FontName", 'Arial');
%xlim([-10 length(L) + 10]);
%xlim([-1 length(L) + 1]);
xticks([0, 50, 100, 150, 200]-4) % subtract 4 since scan ranges from 5% to 195%, so x=1 is 5%
xlim([46 146]);
%xticks([0, length(L)])
%xticklabels([])
xticklabels(xticklab);
ax = gca;
ax.FontSize = fontSize;
ax.FontName = "Arial";
ax.YAxis(2).Color = [128 0 128]/255;
if i_dir == 1 %&& (filename == "plusOnly" || filename == "minusOnly")
    set(ax, 'xdir', 'reverse');
end
if i_dir == 2 && filename == "plusOnly"
    set(ax, 'xdir', 'reverse');
end
pbaspect([1 1 1]);
box off
exportgraphics(fig, sprintf("%s/plot_%s_LandF_%s.pdf", output_folder, label, filename))
end
