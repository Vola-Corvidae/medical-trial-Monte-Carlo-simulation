%% weighted_successrate_3d_plots_conditional.m
% Octave-Sscipt: marginalizing on 2D-grids

clear; clc; close all;

% === 1. read Data ===
filename = 'simulation_results.csv';
data = dlmread(filename, ',', 1, 0);  % omit header

MedianBat       = data(:,1);
MedianGps       = data(:,2);
ShapeBat        = data(:,3);
ShapeGps        = data(:,4);
SuccessRate     = data(:,5); % PoS
ValidRate       = data(:,6);
PosteriorWeight = data(:,7); % pMosPost

valid = PosteriorWeight > 0 & ~isnan(PosteriorWeight) & ~isnan(SuccessRate) & ~isnan(ValidRate);
w   = PosteriorWeight(valid);
sr  = SuccessRate(valid);
vr  = ValidRate(valid);

MedianBat = MedianBat(valid);
MedianGps = MedianGps(valid);
ShapeBat  = ShapeBat(valid);
ShapeGps  = ShapeGps(valid);

w = w / sum(w);


%% === Select a pair to analyze a slice ===
ShapeBat_target = 0.7;
ShapeGps_target = 0.7;


%% === 2. First projection: 2D grid MedianBat × MedianGps (marginalization via shapes) ===

mb_vals = sort(unique(MedianBat));
mg_vals = sort(unique(MedianGps));

Success_2D_Median = zeros(length(mb_vals), length(mg_vals));
Valid_2D_Median   = zeros(length(mb_vals), length(mg_vals));
Posterior_2D_Median = zeros(length(mb_vals), length(mg_vals));

for i = 1:length(mb_vals)
    for j = 1:length(mg_vals)
        mask = (MedianBat == mb_vals(i)) & (MedianGps == mg_vals(j));
        if any(mask)
            % Posterior expected PoS: E[PoS | Median]
            Success_2D_Median(i,j) = sum( sr(mask) .* w(mask) ) / sum(w(mask));
            % Posterior expected ValidRate: E[ValidRate | Median]
            Valid_2D_Median(i,j) = sum( vr(mask) .* w(mask) ) / sum(w(mask));
            % Marginalized Posterior Weight: P(Median | Valid)
            Posterior_2D_Median(i,j) = sum( w(mask) );
        else
            Success_2D_Median(i,j) = NaN;
            Valid_2D_Median(i,j) = NaN;
            Posterior_2D_Median(i,j) = NaN;
        end
    end
end

[MG_grid, MB_grid] = meshgrid(mg_vals, mb_vals);


%% === 3. Second projection: 2D grid ShapeBat × ShapeGps (marginalization via medians) ===

sb_vals = sort(unique(ShapeBat));
sg_vals = sort(unique(ShapeGps));

Success_2D_Shape = zeros(length(sb_vals), length(sg_vals));
Valid_2D_Shape   = zeros(length(sb_vals), length(sg_vals));
Posterior_2D_Shape = zeros(length(sb_vals), length(sg_vals));

for i = 1:length(sb_vals)
    for j = 1:length(sg_vals)
        mask = (ShapeBat == sb_vals(i)) & (ShapeGps == sg_vals(j));
        if any(mask)
            % Posterior expected PoS: E[PoS | Shape]
            Success_2D_Shape(i,j) = sum( sr(mask) .* w(mask) ) / sum(w(mask));
            % Posterior-expected ValidRate: E[ValidRate | Shape]
            Valid_2D_Shape(i,j) = sum( vr(mask) .* w(mask) ) / sum(w(mask));
            % Marginalized Posterior Weight: P(Shape | Valid)
            Posterior_2D_Shape(i,j) = sum( w(mask) );
        else
            Success_2D_Shape(i,j) = NaN;
            Valid_2D_Shape(i,j) = NaN;
            Posterior_2D_Shape(i,j) = NaN;
        end
    end
end

[SG_grid, SB_grid] = meshgrid(sg_vals, sb_vals);


%% === 4. Third projection: 2D grid MedianBat × MedianGps (conditioned on target shape pair) ===

% Find the lines that exactly match the target shape pair.
conditional_mask = (ShapeBat == ShapeBat_target) & (ShapeGps == ShapeGps_target);

if ~any(conditional_mask)
    disp(['FEHLER: Das Ziel-Shape-Paar (', num2str(ShapeBat_target), ', ', num2str(ShapeGps_target), ') wurde nicht in der CSV-Datei gefunden.']);
    Success_2D_Conditional = nan(length(mb_vals), length(mg_vals));
    Valid_2D_Conditional = nan(length(mb_vals), length(mg_vals));
    Posterior_2D_Conditional = nan(length(mb_vals), length(mg_vals));
    % The script does not terminate here, but instead displays empty plots.
else

    % Normalize the posterior weights only for this shape pair.
    w_conditional = w(conditional_mask);
    w_conditional = w_conditional / sum(w_conditional);

    MedianBat_cond = MedianBat(conditional_mask);
    MedianGps_cond = MedianGps(conditional_mask);
    sr_cond = sr(conditional_mask);
    vr_cond = vr(conditional_mask);

    Success_2D_Conditional = zeros(length(mb_vals), length(mg_vals));
    Valid_2D_Conditional   = zeros(length(mb_vals), length(mg_vals));
    Posterior_2D_Conditional = zeros(length(mb_vals), length(mg_vals));

    for i = 1:length(mb_vals)
        for j = 1:length(mg_vals)
            % Find the index of the current 4D scenario in the conditional data.
            idx = find((MedianBat_cond == mb_vals(i)) & (MedianGps_cond == mg_vals(j)), 1);

            if ~isempty(idx)
                % If the scenario exists, take its values.
                Success_2D_Conditional(i,j) = sr_cond(idx);
                Valid_2D_Conditional(i,j) = vr_cond(idx);
                Posterior_2D_Conditional(i,j) = w_conditional(idx);
            else
                Success_2D_Conditional(i,j) = NaN;
                Valid_2D_Conditional(i,j) = NaN;
                Posterior_2D_Conditional(i,j) = NaN;
            end
        end
    end
end

%% === 5. Visualize results (3D plots) ===
colormap(viridis);

% --- Figure 1: Projection via MedianBat × MedianGps ---
figure('Position', [100, 100, 1800, 500]);
%sgtitle('Median parameter projection (marginalized over shape)');

subplot(1,3,1);
surf(MG_grid, MB_grid, Success_2D_Median);
colorbar; title('Posterior expected  PoS');
xlabel('MedianGps'); ylabel('MedianBat'); zlabel('Success rate (PoS)'); axis tight;

subplot(1,3,2);
surf(MG_grid, MB_grid, Valid_2D_Median);
colorbar; title('Posterior expected ValidRate');
xlabel('MedianGps'); ylabel('MedianBat'); zlabel('ValidRate'); axis tight;

subplot(1,3,3);
surf(MG_grid, MB_grid, Posterior_2D_Median);
colorbar; title('Marginalized Posterior Weight P(Median | Valid)');
xlabel('MedianGps'); ylabel('MedianBat'); zlabel('Posterior Weight'); axis tight;


% --- Figure 2: Projection via ShapeBat × ShapeGps ---
figure('Position', [100, 700, 1800, 500]);
%sgtitle('Shape parameter projection (marginalized over median)');

subplot(1,3,1);
surf(SG_grid, SB_grid, Success_2D_Shape);
colorbar; title('Posterior expected PoS');
xlabel('ShapeGps'); ylabel('ShapeBat'); zlabel('SuccessRate (PoS)'); axis tight;

subplot(1,3,2);
surf(SG_grid, SB_grid, Valid_2D_Shape);
colorbar; title('Posterior expected ValidRate');
xlabel('ShapeGps'); ylabel('ShapeBat'); zlabel('ValidRate'); axis tight;

subplot(1,3,3);
surf(SG_grid, SB_grid, Posterior_2D_Shape);
colorbar; title('Marginalized Posterior Weight P(Shape | Valid)');
xlabel('ShapeGps'); ylabel('ShapeBat'); zlabel('Posterior Weight'); axis tight;


% --- Figure 3: Conditioned to target shape pair (MedianBat x MedianGps) ---
figure('Position', [100, 1300, 1800, 500]);
%sgtitle(['Median parameter projection (Conditioned on ShapeBat=', num2str(ShapeBat_target), ', ShapeGps=', num2str(ShapeGps_target), ')']);

subplot(1,3,1);
surf(MG_grid, MB_grid, Success_2D_Conditional);
colorbar; title('PoS with fixed shape');
xlabel('MedianGps'); ylabel('MedianBat'); zlabel('SuccessRate (PoS)'); axis tight;

subplot(1,3,2);
surf(MG_grid, MB_grid, Valid_2D_Conditional);
colorbar; title('ValidRate with fixed shape');
xlabel('MedianGps'); ylabel('MedianBat'); zlabel('ValidRate'); axis tight;

subplot(1,3,3);
surf(MG_grid, MB_grid, Posterior_2D_Conditional);
colorbar; title('Conditioned Posterior Weight P(Median | Valid, Shape)');
xlabel('MedianGps'); ylabel('MedianBat'); zlabel('Posterior Weight (Konditioniert)'); axis tight;

