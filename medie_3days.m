clc


%% controllers walking percentage
clc
perc_walk_CO = walk_tot(1:37);


walk_co_mean = mean(perc_walk_CO);
walk_co_std = std(perc_walk_CO);


fprintf('total CONTROLLERS walk percentage mean: %.2f ± %.2f percent',walk_co_mean,walk_co_std )


%% fallers walking percentage

perc_walk_FL = [walk_tot(38:end) 4.44 5.93 3.95 11.45 3.31 2.146 3.12] ;



walk_fl_mean = mean(perc_walk_FL);
walk_fl_std = std(perc_walk_FL);

fprintf('\ntotal FALLERS walk percentage mean: %.2f ± %.2f percent\n', walk_fl_mean,walk_fl_std);

figure
subplot(121)
boxplot(perc_walk_CO);
title('Controllers % of walk ');

subplot(122)
boxplot(perc_walk_FL)
title('Fallers % of walk')


disp('-----------------------------------------------------------------------------------------');

%% CONTROLLERS number of steps

step_co = step_tot(1:37);

mean_step_co = floor(mean(step_co));
std_step_co = floor(std(step_co));


fprintf('total controllers step mean: %d ± %d',mean_step_co, std_step_co );

%% FALLERS number of steps

step_fl = [step_tot(38:end) 9368 32904 7759 24163 5234 5655 4748];

mean_step_fl = floor(mean(step_fl));
std_step_fl = floor(std(step_fl));

fprintf('\ntotal fallers step mean: %d ± %d\n',mean_step_fl, std_step_fl );


figure
subplot(121)
boxplot(step_co);
title('Controllers number of steps ');

subplot(122)
boxplot(step_fl)
title('Fallers number of steps ')

%% CONTROLLERS average stride duration

median_bout_duration_CO = median_bout_tot(1:37);

mean_median_bout_duration_CO = floor(mean(median_bout_duration_CO));
std_median_bout_duration_CO = floor(std(median_bout_duration_CO));

disp('-----------------------------------------------------------------------------------------');
fprintf('\nControllers average stride duration per walking bout %.2f ± %.2f',mean_median_bout_duration_CO, std_median_bout_duration_CO  );

%% FALLERS average stride duration

median_bout_duration_FL = [median_bout_tot(38:end) 83 101 78.5 92 90.5 80.5 79];

mean_median_bout_duration_FL = mean(median_bout_duration_FL);
std_median_bout_duration_FL = floor(std(median_bout_duration_FL));



fprintf('\nFallers average stride duration per walking bout %.2f ± %.2f',mean_median_bout_duration_FL, std_median_bout_duration_FL );

figure
subplot(121)
boxplot(median_bout_duration_CO);
title('Controllers average stride duration ');

subplot(122)
boxplot(median_bout_duration_FL)
title('Fallers average stride duration ')


