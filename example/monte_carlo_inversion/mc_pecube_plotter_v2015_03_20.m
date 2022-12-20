function mc_pecube_plotter_v141127
% mc_pecube_plotter
% by byron adams, todd ehlers, willi kappler
% beginning 11.2013 (check version number and date @ line 70 of main
% script)
%
% call: matlab -nodisplay -nosplash -noFigureWindows -r "run('mc_pecube_plotter_v141127.m')"
%
%
% this script uses the output data from the pecube monte carlo simulations
% to calculate and plot mean exhumation rate histories for each sample. the
% output data of this script reports the standard arithimitic mean and 2
% standard errors on the mean from all acceptable exhumation rate
% histories.
%
% See the 00readme.rtf file for instructions on how to use this.
%
%-------------------------------------------------------------------------%
%
%
close all
clear all
% load variables from file if they do not exist in memory or as function argument


% max_erate: mc_maximum_erosion_rate
% time_steps: time step years
% ages_gfit(w,cnt2(u),u) = predicted ages
% chisq_gfit(w,cnt2(u),u) = chi squared
% erate_gfit(ct,cnt2(u),u) = monte carlo estimated erosion rates
% cnt2: number of good fits
% u: number of samles (datasets) per line = j
% w: number of chronometers = k
% ct: number of time steps = ntime

max_erate = 0.0;
time_steps = [];

file_id = fopen('Pecube.in');

line = fgets(file_id);

while ischar(line)
    if strncmp(line, 'mc_max_erosion_rate:', 20)
        max_erate = str2double(line(22:end));
        % fprintf('max_erate: "%f"\n', max_erate);
        % fprintf('line: "%s"\n', line(22:end));
    elseif strncmp(line, '$(a) (b) (c)', 12)
        line = fgets(file_id);
        while line(1:2) ~= '$'
            % fprintf('line: "%s"\n', line);
            values = strsplit(line);
            if (numel(values)) > 2
                time_steps(end + 1) = str2double(values(1));
            end
            line = fgets(file_id);
        end
    end
    line = fgets(file_id);
end

fclose(file_id);

%display(max_erate);
%display(time_steps);

file_id = fopen('erate_gfit.txt');

% ignore first line = header
line = fgets(file_id);
line = fgets(file_id);

erate_gfit = [];

while ischar(line)
    values = sscanf(line, '%d %d %d %f');
    erate_gfit(values(2), values(1), values(3)) = values(4);
    line = fgets(file_id);
end

fclose(file_id);

%display(erate_gfit(:,1,:));

file_id = fopen('good_fit_stats.txt');

% ignore first line = header
line = fgets(file_id);
line = fgets(file_id);

ages_gfit = [];
chisq_gfit = [];

while ischar(line)
    values = sscanf(line, '%d %d %d %f %f');
    ages_gfit(values(2), values(1), values(3)) = values(5);
    chisq_gfit(values(2), values(1), values(3)) = values(4);
    line = fgets(file_id);
end

fclose(file_id);

%display(ages_gfit(:,1,:));
%display(chisq_gfit(:,1,:));

file_id = fopen('chrono_data.txt');

% ignore first line = header
line = fgets(file_id);
line = fgets(file_id);

data_cube_names = {};
data_cube_ages = [];
data_cube_sigma = [];

while ischar(line)
    values = textscan(line, '%d %d %s %f %f');
    data_cube_names{values{1}, values{2}} = values{3};
    data_cube_ages(values{1}, values{2}) = values{4};
    data_cube_sigma(values{1}, values{2}) = values{5};
    line = fgets(file_id);
end

fclose(file_id);

display(data_cube_names);
display(data_cube_ages);
display(data_cube_sigma);

% exit(0);

% find the dimension of output files
    [~,~,num_sam] = size(erate_gfit);
    ct = 0;
%
% write the time step array for ploting
%    for y = 1:length(time_steps)
%        times(y + ct) = time_steps(y);
%        if y ~= length(time_steps)
%            times(y*2) = time_steps(y);
%        end
%        ct = ct + 1;
%    end
%    times = transpose(times);

times = transpose(time_steps(1 : end - 1));


%
% loop through all the samples to calculate statistics and plot summary
% exhumation history figures.
    for u = 1:num_sam
        % pull data for each sample
            data = erate_gfit(:,:,u);
            system = data_cube_names(:,u);
            time = data_cube_ages(:,u);
            time_std = data_cube_sigma(:,u);
        %
        % remove zero and NaN values
            real = find(~isnan(time));
            time = time(real);
            time_std = time_std(real);
            system = system(real);
            data(data == 0) = NaN;
            num_sim = sum(~isnan(data(1,:)));
        %
        % calculate the mean and standard error
            mean_rate = nanmean(data,2);
            se2_rate = nanstd(data,0,2)/sqrt(num_sim)*2;
            %
            se_high = mean_rate + se2_rate;
            se_low = mean_rate - se2_rate;
        %
        % concatonate data into a new matrix for storage
            final_data(:,1,u) = mean_rate; %#ok<*AGROW>
            final_data(:,2,u) = se2_rate;
            final_data(1,3,u) = num_sim;
        %
        % plot summary exhumation history
            figure(u)
            hold on

            % plot(times, mean_rate,'r-','LineWidth',2)
            % plot(times, se_high,'b-')
            % plot(times, se_low,'b-')

            stairs(times, mean_rate,'r-','LineWidth',2)
            stairs(times, se_high,'b-')
            stairs(times, se_low,'b-')

            axis([0 time_steps(1) 0 max_erate])
            xlabel('Time (Ma)')
            ylabel('Erosion Rate (km/m.y.)')
            title_name = ['Acceptable Erosion Histories,   ','N = ',num2str(num_sim)];
            title(title_name)
            legend('Mean Rate','2 SE')
        %
        % plot chronometer ages for comparison
            for w = 1:length(system)
                bar(time(w),max_erate,time_std(w),'c')
                text(time(w),0.5,system{w})
            end
            
        %
        % save figures as raster and vector files
            name1 = ['sample',num2str(u),'.jpg'];
            name2 = ['sample',num2str(u),'.eps'];
            h = figure(u);
            print(h,'-djpeg',name1)
            print(h,'-depsc',name2)
    end
%
% save statistical data
    save final_data.mat final_data
%
