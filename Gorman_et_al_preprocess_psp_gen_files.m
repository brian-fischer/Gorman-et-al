clear;

%This script demonstrates the pre-processing of the in-vivo recorded
%membrane potentials used in Gorman et al.


%% Filenames and paths for different neurons

%Data for the neuron in Figure 3 E,F

index = [1 3 4];
filename = '680.02.04.gen.mat';
pathname = 'FrequencyConvergenceData/Apr24-2000/';


%other neurons
%{
index = [1 3 5];
filename = '405.03.00.gen.mat';
pathname = 'FrequencyConvergenceData/Apr13-2000/';
%}

%{
index = [1 3 5];
filename = '405.04.00.gen.mat';
pathname = 'FrequencyConvergenceData/Apr13-2000/';
%}

%{
index = [1 2 5];
filename = '680.02.03.gen.mat';
pathname = 'FrequencyConvergenceData/Apr24-2000/';
%}

%{
index = [2 4 8];
filename = '663.05.03.gen.mat';
pathname = 'FrequencyConvergenceData/Jul8-1999/';
%}

%{
index = [2 4 5];
filename = '665.03.03.gen.mat';
pathname = 'FrequencyConvergenceData/Jul9-1999/';
%}

%{
index = [2 4 8];
filename = '667.02.03.gen.mat';
pathname = 'FrequencyConvergenceData/Jul14-1999/';
%}

    
%% Load the data

d = load(fullfile(pathname, filename));
    
%% Get parameters

if ~isfield(d.params, 'ana_tomv')
    ana_tomv = 1e-1 * 0.306523;
elseif isstruct(d.params.ana_tomv)
    ana_tomv = 1e-1 * d.params.ana_tomv.chan1;
else
    ana_tomv = 1e-1 * d.params.ana_tomv;
end
   
table_params = table();
for k = 1:size(d.trial_params_columns, 1)
    if iscell(d.trial_params)
    % Mixed values (numerical and text)
        try
            table_params(:,k) = num2cell(cell2mat(d.trial_params(:,k)));
        catch
            table_params(:,k) = cellfun(@num2str, d.trial_params(:,k), 'UniformOutput', false);
        end
    else
        table_params(:,k) = num2cell(d.trial_params(:,k));
    end
    table_params.Properties.VariableNames{k} = d.trial_params_columns(k,:);
end
clear k    
    
    
%% Find variables that have been varied

params_columns = cellstr(d.trial_params_columns);
tunable_columns = {'abi', 'itd', 'iid', 'bc'};
tunable_labels = struct('abi', 'ABI [dB SPL]', 'itd', 'ITD [µs]', 'iid', 'IID [dB]', 'bc', 'BC [%]');

% Tunable (numeric) columns with varied values:
var_columns = {};
% Columns with varied values that are not numeric:
grp_columns = {};

disp('Parameter | # unique | Values')
for c = 1:numel(params_columns)
    try
        unique_col_values = unique(table_params.(c));
        had_nan = false;
        if ischar(unique_col_values)
            unique_col_values = cellstr(unique_col_values);
        end
        if isnumeric(unique_col_values)
            had_nan = any(isnan(unique_col_values));
            unique_col_values = unique_col_values(~isnan(unique_col_values));
        end
        if strcmp('stim', params_columns{c}) && isnumeric(unique_col_values)
            tunable_columns{end+1} = 'stim';
            tunable_labels.stim = 'Stim/Freq [Hz]';
        end
        nc = numel(unique_col_values);
        fprintf('%9s | %8d | ', params_columns{c}, nc)
        if iscellstr(unique_col_values) %#ok<ISCLSTR>
            for k = 1:numel(unique_col_values)
                fprintf('''%s''', unique_col_values{k})
                if k < numel(unique_col_values)
                    fprintf(',')
                else
                    fprintf('\n')
                end
            end
        elseif isnumeric(unique_col_values) && ~isempty(unique_col_values)
            for k = 1:numel(unique_col_values)
                fprintf('%.0f', unique_col_values(k))
                if k < numel(unique_col_values)
                    fprintf(',')
                elseif had_nan
                    fprintf(',nan\n')
                else
                    fprintf('\n')
                end
            end
        else
            fprintf('\n')
        end
        if nc > 1
            if any(strcmp(params_columns{c}, tunable_columns))
                var_columns{end+1} = c;
            else
                grp_columns{end+1} = c;
            end
        end
    catch previous_error
        disp(['Error with col ' params_columns{c}])
        rethrow(previous_error)
    end
end

%% Group trials

grp_selectors = {true(size(table_params, 1), 1)};
grp_selectors_labels = {''};

if numel(grp_columns) > 0
    for g = 1:numel(grp_columns)
        new_grp_selectors = {};
        new_grp_selectors_labels = {};
                
        grp_unique_values = unique(table_params.(grp_columns{g}));
        for gg = 1:numel(grp_unique_values)
            grp_value = grp_unique_values{gg};
            grp_col_selectors = strcmp(table_params.(grp_columns{g}), grp_value);
            for ggg = 1:numel(grp_selectors)
                new_grp_selectors{end+1} = grp_selectors{ggg} & grp_col_selectors;
                new_grp_selectors_labels{end+1} = sprintf('%s | %s=%s', grp_selectors_labels{ggg}, grp_columns{g}, grp_value);
            end
        end
        grp_selectors = new_grp_selectors;
        grp_selectors_labels = new_grp_selectors_labels;
    end
    % Crop superfluous ' | ' at the beginning:
    grp_selectors_labels = cellfun(@(s) s(4:end), grp_selectors_labels, 'UniformOutput', 0);
end


%%

%Initialize
y_m = cell(3,1);
x_m = cell(3,1);

if isfield(d, 'trace')
    
   
trace_x = 1000 * ((1:size(d.trace, 2))-1)/(d.params.adfc/double(d.params.ana_decimate));
trace_y = ana_tomv * double(d.trace);
      
stim_y = floor(ana_tomv * min(min(d.trace)) / 10 - 1) * 10;    

during_stim = (trace_x > d.params.delay + 10) & (trace_x < d.params.delay + 60);

during_potentials = median(trace_y(:,during_stim), 2);

%
for g = 1:numel(grp_selectors)
    grp_select = grp_selectors{g};
    grp_label = grp_selectors_labels{g};
    for c = 1:numel(var_columns)
        x_values = unique(table_params.(var_columns{c}));
        x_values = x_values(~isnan(x_values));
        potentials = cell(size(x_values));
        spont_potentials = cell(size(x_values));
        
        for k = 1:numel(x_values)
            select = table_params.(var_columns{c}) == x_values(k);
            potentials{k} = during_potentials(select & grp_select);
        end
        y_means = cellfun(@mean, potentials);
                
         
        if g==index(1) 
            x_m{1} = x_values;
            y_m{1} = y_means;
        elseif g == index(2)
            x_m{2} = x_values;
            y_m{2} = y_means;
        elseif g == index(3)
            x_m{3} = x_values;
            y_m{3} = y_means;
        end
    end
        
end

end % if isfield(d, 'trace')
    

%% Get V_stack and V_add

Nf = length(index);

V_stack = y_m{Nf};

V_add = zeros(size(V_stack));
for n = 1:Nf - 1
    V_add = V_add + y_m{n};
end
min_response = min(cellfun(@min,y_m(1:Nf-1)));
V_add = V_add - (Nf-2)*min_response;
    

%% plot

figure(1);clf;
subplot(211)
plot(trace_x, trace_y, 'Color', [0,0,0,.05]);
set(gca,'box','off','fontsize',12)
xlabel('Time (ms)','fontsize',15)
ylabel('Membrane potential (mV)','fontsize',15)

subplot(223);hold on;
g = linspace(0.1,.5,Nf-1);
for n=1:Nf-1
    plot(x_m{n},y_m{n},'o-','Color',[1 1 1]*g(n));
end
plot(x_m{Nf},y_m{Nf},'ro-')
set(gca,'box','off','fontsize',12)
xlabel('ITD (\mus)','fontsize',15)
ylabel('Membrane potential (mV)','fontsize',15)

subplot(224)
plot(V_add,V_stack,'ko','MarkerFaceColor','k')
set(gca,'box','off','fontsize',12)
xlabel('V_{add} (mV)','fontsize',15)
ylabel('V_{stack} (mV)','fontsize',15)


