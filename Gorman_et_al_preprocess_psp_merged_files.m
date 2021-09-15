clear;

%This script demonstrates the pre-processing of the in-vivo recorded
%membrane potentials used in Gorman et al.

%%

%Data for the neuron in Figure 3 A, B

fnames = {'640.02.03.itd.mat','640.02.04.itd.mat','640.02.05.itd.mat','640.02.07.itd.mat'};
pathname = 'FrequencyConvergenceData/Mar3-1999/';


%other neurons

%{
fnames = {'631.04.03.itd.mat','631.04.04.itd.mat','631.04.06.itd.mat','631.04.07.itd.mat'};
pathname = 'FrequencyConvergenceData/Feb19-1999/';
%}

%{
fnames = {'632.03.03.itd.mat','632.03.04.itd.mat','632.03.05.itd.mat'};
pathname = 'FrequencyConvergenceData/Feb21-1999/';
%}

%{
fnames = {'629.02.04.itd.mat','629.02.05.itd.mat','629.02.06.itd.mat'};
pathname = 'FrequencyConvergenceData/Feb24-1999/';
%}

%{
fnames = {'629.07.01.itd.mat','629.07.02.itd.mat','629.07.05.itd.mat','629.07.06.itd.mat'};
pathname = 'FrequencyConvergenceData/Feb24-1999/';
%}

%{
fnames = {'642.07.03.itd.mat','642.07.04.itd.mat','642.07.05.itd.mat'};
pathname = 'FrequencyConvergenceData/Feb26-1999/';
%}

%{
fnames = {'639.01.05.itd.mat','639.01.06.itd.mat','639.01.07.itd.mat'};
pathname = 'FrequencyConvergenceData/Feb28-1999/';
%}

%{
fnames = {'624.07.04.itd.mat','624.07.05.itd.mat','624.07.06.itd.mat'};
pathname = 'FrequencyConvergenceData/Mar1-1999/';
%}

%{
fnames = {'640.01.04.itd.mat', '640.01.05.itd.mat','640.01.06.itd.mat','640.01.07.itd.mat'};
pathname = 'FrequencyConvergenceData/Mar3-1999/';
%}

%{
fnames = {'640.04.08.itd.mat','640.04.09.itd.mat','640.04.10.itd.mat'};
pathname = 'FrequencyConvergenceData/Mar3-1999/';
%}

%{
fnames = {'645.03.03.itd.mat','645.03.04.itd.mat','645.03.05.itd.mat','645.03.06.itd.mat'};
pathname = 'FrequencyConvergenceData/Mar4-1999/';
%}

%{
fnames = {'628.05.06.itd.mat','628.05.07.itd.mat','628.05.08.itd.mat','628.05.09.itd.mat'};
pathname = 'FrequencyConvergenceData/Mar5-1999/';
%}

%% Initialize
Nf = length(fnames);

y_m = cell(Nf,1);
x_m = cell(Nf,1);

%% Loop over component files

for jj = 1:Nf

%% Get the filename

    filename = fnames{jj};
    
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
    
%% Get traces   
    if isfield(d, 'trace')
        trace_x = 1000 * ((1:size(d.trace, 2))-1)/(d.params.adfc/double(d.params.ana_decimate));
        trace_y = ana_tomv * double(d.trace);

        stim_y = floor(ana_tomv * min(min(d.trace)) / 10 - 1) * 10;
                
        during_stim = (trace_x > d.params.delay + 10) & (trace_x < d.params.delay + 60);
        
        during_potentials = median(trace_y(:,during_stim), 2);
        
    end       
    
    x_values = unique(table_params.(2));
    x_values = x_values(~isnan(x_values));

    potentials = cell(size(x_values));  
    
    for k = 1:numel(x_values)      
        
        potentials{k} = during_potentials(table_params.(2) == x_values(k));        
        
    end  
    
    y_means = cellfun(@mean, potentials);
    
    
    x_m{jj} = x_values;
    y_m{jj} = y_means;
        
       
end    


%% Get V_stack and V_add

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