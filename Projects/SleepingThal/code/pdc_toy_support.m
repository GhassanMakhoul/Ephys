%% Support script for PDC toy example .mlx
addpath('shared_toolboxes/fieldtrip-2019-0819');
addpath("/mnt/ernie_main/000_Data/SEEG/SEEG_Sleep_Staging/data/connectivity/5sDur_1sStride_v2/");
addpath('/home/ghassan/Documents/Research/Ephys/Projects/SleepingThal/code/shared_toolboxes/')
% load("/mnt/ernie_main/000_Data/SEEG/SEEG_Sleep_Staging/data/connectivity/5sDur_1sStride_v2/Spat52/Spat52_N2_imcoh.mat");
% load("/mnt/ernie_main/Ghassan/ephys/data/connectivity/Spat52/Spat52_11A_FIAS_PDC.mat")
subj = 'Spat51';
sleep_conn_dir = "/mnt/ernie_main/000_Data/SEEG/SEEG_Sleep_Staging/data/connectivity/5sDur_1sStride_v2/";
wake_interictal_dir = "/mnt/ernie_main/Ghassan/ephys/data/connectivity/";
%load sleep PDC
load(fullfile(sleep_conn_dir, subj, sprintf("%s_N3_imcoh.mat", subj))); %albeit a bit brittle to assume N3 is there
%load waking interictal PDC
load(fullfile(wake_interictal_dir, subj, sprintf("%s_5_FAS_PDC.mat", subj))) %in truth need to reslect filename each time

REGION_LABELS = ["SOZ", "NIZ", "PZ",];
BAND = 3 % for now investigating alpha
filename = "~/Documents/Research/Ephys/Data/all_pats_bipole_soz_labels.csv";
data = readtable(filename,'ReadVariableNames', false);
% Specify the Subject ID to filter

% Filter rows based on the Subject ID
headers = {'patID','bip','label'};
data.Properties.VariableNames = headers;

filteredData = data(strcmp(data.patID, subj), :);
inputBipoleLabels = cellstr(sleep_struct.bip_labels_used); %converts to cell array of str for direct comparison
regions = strings(size(inputBipoleLabels)); % Preallocate a string array for the results
% Iterate through each bipole_label in the input array
cleanedLabels = cellfun(@(x) regexprep(x, '[\s-]', ''), filteredData.bip, 'UniformOutput', false); %match bipole labels csv
for i = 1:length(inputBipoleLabels)
    % Find the row where the bipole_label matches
    idx = strcmp(cleanedLabels, inputBipoleLabels{i});
    % Get the corresponding region value
    if any(idx) % If a match is found
        bip_label =  filteredData.label(idx);
        regions(i) = map_label(bip_label);
    else
        regions(i) = "Not Found"; % Assign "Not Found" if no match exists
    end
end

%sleep dynamics
% pdc_dynamics = sleep_struct.PDC_all_windowed; state = 'Sleep N3'; %BAND x Time x Node x Node 
%
%waking interictal dynamics
pdc_dynamics = pdc.seizure.PDC_all_windowed;  state = 'Wake Rest FAS'; %BAND x Time x Node x NodeT= size(pdc_dynamics,2);

T =size(pdc_dynamics,2);
inwardConnectivity = struct();
outwardConnectivity = struct();
netConnectivity = struct();
for i = 1:length(unique(REGION_LABELS))
    reg = REGION_LABELS(i);
    reg_inds = strcmp(reg,regions);
    inwardConnectivity.(reg) = squeeze(mean(pdc_dynamics(BAND, :,:,reg_inds),3, 'omitnan')); % T x N_reg, sum of rows, inward connectivity
    outwardConnectivity.(reg) = squeeze(mean(pdc_dynamics(BAND, :,reg_inds,:),4,'omitnan'));
    netConnectivity.(reg) =inwardConnectivity.(reg) - outwardConnectivity.(reg);

end

% Plot inward connectivity
figure;
timePoints = linspace(0, T, T); % Example time range
subplot(3, 1, 1);
hold on;
for i = 1:length(unique(regions))
    reg = REGION_LABELS(i);
    if ismember(1, groupcounts(regions)) && strcmp("SOZ",reg)
        plot(timePoints,inwardConnectivity.(reg)')
    else
        shadedErrorBar(timePoints, inwardConnectivity.(reg)', {@mean, @std}, 'lineprops', '-');
    end
end
ylim(gca,[-.05,0.2]);
title('Peri-Ictal Inward Connectivity');
xlabel('Time (s)');
ylabel('PDC Inward Strength');
legend(REGION_LABELS);
hold off;


% Plot outward connectivity
subplot(3, 1, 2);
hold on;
for i = 1:length(REGION_LABELS)
    reg = REGION_LABELS(i);
    if ismember(1, groupcounts(regions)) && strcmp("SOZ",reg)
        plot(timePoints, outwardConnectivity.(reg)')
    else
        shadedErrorBar(timePoints, outwardConnectivity.(reg)', {@mean, @std}, 'lineprops', '-');
    end
end
ylim(gca, [-.05,0.1])
title('Peri-Ictal Outward Connectivity');
xlabel('Time (s)');
ylabel('PDC Outward Strength');
legend(REGION_LABELS);
hold off;


% Plot net connectivity
subplot(3, 1, 3);
hold on;
for i = 1:length(REGION_LABELS)
    reg = REGION_LABELS(i); 
    if ismember(1, groupcounts(regions)) && strcmp("SOZ",reg)
        plot(timePoints, netConnectivity.(reg)')
    else
        shadedErrorBar(timePoints, netConnectivity.(reg)', {@mean, @std}, 'lineprops', '-');
    end
end
ylim(gca, [-.12,.12])
title('Peri-Ictal Net Connectivity');
xlabel('Time (s)');
ylabel('PDC NET Strength');
legend(REGION_LABELS);
hold off;
sgtitle(sprintf('%s Dynamics', state));
% sgtitle('Waking Interictal Dynamics');

%%

saveas(gcf, sprintf('../viz/%s_dynamics_%s.png', subj,state))
% saveas(gcf, 'spat52_dynamics_wake.png')
function output = map_label(label)
    % Convert label to integer
    label = floor(label); % Ensures the label is an integer

    % Match cases using a switch statement
    switch label
        case 0
            output = "NIZ";
        case 1
            output = "SOZ";
        case 2
            output = "PZ";
        case 3
            output = "NIZ";
        otherwise
            error('Invalid label: %d', label); % Handle unexpected cases
    end
end

function output = plot_time_var_connectivity()

end