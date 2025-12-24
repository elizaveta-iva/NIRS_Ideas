%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%After noisy channels have been pruned - average HbOs and HbRs based on ROI.
%% ! At the end, each participant will have 4 "channels" - in reality, each of those "channels" is averaged signal per ROI


% Define ROI channel groups with the specific source-detector pair for averaging
roiPairs.leftparietal   = [6 6; 6 5; 7 5; 8 5; 8 6];     % Averaged signal will be assigned to source 6, detector 6
roiPairs.rightparietal  = [1 2; 1 1; 3 1; 2 1; 3 2];     % Averaged signal will be assigned to source 3, detector 1
roiPairs.rightfrontal   = [4 3; 4 4; 5 4; 5 3];           % Averaged signal will be assigned to source 5, detector 3
roiPairs.leftfrontal    = [9 8; 9 7; 10 7; 10 8];         % Averaged signal will be assigned to source 9, detector 8

Hb_ROI = Hb_pruned;  % Copy to keep metadata

for i = 1:length(Hb_pruned)
    link = Hb_pruned(i).probe.link;
    data = Hb_pruned(i).data;
    t = Hb_pruned(i).time;

    roiData = [];
    roiLink = table();

    % For each ROI
    roiNames = fieldnames(roiPairs);
    for r = 1:length(roiNames)
        roiName = roiNames{r};
        pairs = roiPairs.(roiName);

        % For each type (e.g., hbo/hbr)
        types = unique(link.type);
        for tIdx = 1:length(types)
            thisType = types{tIdx};

            % Find matching channels for the ROI and type
            matchIdx = false(height(link),1);
            for p = 1:size(pairs,1)
                matchIdx = matchIdx | ...
                    (link.source == pairs(p,1) & ...
                     link.detector == pairs(p,2));
            end
            typeIdx = strcmp(link.type, thisType);
            finalIdx = find(matchIdx & typeIdx);

            if ~isempty(finalIdx)
                if length(finalIdx) > 1
                    % Average the channels if more than one exists for this ROI
                    avgSeries = nanmean(data(:, finalIdx), 2);
                else
                    % Keep the data as is if only one channel
                    avgSeries = data(:, finalIdx);
                end
                roiData = [roiData avgSeries];  % Append to ROI data

                % Update probe.link with the new source-detector pair (averaged)
                if strcmp(roiName, 'leftparietal')
                    newRow = table(6, 6, {thisType}, 'VariableNames', {'source', 'detector', 'type'});
                elseif strcmp(roiName, 'rightparietal')
                    newRow = table(3, 1, {thisType}, 'VariableNames', {'source', 'detector', 'type'});
                elseif strcmp(roiName, 'rightfrontal')
                    newRow = table(5, 3, {thisType}, 'VariableNames', {'source', 'detector', 'type'});
                elseif strcmp(roiName, 'leftfrontal')
                    newRow = table(9, 8, {thisType}, 'VariableNames', {'source', 'detector', 'type'});
                end

                roiLink = [roiLink; newRow];
            end
        end
    end

    % Store averaged data and update probe.link with region names
    Hb_ROI(i).data = roiData;
    Hb_ROI(i).probe.link = roiLink;
    Hb_ROI(i).time = t;
end

job = nirs.modules.RenameStims();
job.listOfChanges = {'LargeNum', 'SmallNum';}
Hb_ROI_combo = job.run(Hb_ROI);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Connectivity 

job = nirs.modules.Connectivity;
job.corrfcn=@(data)nirs.sFC.ar_corr(data,'4x',true); % Whitened correlation (using Pmax 4 x FS)
job.divide_events = 1;  % if true will parse into multiple conditions
job.min_event_duration=5;  % minimum duration of events
job.ignore = 0;  % time at transitions (on/off) to ignore (only valid if dividing events)
ConnStats = job.run(Hb_ROI_combo);

for i = 1:length(ConnStats)
    T = ConnStats(i).table;

    % Filter rows to keep only same-type connections
    keep_idx = strcmp(T.TypeOrigin, T.TypeDest);
    T_filtered = T(keep_idx, :);

    % Create a new object copying ConnStats(i), but with filtered table
    ConnStats_filtered(i) = ConnStats(i); % copy object
    ConnStats_filtered(i).table = T_filtered; % try assign here if allowed
end



job = nirs.modules.MixedEffectsConnectivity();
job.formula = 'R~-1+cond'
%  MixedEffectsConnectivity with properties:
%         formula: 'R ~ -1 + cond'
%     dummyCoding: 'full'
%      centerVars: 1
%            name: ''
%         prevJob: []
GroupConnStats = job.run(ConnStats);


disp(GroupConnStats.conditions);

ContrastStats = GroupConnStats.ttest('LargeNum + SmallNum - Sem');
ContrastStats.table

ContrastStats2 = GroupConnStats.ttest('Sem - 0.5*(LargeNum + SmallNum)');
ContrastStats2.table

Contrast = GroupConnStats_combo.ttest('Sem - SmallNum');
a=Contrast.table
Contrast.draw('z','q<0.05')

job = nirs.modules.DiscardStims;
job.listOfStims = {'SmallNum'};
ConnStats1 = job.run(ConnStats1);

conds = GroupConnStats.conditions;

contrast_vec = zeros(1, length(conds));
contrast_vec(strcmp(conds, 'Sem'))       = 1;
contrast_vec(strcmp(conds, 'SmallNum'))  = -0.5;
contrast_vec(strcmp(conds, 'LargeNum'))  = -0.5;

Contrast = GroupConnStats.ttest(contrast_vec);
Contrast.table


T = table();

for subj = 1:length(ConnStats)
    conds = ConnStats(subj).conditions;
    for c = 1:length(conds)
        % Example: mean connectivity over channel pairs
        meanR = mean(mean(ConnStats(subj).R(:,:,c)));
        newRow = {subj, conds{c}, meanR};
        T = [T; cell2table(newRow, 'VariableNames', {'SubjectID', 'Condition', 'MeanR'})];
    end
end

writetable(T, 'ConnStats_long.csv');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plotting the correlation heatmap for ROI

t = GroupConnStats.table;
conditions = unique(t.condition);

for i = 1:length(conditions)
    conditionName = conditions{i};
    filtered = t(t.qvalue < 0.05 & strcmp(t.condition, conditionName), :);

    % Ensure columns are string
    filtered.SourceOrigin = string(filtered.SourceOrigin); 
    filtered.DetectorOrigin = string(filtered.DetectorOrigin);
    filtered.SourceDest = string(filtered.SourceDest); 
    filtered.DetectorDest = string(filtered.DetectorDest);

    for j = 1:height(filtered)
        if ismember(filtered.SourceOrigin(j), ["9", "10", "6", "8", "7"])
            filtered.SourceOrigin(j) = 'left';
        else
            filtered.SourceOrigin(j) = 'right';
        end
        
        if ismember(filtered.DetectorOrigin(j), ["2", "1", "5", "6"])
            filtered.DetectorOrigin(j) = 'parietal';
        else
            filtered.DetectorOrigin(j) = 'frontal';
        end

        if ismember(filtered.SourceDest(j), ["9", "10", "6", "8", "7"])
            filtered.SourceDest(j) = 'left';
        else
            filtered.SourceDest(j) = 'right';
        end
        
        if ismember(filtered.DetectorDest(j), ["2", "1", "5", "6"])
            filtered.DetectorDest(j) = 'parietal';
        else
            filtered.DetectorDest(j) = 'frontal';
        end
    end

    % Create channel labels
    filtered.OriginChannel = strcat(filtered.SourceOrigin, "-", filtered.DetectorOrigin, "_", filtered.TypeOrigin);
    filtered.DestChannel = strcat(filtered.SourceDest, "-", filtered.DetectorDest, "_", filtered.TypeDest);

    % Unique channels
    channels = unique([filtered.OriginChannel; filtered.DestChannel]);
    n = length(channels);
    connMatrix = NaN(n, n);

    for j = 1:height(filtered)
        src = find(strcmp(channels, filtered.OriginChannel(j)));
        dest = find(strcmp(channels, filtered.DestChannel(j)));
        connMatrix(src, dest) = filtered.R(j);
        connMatrix(dest, src) = filtered.R(j);
    end

    % Plot
    figure;
    h = heatmap(channels, channels, connMatrix, ...
        'Colormap', parula, ...
        'ColorLimits', [-1 1], ...
        'MissingDataColor', [0.8 0.8 0.8], ...
        'MissingDataLabel', 'No connection');
    title(['Functional Connectivity - ', conditionName, ' (q < 0.05)']);
    xlabel('ROI');
    ylabel('ROI');
end

%%%%%%%%%%%%%%%%%%% getting individual z values for contrast
cond = 'SmallNum';
srcOrig = 9;
detOrig = 8;
typeOrig = 'hbo';

srcDest = 5;
detDest = 3;
typeDest = 'hbo';

% Initialise array to store Z values for each recording
Z_values = nan(1, numel(ConnStats));  % Preallocate with NaNs

for subj = 1:numel(ConnStats)
    thisTable = ConnStats(subj).table;  % Access one recording's connectivity table

    for i = 1:height(thisTable)
        if strcmp(thisTable.condition{i}, cond) && ...
           thisTable.SourceOrigin(i) == srcOrig && ...
           thisTable.DetectorOrigin(i) == detOrig && ...
           strcmp(thisTable.TypeOrigin{i}, typeOrig) && ...
           thisTable.SourceDest(i) == srcDest && ...
           thisTable.DetectorDest(i) == detDest && ...
           strcmp(thisTable.TypeDest{i}, typeDest)

            % Store the Z value for this subject
            Z_values(subj) = thisTable.Z(i);
            break  % stop after first match in this recording
        end
    end
end

% Display all collected Z values
disp(Z_values')


