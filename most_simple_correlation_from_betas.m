% Extract beta values
betaTable = SubjStats.table;

% Step 1: Select non-NaN beta values for a specific condition and chromophore
condName = 'Sem';        % your task condition
chromophore = 'hbo';     % or 'hbr'

% Filter rows
validRows = strcmp(betaTable.cond, condName) & ...
            strcmp(betaTable.type, chromophore) & ...
            ~isnan(betaTable.beta);

% Create a new table with only valid rows
filteredTable = betaTable(validRows, :);

% Optional: ensure only one scan per subject or average across scans if needed

% Step 2: Create unique channel labels (source-detector pairs)
filteredTable.Channel = strcat("S", string(filteredTable.source), "-D", string(filteredTable.detector));

% Step 3: Pivot data — rows: subjects, columns: channels, values: beta
% Convert to wide format using unstack
betaWide = unstack(filteredTable, 'beta', 'Channel', 'GroupingVariables', 'subject');

% Step 4: Remove any rows with NaNs
betaMatrix = betaWide{:, 2:end};  % numeric part
validSubjects = all(~isnan(betaMatrix), 2);
betaMatrix = betaMatrix(validSubjects, :);

% Step 5: Compute correlation matrix across channels
connMatrix = corrcoef(betaMatrix);

% Step 6: Visualise
figure;
imagesc(connMatrix);
colorbar;
title(['Functional Connectivity (', condName, ', ', chromophore, ')']);
xlabel('Channels'); ylabel('Channels');
