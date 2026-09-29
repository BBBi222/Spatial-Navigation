%% plot_ERPs_crowdedness_GLOBAL_allParticipants.m
% ERPs time-locked to fixation onset
% Compare "Less crowded" vs "Crowded" using GLOBAL thresholds
% with 2 definitions of crowdedness:
%   (1) person_nr   = number of people
%   (2) person_area = area covered by people
%
% Requires:
%  - per participant: 5b2_epochedFixationRejection_withTargets_vehicle_as_objects_*.set
%  - CSV: filtered_preprocessed_all_participants.csv with columns:
%       participant, event, valid, fixation id, person_nr, person_area

clear; clc; close all;

%% === SETTINGS ===========================================================
dataRoot   = '/MATLAB Drive/EEG_data';

% CSV with person_nr and person_area
csvFile = '/MATLAB Drive/fixations/filtered_preprocessed_all_participants.csv';

% Channel of interest
chanName   = 'Cz';   % 'C3', 'Cz', 'C4'

% Minimum number of trials per condition (per participant)
minTrialsPerCond = 20;

% GLOBAL extremes (bottom q vs top q) across ALL participants
qCrowd = 0.20;   % e.g. 0.20 = bottom 20% vs top 20%

% Separate legend figure?
MAKE_LEGEND_FIG = true;

% Style
STYLE.LINE_W     = 1.6;
STYLE.AX_FONTSZ  = 10;
STYLE.LAB_FONTSZ = 10;
STYLE.GRID_ON    = true;
STYLE.BOX_OFF    = true;
STYLE.X_LIM      = [-300 500];
STYLE.Y_LIM      = [-2 2.5];

% Path to EEGLAB
addpath('/MATLAB Drive/eeglab/eeglab2025.1.0');
eeglab nogui;

%% === OUTPUT FOLDERS =====================================================
runTag = datestr(now,'yyyymmdd_HHMMSS');
outDir = fullfile(dataRoot, 'ERP_plots_crowdedness_GLOBAL', runTag);
if ~exist(outDir, 'dir'), mkdir(outDir); end

%% === LOG ================================================================
logFile = fullfile(outDir, sprintf('log_ERPs_%s_CROWD_GLOBAL.txt', chanName));
if exist(logFile, 'file'), delete(logFile); end
diary(logFile); diary on;

fprintf('=== SETTINGS ===\n');
fprintf('dataRoot: %s\n', dataRoot);
fprintf('csvFile : %s\n', csvFile);
fprintf('chanName: %s\n', chanName);
fprintf('minTrialsPerCond: %d\n', minTrialsPerCond);
fprintf('qCrowd (GLOBAL): %.2f\n', qCrowd);
fprintf('MAKE_LEGEND_FIG: %d\n', MAKE_LEGEND_FIG);
fprintf('Output folder: %s\n\n', outDir);

%% === LOAD CSV ===========================================================
opts = detectImportOptions(csvFile, 'VariableNamingRule','preserve');
opts.SelectedVariableNames = {'participant','event','valid','fixation id','person_nr','person_area'};
T = readtable(csvFile, opts);

% Robust filter: event == fixation AND valid == True
isFix = strcmpi(strtrim(string(T.('event'))), "fixation");
v = lower(strtrim(string(T.('valid'))));
isValid = (v == "true");

T = T(isFix & isValid, :);

fprintf('CSV rows after filter (event==fixation & valid==true): %d\n', height(T));

T.person_nr   = double(T.person_nr);
T.person_area = double(T.person_area);

%% === GLOBAL THRESHOLDS ==================================================
mCountAll = T.person_nr;
mCountAll = mCountAll(~isnan(mCountAll));

mAreaAll  = T.person_area;
mAreaAll  = mAreaAll(~isnan(mAreaAll));

if isempty(mCountAll) || isempty(mAreaAll)
    error('After filtering, person_nr or person_area is empty (all NaN?).');
end

loCount = prctile(mCountAll, qCrowd*100);
hiCount = prctile(mCountAll, (1-qCrowd)*100);

loArea  = prctile(mAreaAll,  qCrowd*100);
hiArea  = prctile(mAreaAll,  (1-qCrowd)*100);

fprintf('\n=== GLOBAL THRESHOLDS (q=%.2f) ===\n', qCrowd);
fprintf('COUNT: Less <= %.6g | Crowded >= %.6g\n', loCount, hiCount);
fprintf('AREA : Less <= %.6g | Crowded >= %.6g\n', loArea,  hiArea);

%% === FIND PARTICIPANTS ==================================================
tmp = dir(fullfile(dataRoot));
participants = [];
inx = 1;

for pId = 1:numel(tmp)
    if tmp(pId).name(1) == '.' || ~contains(tmp(pId).name, 'xdf')
        continue;
    end
    participants(inx).name   = tmp(pId).name;
    participants(inx).folder = tmp(pId).folder;
    participants(inx).date   = tmp(pId).date;
    inx = inx + 1;
end
clear tmp inx;

if isempty(participants)
    error('No participants found in %s', dataRoot);
end

%% === CONTAINERS FOR RESULTS =============================================
subjCount = struct('id', {}, 'times', {}, 'conds', {});
subjArea  = struct('id', {}, 'times', {}, 'conds', {});

%% === LOOP OVER PARTICIPANTS =============================================
for s = 1:numel(participants)

    pid_xdf = participants(s).name;
    pid     = pid_xdf(1:end-4);
    fprintf('\n=== %s ===\n', pid);

    tok = regexp(pid, 'Participant(\d+)', 'tokens', 'once');
    if isempty(tok)
        fprintf('  -> Could not parse participant number from %s, skipping.\n', pid);
        continue;
    end
    pidNum = str2double(tok{1});

    participantFolder = fullfile(dataRoot, ['preproc_' pid]);
    if ~exist(participantFolder, 'dir')
        fprintf('  -> Folder not found, skipping: %s\n', participantFolder);
        continue;
    end

    setFile = fullfile(participantFolder, ...
        sprintf('5b2_epochedFixationRejection_withTargets_vehicle_as_objects_%s.set', pid));

    if ~exist(setFile, 'file')
        fprintf('  -> 5b2 file not found, skipping: %s\n', setFile);
        continue;
    end

    EEG = pop_loadset(setFile);

    chanIdx = find(strcmpi({EEG.chanlocs.labels}, chanName), 1);
    if isempty(chanIdx)
        fprintf('  -> Channel %s not found, skipping participant.\n', chanName);
        continue;
    end

    times  = EEG.times;
    nEpoch = numel(EEG.epoch);

    fixIDeeg = getFixIDsFromEEGepochs(EEG);
    if isempty(fixIDeeg)
        fprintf('  -> Could not find a fixationID field in EEG.epoch for %s, skipping.\n', pid);
        fprintf('     Try printing: fieldnames(EEG.epoch)\n');
        continue;
    end

    Tsub = T(T.participant == pidNum, :);
    if isempty(Tsub)
        fprintf('  -> No CSV rows for participant %d (%s), skipping.\n', pidNum, pid);
        continue;
    end

    fixIDcsv = double(Tsub.('fixation id'));
    [isIn, loc] = ismember(double(fixIDeeg), fixIDcsv);

    personNr   = nan(nEpoch,1);
    personArea = nan(nEpoch,1);
    personNr(isIn)   = double(Tsub.person_nr(loc(isIn)));
    personArea(isIn) = double(Tsub.person_area(loc(isIn)));

    % ===================== COUNT ========================================
    lblCount = repmat({''}, nEpoch, 1);
    lblCount(personNr <= loCount) = {'LessCrowded'};
    lblCount(personNr >= hiCount) = {'Crowded'};

    nLoC = sum(strcmp(lblCount,'LessCrowded'));
    nHiC = sum(strcmp(lblCount,'Crowded'));
    fprintf('  COUNT labels (GLOBAL): LessCrowded=%d, Crowded=%d (min=%d)\n', ...
        nLoC, nHiC, minTrialsPerCond);

    condCount = computeTwoCondERP(EEG, chanIdx, lblCount, minTrialsPerCond);

    if ~isempty(condCount)
        subjCount(end+1).id    = pid; %#ok<AGROW>
        subjCount(end).times   = times;
        subjCount(end).conds   = condCount;

        fprintf('  COUNT accepted: LessCrowded=%d, Crowded=%d\n', ...
            condCount(1).nTrials, condCount(2).nTrials);

        fig = figure('Name', sprintf('ERPs %s - %s (Crowd by COUNT, GLOBAL)', pid, chanName), ...
            'Color','w');
        hold on;
        plot(times, condCount(1).erp, 'LineWidth', STYLE.LINE_W);
        plot(times, condCount(2).erp, 'LineWidth', STYLE.LINE_W);
        xline(0,'--k','HandleVisibility','off');

        xlim(STYLE.X_LIM);
        ylim(STYLE.Y_LIM);

        if STYLE.GRID_ON, grid on; end
        if STYLE.BOX_OFF, box off; end
        set(gca,'FontSize',STYLE.AX_FONTSZ);

        xlabel('Time (ms)', 'FontSize', STYLE.LAB_FONTSZ);
        ylabel('Amplitude (\muV)', 'FontSize', STYLE.LAB_FONTSZ);
        %title(sprintf('%s - %s | Crowd by COUNT (GLOBAL q=%.2f)', pid, chanName, qCrowd), ...
        %    'Interpreter','none');

        baseName = sprintf('ERPs_CrowdCount_GLOBAL_q%.2f_%s_%s', qCrowd, pid, chanName);
        outPng = fullfile(outDir, [baseName '.png']);
        outPdf = fullfile(outDir, [baseName '.pdf']);
        outFig = fullfile(outDir, [baseName '.fig']);
        exportgraphics(gcf, outPng, 'Resolution', 300);
        exportgraphics(gcf, outPdf, 'ContentType', 'vector');
        savefig(gcf, outFig);
        close(fig);

        fprintf('  -> Saved COUNT plot: %s\n', outPng);
    else
        fprintf('  -> COUNT: not enough trials (participant skipped for COUNT).\n');
    end

    % ===================== AREA =========================================
    lblArea = repmat({''}, nEpoch, 1);
    lblArea(personArea <= loArea) = {'LessCrowded'};
    lblArea(personArea >= hiArea) = {'Crowded'};

    nLoA = sum(strcmp(lblArea,'LessCrowded'));
    nHiA = sum(strcmp(lblArea,'Crowded'));
    fprintf('  AREA labels (GLOBAL): LessCrowded=%d, Crowded=%d (min=%d)\n', ...
        nLoA, nHiA, minTrialsPerCond);

    condArea = computeTwoCondERP(EEG, chanIdx, lblArea, minTrialsPerCond);

    if ~isempty(condArea)
        subjArea(end+1).id    = pid; %#ok<AGROW>
        subjArea(end).times   = times;
        subjArea(end).conds   = condArea;

        fprintf('  AREA accepted: LessCrowded=%d, Crowded=%d\n', ...
            condArea(1).nTrials, condArea(2).nTrials);

        fig = figure('Name', sprintf('ERPs %s - %s (Crowd by AREA, GLOBAL)', pid, chanName), ...
            'Color','w');
        hold on;
        plot(times, condArea(1).erp, 'LineWidth', STYLE.LINE_W);
        plot(times, condArea(2).erp, 'LineWidth', STYLE.LINE_W);
        xline(0,'--k','HandleVisibility','off');

        xlim(STYLE.X_LIM);
        ylim(STYLE.Y_LIM);

        if STYLE.GRID_ON, grid on; end
        if STYLE.BOX_OFF, box off; end
        set(gca,'FontSize',STYLE.AX_FONTSZ);

        xlabel('Time (ms)', 'FontSize', STYLE.LAB_FONTSZ);
        ylabel('Amplitude (\muV)', 'FontSize', STYLE.LAB_FONTSZ);
        %title(sprintf('%s - %s | Crowd by AREA (GLOBAL q=%.2f)', pid, chanName, qCrowd), ...
        %    'Interpreter','none');

        baseName = sprintf('ERPs_CrowdArea_GLOBAL_q%.2f_%s_%s', qCrowd, pid, chanName);
        outPng = fullfile(outDir, [baseName '.png']);
        outPdf = fullfile(outDir, [baseName '.pdf']);
        outFig = fullfile(outDir, [baseName '.fig']);
        exportgraphics(gcf, outPng, 'Resolution', 300);
        exportgraphics(gcf, outPdf, 'ContentType', 'vector');
        savefig(gcf, outFig);
        close(fig);

        fprintf('  -> Saved AREA plot: %s\n', outPng);
    else
        fprintf('  -> AREA: not enough trials (participant skipped for AREA).\n');
    end
end

%% === GRAND AVERAGES =====================================================
fprintf('\n=== GRAND AVERAGES (weighted by trial count) ===\n');

if ~isempty(subjCount)
    grandCount = grandAverageTwoCond(subjCount);

    fprintf('Grand COUNT summary:\n');
    fprintf('  LessCrowded: trials=%d, subjects=%d\n', ...
        grandCount.less.nTrials, grandCount.less.nSubj);
    fprintf('  Crowded    : trials=%d, subjects=%d\n', ...
        grandCount.more.nTrials, grandCount.more.nSubj);

    fig = figure('Name', sprintf('Grand Avg - %s (Crowd by COUNT, GLOBAL)', chanName), ...
        'Color','w');
    hold on;
    p1 = plot(grandCount.times, grandCount.less.erp, 'LineWidth', STYLE.LINE_W);
    p2 = plot(grandCount.times, grandCount.more.erp, 'LineWidth', STYLE.LINE_W);
    xline(0,'--k','HandleVisibility','off');

    xlim(STYLE.X_LIM);
    ylim(STYLE.Y_LIM);

    if STYLE.GRID_ON, grid on; end
    if STYLE.BOX_OFF, box off; end
    set(gca,'FontSize',STYLE.AX_FONTSZ);

    xlabel('Time (ms)', 'FontSize', STYLE.LAB_FONTSZ);
    ylabel('Amplitude (\muV)', 'FontSize', STYLE.LAB_FONTSZ);
    %title(sprintf('Grand Average - %s | Crowd by COUNT (GLOBAL q=%.2f)', chanName, qCrowd), ...
    %    'Interpreter','none');

    baseName = sprintf('ERPs_GrandAverage_CrowdCount_GLOBAL_q%.2f_%s', qCrowd, chanName);
    outPng = fullfile(outDir, [baseName '.png']);
    outPdf = fullfile(outDir, [baseName '.pdf']);
    outFig = fullfile(outDir, [baseName '.fig']);
    exportgraphics(gcf, outPng, 'Resolution', 300);
    exportgraphics(gcf, outPdf, 'ContentType', 'vector');
    savefig(gcf, outFig);
    close(fig);

    fprintf('  -> Saved GRAND COUNT: %s\n', outPng);

    if MAKE_LEGEND_FIG
        figLeg = figure('Color','w','Units','centimeters','Position',[2 2 6 4]);
        ax = axes(figLeg);
        hold(ax,'on');
        h1 = plot(ax, nan, nan, 'LineWidth', STYLE.LINE_W);
        h2 = plot(ax, nan, nan, 'LineWidth', STYLE.LINE_W);
        axis(ax,'off');
        %legend(ax, [h1 h2], ...
        %    {sprintf('Less crowded (trials=%d, subjects=%d)', grandCount.less.nTrials, grandCount.less.nSubj), ...
        %     sprintf('Crowded (trials=%d, subjects=%d)', grandCount.more.nTrials, grandCount.more.nSubj)}, ...
        %    'Location', 'northwest', 'Box', 'off');
        legend(ax, [h1 h2], ...
            {'Not crowded', 'Crowded'}, ...
            'Location', 'northwest', 'Box', 'off');

        legBase = sprintf('Legend_GrandAverage_CrowdCount_GLOBAL_q%.2f_%s', qCrowd, chanName);
        exportgraphics(figLeg, fullfile(outDir, [legBase '.png']), 'Resolution', 300);
        exportgraphics(figLeg, fullfile(outDir, [legBase '.pdf']), 'ContentType', 'vector');
        savefig(figLeg, fullfile(outDir, [legBase '.fig']));
        close(figLeg);

        fprintf('  -> Saved separate COUNT legend figure.\n');
    end
else
    fprintf('No usable subjects for COUNT grand average.\n');
end

if ~isempty(subjArea)
    grandArea = grandAverageTwoCond(subjArea);

    fprintf('\nGrand AREA summary:\n');
    fprintf('  LessCrowded: trials=%d, subjects=%d\n', ...
        grandArea.less.nTrials, grandArea.less.nSubj);
    fprintf('  Crowded    : trials=%d, subjects=%d\n', ...
        grandArea.more.nTrials, grandArea.more.nSubj);

    fig = figure('Name', sprintf('Grand Avg - %s (Crowd by AREA, GLOBAL)', chanName), ...
        'Color','w');
    hold on;
    p1 = plot(grandArea.times, grandArea.less.erp, 'LineWidth', STYLE.LINE_W);
    p2 = plot(grandArea.times, grandArea.more.erp, 'LineWidth', STYLE.LINE_W);
    xline(0,'--k','HandleVisibility','off');

    xlim(STYLE.X_LIM);
    ylim(STYLE.Y_LIM);

    if STYLE.GRID_ON, grid on; end
    if STYLE.BOX_OFF, box off; end
    set(gca,'FontSize',STYLE.AX_FONTSZ);

    xlabel('Time (ms)', 'FontSize', STYLE.LAB_FONTSZ);
    ylabel('Amplitude (\muV)', 'FontSize', STYLE.LAB_FONTSZ);
    %title(sprintf('Grand Average - %s | Crowd by AREA (GLOBAL q=%.2f)', chanName, qCrowd), ...
    %    'Interpreter','none');

    baseName = sprintf('ERPs_GrandAverage_CrowdArea_GLOBAL_q%.2f_%s', qCrowd, chanName);
    outPng = fullfile(outDir, [baseName '.png']);
    outPdf = fullfile(outDir, [baseName '.pdf']);
    outFig = fullfile(outDir, [baseName '.fig']);
    exportgraphics(gcf, outPng, 'Resolution', 300);
    exportgraphics(gcf, outPdf, 'ContentType', 'vector');
    savefig(gcf, outFig);
    close(fig);

    fprintf('  -> Saved GRAND AREA: %s\n', outPng);

    if MAKE_LEGEND_FIG
        figLeg = figure('Color','w','Units','centimeters','Position',[2 2 6 4]);
        ax = axes(figLeg);
        hold(ax,'on');
        h1 = plot(ax, nan, nan, 'LineWidth', STYLE.LINE_W);
        h2 = plot(ax, nan, nan, 'LineWidth', STYLE.LINE_W);
        axis(ax,'off');
        %legend(ax, [h1 h2], ...
        %    {sprintf('Less crowded (trials=%d, subjects=%d)', grandArea.less.nTrials, grandArea.less.nSubj), ...
        %     sprintf('Crowded (trials=%d, subjects=%d)', grandArea.more.nTrials, grandArea.more.nSubj)}, ...
        %    'Location', 'northwest', 'Box', 'off');
        legend(ax, [h1 h2], ...
            {'Not crowded', 'Crowded'}, ...
            'Location', 'northwest', 'Box', 'off');

        legBase = sprintf('Legend_GrandAverage_CrowdArea_GLOBAL_q%.2f_%s', qCrowd, chanName);
        exportgraphics(figLeg, fullfile(outDir, [legBase '.png']), 'Resolution', 300);
        exportgraphics(figLeg, fullfile(outDir, [legBase '.pdf']), 'ContentType', 'vector');
        savefig(figLeg, fullfile(outDir, [legBase '.fig']));
        close(figLeg);

        fprintf('  -> Saved separate AREA legend figure.\n');
    end
else
    fprintf('No usable subjects for AREA grand average.\n');
end

%% === SAVE RESULTS TO .MAT ===============================================
outMat = fullfile(outDir, sprintf('ERPs_Crowdedness_GLOBAL_q%.2f_%s.mat', qCrowd, chanName));
save(outMat, ...
    'subjCount','subjArea','chanName','minTrialsPerCond','qCrowd', ...
    'loCount','hiCount','loArea','hiArea');
fprintf('\nSaved results .mat to %s\n', outMat);

diary off;
fprintf('Saved Command Window log to %s\n', logFile);

%% ===================== LOCAL FUNCTIONS =================================
function fixIDs = getFixIDsFromEEGepochs(EEG)
    candidates = {'fixationID','fixation_id','fixationId','fix_id','fixID','fixation','fixationid'};
    fixIDs = [];

    for i = 1:numel(candidates)
        fn = candidates{i};
        if isfield(EEG.epoch, fn)
            tmp = {EEG.epoch.(fn)};
            fixIDs = cellfun(@toNum, tmp);
            return;
        end
    end
end

function x = toNum(v)
    if iscell(v), v = v{1}; end
    if isstring(v), v = char(v); end
    if ischar(v)
        x = str2double(v);
    else
        x = double(v);
    end
end

function condStruct = computeTwoCondERP(EEG, chanIdx, labels, minTrials)
    condStruct = struct('name', {}, 'erp', {}, 'nTrials', {});
    T = size(EEG.data,2);

    idxLo = strcmp(labels, 'LessCrowded');
    idxHi = strcmp(labels, 'Crowded');

    nLo = sum(idxLo);
    nHi = sum(idxHi);

    if nLo < minTrials || nHi < minTrials
        condStruct = [];
        return;
    end

    dataLo = EEG.data(chanIdx, :, idxLo);
    dataHi = EEG.data(chanIdx, :, idxHi);

    erpLo = mean(dataLo, 3);
    erpHi = mean(dataHi, 3);

    erpLo = reshape(erpLo, [1, T]);
    erpHi = reshape(erpHi, [1, T]);

    condStruct(1).name    = 'LessCrowded';
    condStruct(1).erp     = erpLo;
    condStruct(1).nTrials = nLo;

    condStruct(2).name    = 'Crowded';
    condStruct(2).erp     = erpHi;
    condStruct(2).nTrials = nHi;
end

function G = grandAverageTwoCond(subjData)
    times = subjData(1).times;

    numerLo = zeros(1, numel(times));
    denomLo = 0; nSubjLo = 0;

    numerHi = zeros(1, numel(times));
    denomHi = 0; nSubjHi = 0;

    for s = 1:numel(subjData)
        conds = subjData(s).conds;
        if numel(conds) < 2, continue; end

        numerLo = numerLo + conds(1).erp * conds(1).nTrials;
        denomLo = denomLo + conds(1).nTrials;
        nSubjLo = nSubjLo + 1;

        numerHi = numerHi + conds(2).erp * conds(2).nTrials;
        denomHi = denomHi + conds(2).nTrials;
        nSubjHi = nSubjHi + 1;
    end

    G.times = times;
    G.less.erp     = numerLo / denomLo;
    G.less.nTrials = denomLo;
    G.less.nSubj   = nSubjLo;

    G.more.erp     = numerHi / denomHi;
    G.more.nTrials = denomHi;
    G.more.nSubj   = nSubjHi;
end
