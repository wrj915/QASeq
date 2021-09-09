% CNV analysis example code

close all
clear all

% plotswitch = 1, 175-plex; plotswitch = 2, 223-plex; other, no plotting
plotswitch = 2;


% Alpha value; for 175-plex alpha = 0.05/12; for 223-plex alpha = 0.05/22
if plotswitch  == 1
    alphaval = 0.05/12;
elseif plotswitch  == 2
    alphaval = 0.05/22;
end

%% Read test sample data
% Need to specify filename, Sheet, and range
% Make sure this new file has the same order of amplicons as All_StdSamples.xlsx
readdata = readtable('example_UMI_count_input.xlsx','Sheet','group',...
    'Range','G2:J221',...
    'ReadVariableNames',false);
alldata = zeros(size(readdata,1),size(readdata,2));
for i = 1:size(readdata,1)
    for j = [1 3:size(readdata,2)]
        alldata(i,j) = readdata{i,j};
    end
end


%% Read standard sample data for calibration
% Need to specify Sheet and range for the correct sample type
% Here sample type is cfDNA + 223-plex
readdata = readtable('All_StdSamples.xlsx','Sheet','cf_StdSample_223',...
    'Range','G2:P221',...
    'ReadVariableNames',false);
stddata = zeros(size(readdata,1),size(readdata,2));
plex_group = cell(size(readdata,1),1);
for i = 1:size(readdata,1)
    for j = [1 3:size(readdata,2)]
        stddata(i,j) = readdata{i,j};
    end
    plex_group{i} = strsplit(readdata{i,2}{1},',');
end
refavg = mean(stddata(:,3:end),2);


all_norm_final = zeros(size(alldata,1), size(alldata,2)-2);
all_utest_info = {};

%% Example samples
for libnum = [1:2]
    std_count = refavg;
    
    % Find sample indexes you'd like to test (may need to change)
    sample_count = alldata(:,2+libnum);
    
    % U-test main function
    [Utestinfo plex_refornot] = UtestQASeq_allPvalv2(sample_count, std_count, plex_group, alphaval);
    
    % Re-calculate normalized yield for each plex
    norm_yield = sample_count./std_count * 2;
    
    ref_ind = find(plex_refornot == 1);
    ref_count = norm_yield(ref_ind);
    
    norm_final = norm_yield ./ median(ref_count) * 2;
    
    all_norm_final(:,libnum) = norm_final;
    
    %% Plotting section for 175-plex
    if plotswitch == 1
        figure
        hold all
        for tarnum = 1:length(norm_final)
            plot([tarnum tarnum], log2([2.00 norm_final(tarnum)]),'-b');
        end
        
        plot([0 180], [1 1], '-k');
        
        for tarnum = 2:length(norm_final)
            if ~strcmp(readdata{tarnum-1,2}{1},readdata{tarnum,2}{1})
                plot([tarnum-0.5 tarnum-0.5], [min([min(log2(norm_final))-0.2 -1]) max([max(log2(norm_final))+0.2 3])], '-', 'Color', [0.5 0.5 0.5]);
            end
        end
        
        set(gcf, 'Position', [0 200 480 70]);
        set(gca, 'XTick', []);
        set(gca, 'YTick', [-1:1:6]);
        set(gca, 'YTickLabel', 2.^[-1:1:6]);
        set(gca,'FontSize',8);
        
        axis([0 173 min([min(log2(norm_final))-0.2 -1]) max([max(log2(norm_final))+0.2 3])])
        
        filename = sprintf('Utest_norm_ploidy_lib%d.eps',libnum);
        saveas(gcf,filename,'epsc');
        
        
        hold all
        logpval = log10([Utestinfo{:,4}]);
        logpval([5 7 10]) = [];
        plot([1:12], logpval, '.b', 'MarkerSize',7);
        plot([0 12+1], log10(alphaval)*ones(1,2),'--k');
        axis([0 12+1 -22 0]);
        set(gcf, 'Position', [0 200 230 69.4]);
        
        ax = gca;
        ax.XTick = [1:12];
        ax.XTickLabel = {'ERBB4' 'PIK3CA' 'ESR1' 'EGFR' 'PTEN' 'ERBB3' 'BRCA2' 'TP53' 'ERBB2' 'BRCA1' 'Chr17p' 'Chr17'};
        ax.XTickLabelRotation = 90;
        set(gca,'FontSize',8);
        
        filename = sprintf('Utest_pval_lib%d.eps',libnum);
        saveas(gcf,filename,'epsc');
    end
    
    
    
    %% Plotting section for 223-plex
    if plotswitch == 2
        % Plot normalized ploidy at each amplicon
        figure
        hold all
        for tarnum = 1:length(norm_final)
            plot([tarnum tarnum], log2([2.00 norm_final(tarnum)]),'-b');
        end
        
        plot([0 240], [1 1], '-k');
        
        for tarnum = 2:length(norm_final)
            if ~strcmp(readdata{tarnum-1,2}{1},readdata{tarnum,2}{1})
                plot([tarnum-0.5 tarnum-0.5], [min([min(log2(norm_final))-0.2 -1]) max([max(log2(norm_final))+0.2 3])], '-', 'Color', [0.5 0.5 0.5]);
            end
        end
        
        set(gcf, 'Position', [0 200 480 70]);
        set(gca, 'XTick', []);
        set(gca, 'YTick', [-1:1:6]);
        set(gca, 'YTickLabel', 2.^[-1:1:6]);
        set(gca,'FontSize',8);
        
        axis([0 221 min([min(log2(norm_final))-0.2 -1]) max([max(log2(norm_final))+0.2 3])])
        
        filename = sprintf('Utest_norm_ploidy_lib%d.eps',libnum);
        saveas(gcf,filename,'epsc');
        
        figure
        hold all
        logpval = log10([Utestinfo{:,4}]);
        logpval([9 15 18]) = [];% 9, 15, 18 are not effective groups
        plot([1:22], logpval, '.b', 'MarkerSize',7);
        plot([0 22+1], log10(alphaval)*ones(1,2),'--k');
        axis([0 22+1 min([-22 min(logpval)-0.2]) 0]);
        set(gcf, 'Position', [0 200 230 69.4]);
        
        ax = gca;
        ax.XTick = [1:22];
        ax.XTickLabel = {'NCSTN' 'PTGS2' 'AKT3' 'Chr1' 'ERBB4' 'PIK3CA' 'ESR1' 'EGFR' 'FGFR1' 'NBN' 'MYC' 'Chr8' 'PTEN' 'ERBB3' 'BRCA2' 'TP53' 'ERBB2' 'BRCA1' 'RAD51C' 'BRIP1' 'Chr17p' 'Chr17'};
        ax.XTickLabelRotation = 90;
        set(gca,'FontSize',8);
        
        filename = sprintf('Utest_pval_lib%d.eps',libnum);
        saveas(gcf,filename,'epsc');
    end
    
    all_utest_info = [all_utest_info; Utestinfo];
end

% Summary of ploidy values
allploidy = reshape([all_utest_info{:,6}],size(Utestinfo,1),size(all_utest_info,1)/size(Utestinfo,1));

% Summary of p-values
allpval = reshape([all_utest_info{:,4}],size(Utestinfo,1),size(all_utest_info,1)/size(Utestinfo,1));