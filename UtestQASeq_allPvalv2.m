function [Utestinfo2 plex_refornot] = UtestQASeq_allPvalv2(sample_count, std_count, plex_group, alphaval)
% Similar to UtestQASeq, but calculate all p-values

% sample_count: raw counts at each plex
% std_count: counts/yield/ploidy at each plex, mean of standard samples
% plex_group: which group does each plex fit in
% column1: group names (need to be numbers)
% sample_count, std_count, and plex_group need to be in the same order

norm_yield = sample_count./std_count * 2;

% columns: group name; group size; list of plex #; current p value; ref or not (1 for ref, 0 for not); ploidy
allgroupinfo = cell(1,6);
allgroupinfo{1} = plex_group{1}{1};
allgroupinfo{2} = 0;
allgroupinfo{3} = [];
allgroupinfo{4} = -1;
allgroupinfo{5} = 1;
curtotalgroup = 1;

for curplex = 1:size(plex_group,1)
    for temp1 = 1:length(plex_group{curplex})
        findflag = 0;
        for temp2 = 1:size(allgroupinfo,1)
            if strcmp(allgroupinfo{temp2,1},plex_group{curplex}{temp1})
                findflag = 1;
                allgroupinfo{temp2,2} = allgroupinfo{temp2,2}+1;
                allgroupinfo{temp2,3} = [allgroupinfo{temp2,3} curplex];
                break
            end
        end
        
        if findflag == 0
            curtotalgroup = curtotalgroup+1;
            allgroupinfo{curtotalgroup,1} = plex_group{curplex}{temp1};
            allgroupinfo{curtotalgroup,2} = 1;
            allgroupinfo{curtotalgroup,3} = curplex;
            allgroupinfo{curtotalgroup,4} = -1;
            allgroupinfo{curtotalgroup,5} = 1;
        end
    end
end
allgroupinfo = sortrows(allgroupinfo,-2);

% record which plexes are used as ref (1 for ref, 0 for not)
plex_refornot = ones(size(plex_group,1),1);


curtestgroup = 1; % current group rank in allgroupinfo
while curtestgroup <= size(allgroupinfo,1) && allgroupinfo{curtestgroup,2} >= 3 % cannot U-test if group size <3!
    
    test_data = norm_yield(allgroupinfo{curtestgroup,3});
    
    plex_refornot_temp = plex_refornot;
    plex_refornot_temp(allgroupinfo{curtestgroup,3}) = zeros(allgroupinfo{curtestgroup,2},1);
    ref_ind = find(plex_refornot_temp == 1);
    ref_data = norm_yield(ref_ind);
    
    [pval,hstat] = ranksum(test_data,ref_data,'alpha',alphaval);
    allgroupinfo{curtestgroup,4} = pval;
    
    if hstat == 1 % test set is different from ref
        % remove current test group from pool
        allgroupinfo{curtestgroup,5} = 0;
        plex_refornot = plex_refornot_temp;
        
        % next test group is largest group in the ref pool
        %   if part of the group is ref, can still test this group
        for curtestgroup = 1:size(allgroupinfo,1)
            if allgroupinfo{curtestgroup,5} == 1
                break
            end
        end
        
        
    elseif hstat == 0 % test set not different from ref
        % next test group is next largest group in the ref pool
        for curtestgroup = curtestgroup+1:size(allgroupinfo,1)
            if allgroupinfo{curtestgroup,5} == 1
                break
            end
        end
    else
        error('Wrong hstat\n');
    end
    
end

for curgroup = 1:size(allgroupinfo,1)
    if allgroupinfo{curgroup,5} == 1 % not significant groups
        allgroupinfo{curgroup,6} = 2; % ploidy = 2
    elseif allgroupinfo{curgroup,5} == 0 % significant groups
        test_data = norm_yield(allgroupinfo{curgroup,3});
        ref_ind = find(plex_refornot == 1);
        ref_data = norm_yield(ref_ind);
        [pval,hstat] = ranksum(test_data,ref_data,'alpha',alphaval); % recalculate pval
        allgroupinfo{curgroup,4} = pval;
        
        if hstat == 1 % test set is still different from ref
            allgroupinfo{curgroup,6} = median(test_data)/median(ref_data)*2;
            
        else
            allgroupinfo{curgroup,6} = 2;
            
        end
    else
        error('Wrong ref info\n');
    end
end

Utestinfo2 = sortrows(allgroupinfo,1);