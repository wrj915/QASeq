% UMI filter and plot
% remove UMIs w/ small number of counts; threshold determined by dynamic cutoff, top3 UMI family size does not distinguish different genotypes

clear all;
close all;
uniqUMInum1 = zeros(1,226);
uniqUMInum2 = zeros(1,226);
uniqUMInum3 = zeros(1,226);
uniqUMInum4 = zeros(1,226);
alignedreads = zeros(27,226);

for libnum = 1
    for tarnum = [1 2] %may change to 1:179 or 1:226
        fprintf('lib %d target %d\n',libnum,tarnum);
        umifilename = sprintf('UMIquantlib%d_sortedtar%d.txt',libnum,tarnum);
        fid = fopen(umifilename);
        if fid ~= -1
            tline = fgets(fid);
            linecount = 1;
            UMIinfo = {};
            while ischar(tline)
                spaceind = find(tline == '	', 1);
                UMIinfo{linecount,1} = tline(1:spaceind-1);
                UMIinfo{linecount,2} = str2num(tline(spaceind+1:end-1));
                tline = fgets(fid);
                linecount = linecount + 1;
            end
            fclose(fid);
            
            alignedreads(libnum,tarnum) = sum([UMIinfo{:,2}]);
            uniqUMInum1(libnum,tarnum) = size(UMIinfo,1);
            
            UMIinfo2 = UMIinfo;
            uc1 = 1;
            while uc1 <= size(UMIinfo2,1) % remove G
                umi1 = UMIinfo2{uc1,1};
                if ~isempty(find(umi1=='G'))
                    UMIinfo2(uc1,:) = [];
                else
                    uc1 = uc1 + 1;
                end
            end
            uniqUMInum2(libnum,tarnum) = size(UMIinfo2,1);
            
            UMIinfo3 = UMIinfo2;
            uc1 = 1;
            while uc1 <= size(UMIinfo3,1) % remove <=3 count UMIs
                umict = UMIinfo3{uc1,2};
                if umict <= 3
                    UMIinfo3(uc1,:) = [];
                else
                    uc1 = uc1 + 1;
                end
            end
            uniqUMInum3(libnum,tarnum) = size(UMIinfo3,1);
            
            % ADDED 20200413
            UMIinfo4 = UMIinfo3;
            uc1 = 1;
            if size(UMIinfo4,1) > 2 % if less than 3 UMI left, dynamic cutoff will not apply
                Top3 = [UMIinfo4{1:3,2}];
                
                % 5% of TOP average3, can be adjusted
                Threshold = floor(mean(Top3)*0.05);
                
                Threshold2 = max(3,Threshold);
                
                while uc1 <= size(UMIinfo4,1) % remove <=Threshold2 count UMIs
                    umict = UMIinfo4{uc1,2};
                    if umict <= Threshold2
                        UMIinfo4(uc1,:) = [];
                    else
                        uc1 = uc1 + 1;
                    end
                end
            end
            uniqUMInum4(libnum,tarnum) = size(UMIinfo4,1);
        end
        
%         % plotting histogram
%         figure
%         hold all
%         histogram(log10([UMIinfo2{:,2}]));
%         plot(log10([Threshold2 Threshold2]), [0 2000],'-r','LineWidth',1);
%         set(gcf,'Position',[200 200 200 170])
%         set(gca,'XTick',log10([1 3 10 30 100 300]))
%         xlim([0 2.4])
    end
    
end

csvwrite('allcountsDynamic.csv',uniqUMInum4);
