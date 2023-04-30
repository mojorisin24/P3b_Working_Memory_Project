clear; clc; 
addpath([pwd, '/Decoding_Oddball_Prepped'])
load time.mat 
% Plotting Data 
subList = 1:22; 


Nsub = length(subList);

Nblock = 3; % cross-validation
Nitr = 100; % iteration
Ntp = 100; % # of time points
NBins = 4; % # of location bin
tm = time(1,1):20:time(1,end);

%create empty matrix
AverageAccuracy = nan(Nsub,Ntp);
    
for cond = 1:4
%1 = rare cardinal
%2 = frequent cardinal 
%3 = rare diagonal
%4 = frequent diagonal 

for sub = 1:Nsub
    DecodingAccuracy = nan(Ntp,Nblock,Nitr);
    
     %% load SVM_ECOC output files
    fileLocation = ([pwd,'/Disc_Results_Angles']);
    
    if cond == 1 
        readThis =strcat(fileLocation,'/Disc_Rare_Cardinal_Decoding_',num2str(subList(sub)),'.mat');
    elseif cond == 2 
        readThis =strcat(fileLocation,'/Disc_Freq_Cardinal_Decoding_',num2str(subList(sub)),'.mat');
    elseif cond == 3 
        readThis =strcat(fileLocation,'/Disc_Rare_Diagonal_Decoding_',num2str(subList(sub)),'.mat');
    elseif cond == 4 
        readThis =strcat(fileLocation,'/Disc_Freq_Diagonal_Decoding_',num2str(subList(sub)),'.mat');
    end 
 
    load(readThis)
     

    % prediciton from SVM-ECOC model
    svmPrediction = squeeze(svmECOC.modelPredict);
    tstTargets = squeeze(svmECOC.targets);
    clear svmECOC
    
    % compute decoding accuracy of each decoding trial
    for block = 1:Nblock
        for itr = 1:Nitr
            for tp = 1:Ntp  

                prediction = squeeze(svmPrediction(itr,tp,block,:)); % this is predictions from models
                TrueAnswer = squeeze(tstTargets(itr,tp,block,:)); % this is predictions from models
                Err = TrueAnswer - prediction;
                ACC = mean(Err==0);
                DecodingAccuracy(tp,block,itr) = ACC; % average decoding accuracy

            end
        end
    end
      
     %average across block and iterations
     grandAvg = squeeze(mean(mean(DecodingAccuracy,2),3));
     
     % Save smoothed data
     AverageAccuracy(sub,:) = grandAvg; % average across iteration and block
   
end %End of subject

if cond == 1 
    rare_card = AverageAccuracy; 
    %std dev
    rare_c_std = std(rare_card); 
    %std err
    rare_c_se = rare_c_std/sqrt(size(rare_card,1)); 

elseif cond == 2 
    freq_card = AverageAccuracy; 
    %std dev
    freq_c_std = std(freq_card); 
    %std err
    freq_c_se = freq_c_std/sqrt(size(freq_card,1)); 
elseif cond == 3 
    rare_diag = AverageAccuracy; 
    %std dev
    rare_d_std = std(rare_diag); 
    %std err
    rare_d_se = rare_d_std/sqrt(size(rare_diag,1)); 
elseif cond == 4 
    freq_diag = AverageAccuracy; 
    %std dev
    freq_d_std = std(freq_diag); 
    %std err
    freq_d_se = freq_d_std/sqrt(size(freq_diag,1)); 
    
end 
    

end 

%%fig

% %% Stats 

rare_c_pval = FDR_Correction(rare_card); 
freq_c_pval = FDR_Correction(freq_card); 


rare_d_pval = FDR_Correction(rare_diag); 
freq_d_pval = FDR_Correction(freq_diag); 

[C_between_pval,C_between_pval_uncorrected ] = FDR_Groups_Test(rare_card,freq_card); 
 
[D_between_pval, D_between_pval_uncorrected] = FDR_Groups_Test(rare_diag,freq_diag); 


figure()
cl=colormap(hot(50));
hold on 


%Drawing Significance

%uncorrected Pvalues 
accEst = mean(rare_card); 
a = area(1:(size(rare_card,2)),100.*C_between_pval_uncorrected');
a.EdgeColor = 'none';
a.FaceColor = [0.9,0.9,0.9];
child = get(a,'Children');
set(a,'FaceAlpha',0.9)

%adding shading to bottom of plot 
b = area(1:(size(rare_card,2)),-100.*C_between_pval_uncorrected');
b.EdgeColor = 'none';
b.FaceColor = [0.9,0.9,0.9];
child = get(b,'Children');
set(child,'FaceAlpha',0.9)

% Corrected P values
accEst = mean(rare_card); 
a = area(1:(size(rare_card,2)),100.*C_between_pval');
a.EdgeColor = 'none';
a.FaceColor = [0.3010 0.7450 0.9330];
child = get(a,'Children');
set(a,'FaceAlpha',0.4)

%adding shading to bottom of plot 
b = area(1:(size(rare_card,2)),-100.*C_between_pval');
b.EdgeColor = 'none';
b.FaceColor = [0.3010 0.7450 0.9330];
child = get(b,'Children');
set(child,'FaceAlpha',0.4)



rare_cardinal = boundedline(1:(size(rare_card,2)),((mean(rare_card))),rare_c_se, 'cmap',cl(15,:),'alpha','transparency',0.35);
%signficance
plot(rare_c_pval' * .125,'.','MarkerSize',20,...
    'MarkerEdgeColor','r',...
    'MarkerFaceColor',[1 .6 .6])

freq_cardinal = boundedline(1:(size(freq_card,2)),((mean(freq_card))),freq_c_se, 'cmap',cl(1,:),'alpha','transparency',0.35);
%significance
plot(freq_c_pval' * .145,'.','MarkerSize',20,...
    'MarkerEdgeColor','k',...
    'MarkerFaceColor',[1 .6 .6])
set(0,'DefaultAxesTitleFontWeight','normal');
xlabel('Time (ms)');ylabel('Decoding Accuracy')
yline(0.25, '--k','LineWidth',2); 
xline(26, '--k','LineWidth',2);
set(gca,'FontSize',20); 
set(findall(gca, 'Type', 'Line'),'LineWidth',2);
set(gca,'linewidth',2)
ylim([.1 .6])
xticks ([0 26 51 76 100]);
xticklabels ({'-500','0','500','1000','1500'});
xlim([0 100])
title('Cardinal')
set(gca,'FontWeight','normal')
box off





ylim([.1 .6])
set(gca,'TickDir','out'); % The only other option is 'in'


h = legend; 
set(h,'FontSize',13);
legend('', '' , '', '','','Rare','p < 0.05, FDR Corrected','','Frequent','p < 0.05, FDR Corrected','','')
legend box off
saveas(figure(1),['DO_Decoding_Angles_Cardinal_SigDiff_Between'],'svg')

%% Diagonal 
figure()
cl=colormap(hot(50));
hold on 
%uncorrected Pvalues 
accEst = mean(rare_card); 
a = area(1:(size(rare_card,2)),100.*D_between_pval_uncorrected');
a.EdgeColor = 'none';
a.FaceColor = [0.9,0.9,0.9];
child = get(a,'Children');
set(a,'FaceAlpha',0.9)

%adding shading to bottom of plot 
b = area(1:(size(rare_card,2)),-100.*D_between_pval_uncorrected');
b.EdgeColor = 'none';
b.FaceColor = [0.9,0.9,0.9];
child = get(b,'Children');
set(child,'FaceAlpha',0.9)

% Corrected P values
accEst = mean(rare_card); 
a = area(1:(size(rare_card,2)),100.*D_between_pval');
a.EdgeColor = 'none';
a.FaceColor = [0.3010 0.7450 0.9330];
child = get(a,'Children');
set(a,'FaceAlpha',0.4)

%adding shading to bottom of plot 
b = area(1:(size(rare_card,2)),-100.*D_between_pval');
b.EdgeColor = 'none';
b.FaceColor = [0.3010 0.7450 0.9330];
child = get(b,'Children');
set(child,'FaceAlpha',0.4)


rare_diagonal = boundedline(1:(size(rare_diag,2)),((mean(rare_diag))),rare_d_se, 'cmap',cl(15,:),'alpha','transparency',0.35);
%significance
plot(rare_d_pval' * .125,'.','MarkerSize',20,...
    'MarkerEdgeColor','r',...
    'MarkerFaceColor',[1 .6 .6])

freq_diagonal = boundedline(1:(size(freq_diag,2)),((mean(freq_diag))),freq_d_se, 'cmap',cl(1,:),'alpha','transparency',0.35);
%significance
plot(freq_d_pval' * .145,'.','MarkerSize',20,...
    'MarkerEdgeColor','k',...
    'MarkerFaceColor',[1 .6 .6])

xlabel('Time (ms)');ylabel('Decoding Accuracy')
yline(0.25, '--k','LineWidth',2); 
xline(26, '--k','LineWidth',2);
set(gca,'FontSize',20); 
set(findall(gca, 'Type', 'Line'),'LineWidth',2);
set(gca,'linewidth',2)
ylim([.1 .6])
xticks ([0 26 51 76 100]);
xticklabels ({'-500','0','500','1000','1500'});
xlim([0 100])
title('Diagonal')
set(gca,'FontWeight','normal')
box off

ylim([0.1 .6])
h = legend; 
set(h,'FontSize',13);
legend('', '' , '', '','','Rare','p < 0.05, FDR Corrected','','Frequent','p < 0.05, FDR Corrected','','')
legend box off 
set(gca,'TickDir','out'); % The only other option is 'in'

saveas(figure(2),['DO_Decoding_Angles_Diagonal_SigDiff'],'svg')

close all; 

%% Mean Decoding Accuracy SigDiff calculation 
rare_card_mean = mean(rare_card(:,44:60),2); %starting where signals deviate in univariate 360 ms 
freq_card_mean = mean(freq_card(:,44:60),2); 

rare_diag_mean = mean(rare_diag(:,44:60),2); %starting where signals deviate in univariate 360 ms 
freq_diag_mean = mean(freq_diag(:,44:60),2); 

%% Anova test
%% anova 
%matrix manipulation 
rare_cd = [rare_card_mean;rare_diag_mean]; 
freq_cd = [freq_card_mean;freq_diag_mean]; 

anova_matrix = [rare_cd,freq_cd]; 

[p,tbl,stats] = anova2(anova_matrix,22);

c1 = multcompare(stats);
title 'Rare vs Frequent' 


tbl1 = array2table(c1,"VariableNames", ...
    ["Rare","Frequent","Lower Limit","Estimated Diff","Upper Limit","P-value"]); 

c2 = multcompare(stats,"Estimate","row");
title 'Cardinal vs Diagonal' 
tbl2 = array2table(c2,"VariableNames", ...
    ["Cardinal","Diagonal","Lower Limit","Estimated Diff","Upper Limit","P-value"]);

close all; 
%% T Tests 

[h,pc,ci,stats] = ttest(rare_card_mean,freq_card_mean)

[hd,pd,cid,statsd] = ttest(rare_diag_mean,freq_diag_mean)

% chance tests 
[h,prc,ci,stats] = ttest(rare_card_mean,0)
[h,pfc,ci,stats] = ttest(freq_card_mean,0)

[h,prd,ci,stats] = ttest(rare_diag_mean,0)
[h,prf,ci,stats] = ttest(freq_diag_mean,0)

%FDR Correction 
pvals = [pc,pd,prc,pfc,prd,prf]; 

 [~, ~, ~,  ]=fdr_bh(pvals,0.05,'dep','yes');





%plots 

%cardinal 
cardinal_means = [rare_card_mean,freq_card_mean]; 


rare_mean_c = [mean(cardinal_means(:,1)),0]; 
freq_mean_c = [0,mean(cardinal_means(:,2))];
figure()
hold on
bar(freq_mean_c,.9,'k')
bar(rare_mean_c,.9,'r')
dec_err = std(cardinal_means)/ sqrt(size(cardinal_means,1));
er = errorbar(mean(cardinal_means(:,[1 2]),1),dec_err([1 2]),'b','LineWidth',2);
er.LineStyle = 'none';
xticks ([1:2]);
xticklabels ({'Rare', 'Frequent'});
ylabel('Decoding Accuracy')
title('Cardinal')
set(gca,'FontSize',20); 
set(findall(gca, 'Type', 'Line'),'LineWidth',2);
set(gca,'LineWidth',2); 
box off 
yline(0.25, '--b','LineWidth',2);
ylim([0,.5])
set(gca,'TickDir','out'); % The only other option is 'in'

saveas(figure(1),['DO_C_Mean'],'svg')
close all; 

%diagonal
diagonal_means = [rare_diag_mean,freq_diag_mean]; 


rare_mean_d = [mean(diagonal_means(:,1)),0]; 
freq_mean_d = [0,mean(diagonal_means(:,2))];
figure()
hold on
bar(freq_mean_d,.9,'k')
bar(rare_mean_d,.9,'r')
dec_err_d = std(diagonal_means)/ sqrt(size(diagonal_means,1));
er = errorbar(mean(diagonal_means(:,[1 2]),1),dec_err_d([1 2]),'b','LineWidth',2);
er.LineStyle = 'none';
xticks ([1:2]);
xticklabels ({'Rare', 'Frequent'});
ylabel('Decoding Accuracy')
title('Diagonal')
set(gca,'FontSize',20); 
set(findall(gca, 'Type', 'Line'),'LineWidth',2);
set(gca,'LineWidth',2); 
box off 
yline(0.25, '--b','LineWidth',2); 
ylim([0,.5])
set(gca,'TickDir','out'); % The only other option is 'in'

saveas(figure(1),['DO_D_Mean'],'svg')
close all; 

