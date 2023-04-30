load 'group_ar_teased.mat'

set(groot,'defaultLineLineWidth',2.0)
figure()
xline(0,'k--','LineWidth',1)
xline(-5,'k--','LineWidth',1)
xline(5,'k--','LineWidth',1)
xline(-10,'k--','LineWidth',1)
xline(10,'k--','LineWidth',1)
xline(-15,'k--','LineWidth',1)
xline(15,'k--','LineWidth',1)
yline(.05,'k--','LineWidth',1)
yline(.1,'k--','LineWidth',1)
yline(.025,'k--','LineWidth',1)
yline(.125,'k--','LineWidth',1)
yline(.075,'k--','LineWidth',1)
hold on 
h1 = histogram(rare_card,'FaceColor','r') ; 
h1.Normalization = 'probability';
h1.BinWidth = .5;
h1.FaceAlpha = .9; 
h2 = histogram(freq_card,'FaceColor','k'); %have border between adjacent bins right at 0 same bin width 
h2.Normalization = 'probability';
h2.BinWidth = .5;
h2.FaceAlpha = .7; 
box off
ylabel('Probability')
xlabel('←repulsion  Signed Error  attraction→')
xticks([-20 -10 0 10 20 ])
xticklabels({'-20°','-10°','0°','10°','20°'})
title('Cardinal')
set(gca,'FontSize',20); 

set(gca,'LineWidth',2); 
xlim([-20 20])
ylim([0 .15])
set(gca,'TickDir','out'); % The only other option is 'in'
legend
legend box on
legend edgecolor 'none'
h = legend('', '' , '','','','','','','','','','','Rare','Frequent','FontSize',12);
% legend off 

saveas(figure(1),['DO_Card_RF'],'svg')

figure()
xline(0,'k--','LineWidth',1)
xline(-5,'k--','LineWidth',1)
xline(5,'k--','LineWidth',1)
xline(-10,'k--','LineWidth',1)
xline(10,'k--','LineWidth',1)
xline(-15,'k--','LineWidth',1)
xline(15,'k--','LineWidth',1)
yline(.05,'k--','LineWidth',1)
yline(.1,'k--','LineWidth',1)
yline(.025,'k--','LineWidth',1)
yline(.125,'k--','LineWidth',1)
yline(.075,'k--','LineWidth',1)
hold on 
h3 = histogram(rare_diag,'FaceColor','r') ; 
h3.Normalization = 'probability';
h3.BinWidth = .5;
h3.FaceAlpha = .8; 
h4 = histogram(freq_diag,'FaceColor','k'); %have border between adjacent bins right at 0 same bin width 
h4.Normalization = 'probability';
h4.BinWidth = .5;
h4.FaceAlpha = .7; 
box off
ylabel('Probability')
xlabel('←repulsion  Signed Error  attraction→')
xticks([-20 -10 0 10 20 ])
xticklabels({'-20°','-10°','0°','10°','20°'})
title('Diagonal')
set(gca,'FontSize',20); 
set(gca,'LineWidth',2); 
xlim([-20 20])
ylim([0 .15])
set(gca,'TickDir','out'); % The only other option is 'in'
legend
legend('', '' , '','','','','','','','','','','Rare','Frequent','FontSize',12)
legend box on 
legend edgecolor 'none'
% legend off 
saveas(figure(2),['DO_Diag_RF'],'svg')
clear; close all; 

%% attraction repulsion plots 

subList = 102:123; 

% freq first column rare second column 

degree_of_error = zeros(length(subList),4); 

for i = 1:length(subList)

    %% Parameters to set per subject 
    currentSub = num2str(subList(i));
    loadThis = strcat(currentSub,'_card_diag_attraction_repulsion.mat');
    load(loadThis)
    
    degree_of_error(i,1) = mean(freq_ar_card); 
    degree_of_error(i,2) = mean(rare_ar_card); 
    degree_of_error(i,3) = mean(freq_ar_diag); 
    degree_of_error(i,4) = mean(rare_ar_diag); 

    
end 

% group accuracys now 
Group_degree = (mean(degree_of_error,1)); 

Group_std_error = std(degree_of_error,1)/ sqrt(size(degree_of_error,1));

% breaking things down for nicer figures cardinal 
freq_card = [0,Group_degree(1)]; 
rare_card = [Group_degree(2),0];
% diagonal 
freq_diag = [0,Group_degree(3)]; 
rare_diag = [Group_degree(4),0];


% figure attraction repulsion 
figure()
hold on; 
bar(rare_card,.9,'r')
bar(freq_card,.9,'k')
yline(0,'k','LineWidth',2); 
er = errorbar(mean(degree_of_error(:,[2 1]),1),Group_std_error([2 1]),'b','LineWidth',2);
er.LineStyle = 'none';
title('Cardinal')
set(gca,'FontSize',20); 
xticks ([1:2]);
set(gca,'FontSize',20); 
ylabel('←repulsion  Signed Error  attraction→','FontSize',18)
xticks ([1:2]);
ylim([-1 1])
yticks([-1 -.5 0 .5 1])
yticklabels({'-1°','-0.5°','0°','0.5°','1°'})
set(gca,'LineWidth',2); 
set(findall(gca, 'Type', 'Line'),'LineWidth',2);
set(gca,'xtick',[])
set(gca,'TickDir','out'); % The only other option is 'in'
legend({'Rare','Frequent'},'FontSize',12)
legend box off 
saveas(figure(1),['DO_AR_Cardinal'],'svg')

% figure attraction repulsion 
figure()
hold on; 
bar(rare_diag,.9,'r')
bar(freq_diag,.9,'k')
yline(0,'k','LineWidth',2); 
er = errorbar(mean(degree_of_error(:,[4 3]),1),Group_std_error([4 3]),'b','LineWidth',2);
er.LineStyle = 'none';
title('Diagonal')
set(gca,'FontSize',20); 
ylabel('←repulsion  Signed Error  attraction→','FontSize',18)
xticks ([1:2]);
ylim([-1 1])
yticks([-1 -.5 0 .5 1])
yticklabels({'-1°','-0.5°','0°','0.5°','1°'})
set(gca,'LineWidth',2); 
set(findall(gca, 'Type', 'Line'),'LineWidth',2);
set(gca,'xtick',[])
legend({'Rare','Frequent'},'FontSize',12)
legend box off 
set(gca,'TickDir','out'); % The only other option is 'in'

saveas(figure(2),['DO_AR_Diagonal'],'svg')


%% 2 way anova statistical test 

%% anova 
%matrix manipulation 
rare_cd = [degree_of_error(:,2);degree_of_error(:,4)]; 
freq_cd = [degree_of_error(:,1);degree_of_error(:,3)]; 

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

%ttest 
[h,pc,ci,stats] = ttest(degree_of_error(:,2),degree_of_error(:,1)) 
[h,pd,ci,stats] = ttest(degree_of_error(:,4),degree_of_error(:,3)) 


%FDR Correction 
pvals = [pc,pd]; 

[~, ~, ~, fdr_test]=fdr_bh(pvals,0.05,'dep','yes');


