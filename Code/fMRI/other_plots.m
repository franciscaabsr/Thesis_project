% Some plots


figure('defaultaxesfontsize',12)
%subplot(2,1,1)
plot(1:470,0.5+double(dFC_labels(2,1:470) == 1),1:470,2+double(dFC_labels(2,1:470) == 2),1:470,3.5+double(dFC_labels(2,1:470) == 3),1:470,5+double(dFC_labels(2,1:470) == 4),'linewidth',2)
ylim([0 6.5])
xlim([0 470])
set(gca,'ytick',[1 2.5 4 5.5],'yticklabel',{'FC State 1','FC State 2','FC State 3','FC State 4'})
title('Cluster time-course after label smoothing','fontsize',14)
xlabel('Time (TR)','fontsize',12)

dFC_labels2 = Kmeans_results{2}.IDX;

%subplot(2,1,2)
plot(1:470,0.5+double(dFC_labels2(1:470) == 1),1:470,2+double(dFC_labels2(1:470) == 2),1:470,3.5+double(dFC_labels2(1:470) == 3),1:470,5+double(dFC_labels2(1:470) == 4),'linewidth',2)
ylim([0 6.5])
xlim([0 470])
set(gca,'ytick',[1 2.5 4 5.5],'yticklabel',{'FC State 1','FC State 2','FC State 3','FC State 4'})
title('Cluster time-course before label smoothing','fontsize',14)
xlabel('Time (TR)','fontsize',12)

%%

figure;
subplot(2,2,1)
plot(relative_delta_allTR,'*')
hold on
line([470 470],[0 1],'Color','k','linewidth',2,'LineStyle','--')
line([2*470 2*470],[0 1],'Color','k','linewidth',2,'LineStyle','--')
line([3*470 3*470],[0 1],'Color','k','linewidth',2,'LineStyle','--')
line([4*470 4*470],[0 1],'Color','k','linewidth',2,'LineStyle','--')
line([5*470 5*470],[0 1],'Color','k','linewidth',2,'LineStyle','--')
line([6*470 6*470],[0 1],'Color','k','linewidth',2,'LineStyle','--')
line([7*470 7*470],[0 1],'Color','k','linewidth',2,'LineStyle','--')
line([8*470 8*470],[0 1],'Color','k','linewidth',2,'LineStyle','--')
hold off
xlabel('Time (TR)','fontsize',18)
ylabel('Relative power','fontsize',18)
title('Delta','fontsize',20)
set(gca,'fontsize',16)
xlim([0 4230])
ylim([0 1])
subplot(2,2,2)
plot(relative_theta_allTR,'*')
hold on
line([470 470],[0 1],'Color','k','linewidth',2,'LineStyle','--')
line([2*470 2*470],[0 1],'Color','k','linewidth',2,'LineStyle','--')
line([3*470 3*470],[0 1],'Color','k','linewidth',2,'LineStyle','--')
line([4*470 4*470],[0 1],'Color','k','linewidth',2,'LineStyle','--')
line([5*470 5*470],[0 1],'Color','k','linewidth',2,'LineStyle','--')
line([6*470 6*470],[0 1],'Color','k','linewidth',2,'LineStyle','--')
line([7*470 7*470],[0 1],'Color','k','linewidth',2,'LineStyle','--')
line([8*470 8*470],[0 1],'Color','k','linewidth',2,'LineStyle','--')
hold off
xlabel('Time (TR)','fontsize',18)
ylabel('Relative power','fontsize',18)
title('Theta','fontsize',20)
set(gca,'fontsize',16)
xlim([0 4230])
ylim([0 1])
subplot(2,2,3)
plot(relative_alpha_allTR,'*')
hold on
line([470 470],[0 1],'Color','k','linewidth',2,'LineStyle','--')
line([2*470 2*470],[0 1],'Color','k','linewidth',2,'LineStyle','--')
line([3*470 3*470],[0 1],'Color','k','linewidth',2,'LineStyle','--')
line([4*470 4*470],[0 1],'Color','k','linewidth',2,'LineStyle','--')
line([5*470 5*470],[0 1],'Color','k','linewidth',2,'LineStyle','--')
line([6*470 6*470],[0 1],'Color','k','linewidth',2,'LineStyle','--')
line([7*470 7*470],[0 1],'Color','k','linewidth',2,'LineStyle','--')
line([8*470 8*470],[0 1],'Color','k','linewidth',2,'LineStyle','--')
hold off
xlabel('Time (TR)','fontsize',18)
ylabel('Relative power','fontsize',18)
title('ASI','fontsize',20)
set(gca,'fontsize',16)
xlim([0 4230])
ylim([0 1])
subplot(2,2,4)
plot(relative_beta_allTR,'*')
hold on
line([470 470],[0 1],'Color','k','linewidth',2,'LineStyle','--')
line([2*470 2*470],[0 1],'Color','k','linewidth',2,'LineStyle','--')
line([3*470 3*470],[0 1],'Color','k','linewidth',2,'LineStyle','--')
line([4*470 4*470],[0 1],'Color','k','linewidth',2,'LineStyle','--')
line([5*470 5*470],[0 1],'Color','k','linewidth',2,'LineStyle','--')
line([6*470 6*470],[0 1],'Color','k','linewidth',2,'LineStyle','--')
line([7*470 7*470],[0 1],'Color','k','linewidth',2,'LineStyle','--')
line([8*470 8*470],[0 1],'Color','k','linewidth',2,'LineStyle','--')
hold off
xlabel('Time (TR)','fontsize',18)
ylabel('Relative power','fontsize',18)
title('Beta','fontsize',20)
set(gca,'fontsize',16)
xlim([0 4230])
ylim([0 1])
    
%%

figure('defaultaxesfontsize',16)
plot(ASI_allTR_norm,'*')
hold on
line([470 470],[0 1],'Color','k','linewidth',2,'LineStyle','--')
line([2*470 2*470],[0 1],'Color','k','linewidth',2,'LineStyle','--')
line([3*470 3*470],[0 1],'Color','k','linewidth',2,'LineStyle','--')
line([4*470 4*470],[0 1],'Color','k','linewidth',2,'LineStyle','--')
line([5*470 5*470],[0 1],'Color','k','linewidth',2,'LineStyle','--')
line([6*470 6*470],[0 1],'Color','k','linewidth',2,'LineStyle','--')
line([7*470 7*470],[0 1],'Color','k','linewidth',2,'LineStyle','--')
line([8*470 8*470],[0 1],'Color','k','linewidth',2,'LineStyle','--')
hold off
xlabel('Time (TR)','fontsize',18)
ylabel('Normalized ASI','fontsize',18)
xlim([0 4230])
ylim([0 1])