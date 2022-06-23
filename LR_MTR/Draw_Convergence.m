load error_GTDA.mat;
load error_TLDE.mat;
load error_TLPP.mat;
load error_TNPE.mat;
load error_TPCA.mat;


M=(1:15);
MM=(1:2:15);

figure(1);
semilogy(M,error_GTDA,'-or','LineWidth',2,...
    'MarkerSize',5,...
    'MarkerFaceColor','r');
xlabel('Iteration Number  ','FontName','Times New Rom an','FontSize',30);
ylabel('Error','FontName','Times New Roman','FontSize',30);
% legend('\it N_P','\it N_N','\it OF','\it N_E','FontName','Times New Roman','Location','North','Orientation','horizontal');
set(gca,'xtick',MM,'xticklabel',MM) 
set(gca,'FontName','Times New Roman','FontSize',20);
grid on
grid minor

% figure(2);
% plot(M,error_TLDE,'-or','LineWidth',2,...
%     'MarkerSize',5,...
%     'MarkerFaceColor','r');
% xlabel('Iteration Number  ','FontName','Times New Rom an','FontSize',30);
% ylabel('Error','FontName','Times New Roman','FontSize',30);
% % legend('\it N_P','\it N_N','\it OF','\it N_E','FontName','Times New Roman','Location','North','Orientation','horizontal');
% set(gca,'xtick',M,'xticklabel',M) 
% set(gca,'FontName','Times New Roman','FontSize',15);
% grid on
% grid minor
% 
% figure(3);
% plot(M,error_TLPP,'-or','LineWidth',2,...
%     'MarkerSize',5,...
%     'MarkerFaceColor','r');
% xlabel('Iteration Number  ','FontName','Times New Rom an','FontSize',30);
% ylabel('Error','FontName','Times New Roman','FontSize',30);
% % legend('\it N_P','\it N_N','\it OF','\it N_E','FontName','Times New Roman','Location','North','Orientation','horizontal');
% set(gca,'xtick',M,'xticklabel',M) 
% set(gca,'FontName','Times New Roman','FontSize',15);
% grid on
% grid minor

% figure(4);
% plot(M,error_TNPE,'-or','LineWidth',2,...
%     'MarkerSize',5,...
%     'MarkerFaceColor','r');
% xlabel('Iteration Number  ','FontName','Times New Rom an','FontSize',30);
% ylabel('Error','FontName','Times New Roman','FontSize',30);
% % legend('\it N_P','\it N_N','\it OF','\it N_E','FontName','Times New Roman','Location','North','Orientation','horizontal');
% set(gca,'xtick',M,'xticklabel',M) 
% set(gca,'FontName','Times New Roman','FontSize',15);
% grid on
% grid minor

figure(5);
semilogy(M,error_TPCA,'-ob','LineWidth',2,...
    'MarkerSize',5,...
    'MarkerFaceColor','b');
xlabel('Iteration Number  ','FontName','Times New Rom an','FontSize',30);
ylabel('Error','FontName','Times New Roman','FontSize',30);
% legend('\it N_P','\it N_N','\it OF','\it N_E','FontName','Times New Roman','Location','North','Orientation','horizontal');
set(gca,'xtick',MM,'xticklabel',MM) 
set(gca,'FontName','Times New Roman','FontSize',20);
grid on
grid minor