
% Get the 1D cuts to see how the terms releates each other (keep just the energy terms)
%--------------------------------------------------------------------------
%
%----------------------------------------------------------------------
ypc_i=127; % <---- Here is where I choose the 1D line that crosses the reconnection event 
xpc_i = 93; xll(xpc_i)
% For electrons
%--------------------------------------------------------------------------
cut_1d=ve_Gp{1,1}'; ve_1_1d=cut_1d(xpc_i,:); 
cut_1d=ve_Gp{1,2}'; ve_2_1d=cut_1d(xpc_i,:); 
cut_1d=ve_Gp{1,3}'; ve_3_1d=cut_1d(xpc_i,:); 

% Time derivatives
cut_1d=dtKe_in_Gp_2{1,1}'; dtKe_in_Gp_1d = cut_1d(xpc_i,:);
cut_1d=dtUe_in_Gp_2{1,1}'; dtUe_in_Gp_1d = cut_1d(xpc_i,:);

% Kinetic
%cut_1d=partialKe_Gp{1,1}'; %partialKe_1d = cut_1d(xpc_i,:);
cut_1d=ve_grad_Ke_Gp{1,1}'; ve_grad_Ke_1d = cut_1d(xpc_i,:);
cut_1d=ve_Div_Pet_Gp{1,1}'; ve_Div_Pet_1d = cut_1d(xpc_i,:);
cut_1d=Ke_Div_ve_Gp{1,1}'; Ke_Div_ve_1d = cut_1d(xpc_i,:);
cut_1d=qe_ne_ve_E_Gp{1,1}'; qe_ne_ve_E_Gp_1d = cut_1d(xpc_i,:);

% Thermal
%cut_1d=partialUe_Gp{1,1}'; partialUe_1d = cut_1d(xpc_i,:);
cut_1d=ve_grad_Ue_Gp{1,1}'; ve_grad_Ue_1d = cut_1d(xpc_i,:);
cut_1d=dkQiik_05_e_Gp{1,1}'; dkQiik_05_e_Gp_1d = cut_1d(xpc_i,:);
cut_1d=Pet_grad_ve_Gp{1,1}'; Pet_grad_ve_Gp_1d = cut_1d(xpc_i,:);
cut_1d=Pet_grad_ve_ii_Gp{1,1}'; Pet_grad_ve_ii_Gp_1d = cut_1d(xpc_i,:);
cut_1d=Pet_grad_ve_ij_Gp{1,1}'; Pet_grad_ve_ij_Gp_1d = cut_1d(xpc_i,:);
cut_1d=Ue_Div_ve_Gp{1,1}'; Ue_Div_ve_Gp_1d = cut_1d(xpc_i,:);
cut_1d=Div_Pes_dot_ve_Gp{1,1}'; Div_Pes_dot_ve_Gp_1d = cut_1d(xpc_i,:);

%Total derivatives
total_dKedt_1d = dtKe_in_Gp_1d + ve_grad_Ke_1d;
total_dUedt_1d = dtUe_in_Gp_1d + ve_grad_Ue_1d;

% pressure temrs
cut_1d=(Pij_e_Gp{1,1}' + Pij_e_Gp{2,2}' + Pij_e_Gp{3,3}')/3; pe_e_Gp_1d = cut_1d(xpc_i,:);
cut_1d=(Pij_e_Gp{1,2}' + Pij_e_Gp{2,1}')/2; Pi_era_Gp_1d = cut_1d(xpc_i,:);
cut_1d=(Pij_e_Gp{3,2}' + Pij_e_Gp{2,3}')/2; Pi_epa_Gp_1d = cut_1d(xpc_i,:);
cut_1d=(Pij_e_Gp{3,1}' + Pij_e_Gp{1,3}')/2; Pi_epr_Gp_1d = cut_1d(xpc_i,:);

% Parameter proxies
cut_1d= Zernitani_e_Gp{1,1}'; Zernitani_e_Gp_1d = cut_1d(xpc_i,:);
cut_1d= p_theta_e_Gp{1,1}'; p_theta_e_Gp_1d = cut_1d(xpc_i,:);
cut_1d= PiD_e_Gp{1,1}'; PiD_e_Gp_1d = cut_1d(xpc_i,:);

% collisionale terms
cut_1d=dtKe_in_Gp_2{1,1}' + ve_grad_Ke_Gp{1,1}' + ve_Div_Pet_Gp{1,1}' + ...
    Ke_Div_ve_Gp{1,1}' -qe_ne_ve_E_Gp{1,1}'; dum_pXIe_Gp_1d = cut_1d(xpc_i,:);

cut_1d=dtUe_in_Gp_2{1,1}' + ve_grad_Ue_Gp{1,1}' + dkQiik_05_e_Gp{1,1}' + ...
    Pet_grad_ve_Gp{1,1}' + Ue_Div_ve_Gp{1,1}' +...
    (dtKe_in_Gp_2{1,1}' + ve_grad_Ke_Gp{1,1}' + ve_Div_Pet_Gp{1,1}' + ...
    Ke_Div_ve_Gp{1,1}' -qe_ne_ve_E_Gp{1,1}'); dum_pXI2e_Gp_1d=cut_1d(xpc_i,:);

% Density
cut_1d=ne_Gp{1,1}'; ne_Gp_1d=cut_1d(xpc_i,:);

%--------------------------------------------------------------------------

% For Ions
%--------------------------------------------------------------------------
cut_1d=vi_Gp{1,1}'; vi_1_1d=cut_1d(xpc_i,:); 
cut_1d=vi_Gp{1,2}'; vi_2_1d=cut_1d(xpc_i,:); 
cut_1d=vi_Gp{1,3}'; vi_3_1d=cut_1d(xpc_i,:); 

% Time derivatives
cut_1d=dtKi_in_Gp_2{1,1}'; dtKi_in_Gp_1d = cut_1d(xpc_i,:);
cut_1d=dtUi_in_Gp_2{1,1}'; dtUi_in_Gp_1d = cut_1d(xpc_i,:);

% Kinetic
%cut_1d=partialKi_Gp{1,1}'; %partialKi_1d = cut_1d(xpc_i,:);
cut_1d=vi_grad_Ki_Gp{1,1}'; vi_grad_Ki_1d = cut_1d(xpc_i,:);
cut_1d=vi_Div_Pit_Gp{1,1}'; vi_Div_Pit_1d = cut_1d(xpc_i,:);
cut_1d=Ki_Div_vi_Gp{1,1}'; Ki_Div_vi_1d = cut_1d(xpc_i,:);
cut_1d=qi_ni_vi_E_Gp{1,1}'; qi_ni_vi_E_Gp_1d = cut_1d(xpc_i,:);

% Thermal
%cut_1d=partialUe_Gp{1,1}'; partialUe_1d = cut_1d(xpc_i,:);
cut_1d=vi_grad_Ui_Gp{1,1}'; vi_grad_Ui_1d = cut_1d(xpc_i,:);
cut_1d=dkQiik_05_i_Gp{1,1}'; dkQiik_05_i_Gp_1d = cut_1d(xpc_i,:);
cut_1d=Pit_grad_vi_Gp{1,1}'; Pit_grad_vi_Gp_1d = cut_1d(xpc_i,:);
cut_1d=Pit_grad_vi_ii_Gp{1,1}'; Pit_grad_vi_ii_Gp_1d = cut_1d(xpc_i,:);
cut_1d=Pit_grad_vi_ij_Gp{1,1}'; Pit_grad_vi_ij_Gp_1d = cut_1d(xpc_i,:);
cut_1d=Ui_Div_vi_Gp{1,1}'; Ui_Div_vi_Gp_1d = cut_1d(xpc_i,:);
cut_1d=Div_Pis_dot_vi_Gp{1,1}'; Div_Pis_dot_vi_Gp_1d = cut_1d(xpc_i,:);

%Total derivatives
total_dKidt_1d = dtKi_in_Gp_1d + vi_grad_Ki_1d;
total_dUidt_1d = dtUi_in_Gp_1d + vi_grad_Ui_1d;

% pressure temrs
cut_1d=(Pij_i_Gp{1,1}' + Pij_i_Gp{2,2}' + Pij_i_Gp{3,3}')/3; pe_i_Gp_1d = cut_1d(xpc_i,:);
cut_1d=(Pij_i_Gp{1,2}' + Pij_i_Gp{2,1}')/2; Pi_ira_Gp_1d = cut_1d(xpc_i,:);
cut_1d=(Pij_i_Gp{3,2}' + Pij_i_Gp{2,3}')/2; Pi_ipa_Gp_1d = cut_1d(xpc_i,:);
cut_1d=(Pij_i_Gp{3,1}' + Pij_i_Gp{1,3}')/2; Pi_ipr_Gp_1d = cut_1d(xpc_i,:);

% Parameter proxies
cut_1d= Zernitani_i_Gp{1,1}'; Zernitani_i_Gp_1d = cut_1d(xpc_i,:);
cut_1d= p_theta_i_Gp{1,1}'; p_theta_i_Gp_1d = cut_1d(xpc_i,:);
cut_1d= PiD_i_Gp{1,1}'; PiD_i_Gp_1d = cut_1d(xpc_i,:);

% collisionale terms
cut_1d=dtKi_in_Gp_2{1,1}' + vi_grad_Ki_Gp{1,1}' + vi_Div_Pit_Gp{1,1}' + ...
    Ki_Div_vi_Gp{1,1}' -qi_ni_vi_E_Gp{1,1}'; dum_pXIi_Gp_1d = cut_1d(xpc_i,:);

cut_1d=dtUi_in_Gp_2{1,1}' + vi_grad_Ui_Gp{1,1}' + dkQiik_05_i_Gp{1,1}' + ...
    Pit_grad_vi_Gp{1,1}' + Ui_Div_vi_Gp{1,1}' +...
    (dtKi_in_Gp_2{1,1}' + vi_grad_Ki_Gp{1,1}' + vi_Div_Pit_Gp{1,1}' + ...
    Ki_Div_vi_Gp{1,1}' -qi_ni_vi_E_Gp{1,1}'); dum_pXI2i_Gp_1d=cut_1d(xpc_i,:);

% Density
cut_1d=ni_Gp{1,1}'; ni_Gp_1d=cut_1d(xpc_i,:);
%--------------------------------------------------------------------------

% Magnetic field
cut_1d=B_Gp{1,1}'; B_1_1d=cut_1d(xpc_i,:); 
cut_1d=B_Gp{1,2}'; B_2_1d=cut_1d(xpc_i,:); 
cut_1d=B_Gp{1,3}'; B_3_1d=cut_1d(xpc_i,:); 
%--------------------------------------------------------------------------

%}
%--------------------------------------------------------------------------

%  <--- internal {} (not internal energy) one
%{
f213 = figure(213);
%subplot(1,2,1)
%plot(xll,partialKe_1d + ve_grad_Ke_1d,'k')
%hold on
plot(xll,ve_Div_Pet_1d,'-*k')
hold on
plot(xll,Ke_Div_ve_1d,'b')
plot(xll,-qe_ne_ve_E_Gp_1d,'m')
xline(4,'--k'); xline(8,'--k'); xline(5.7,'--b')
%legend('dKe/dt','vedivPe','Kedivve','-qEe')
legend('u_{e} \cdot \nabla \cdot P_{e}','K_{e}\nabla \cdot u_{e}','-qE\cdot v_{e}')
%plot(xll,ve_2_1d)
xlabel('r / d_{i}')
ylabel('Energy rate')
hold off
%subplot(1,2,2)
%dum_p=ve_Gp{1,2};
%hc = pcolor(XLL,YLL,dum_p);
%set(hc,'edgecolor','none')
%}

fs=18; lw=1.5;

% 1D plots for electrons no smoothing 
%--------------------------------------------------------------------------
%{
f214 = figure(214);
plot(xll,ve_1_1d,'Color',colorblind(1,:), 'LineWidth', lw)
hold on
plot(xll,ve_2_1d,'Color',colorblind(2,:), 'LineWidth', lw)
plot(xll,ve_3_1d,'Color',colorblind(5,:), 'LineWidth', lw)
xline(4,'--k'); xline(8,'--k'); xline(5.7,'--b')
ax = gca;
ax.FontSize = fs;
ax.TickLabelInterpreter ='latex';
legend('$u_{p}$','$u_{a}$','$u_{r}$','Interpreter','latex','FontSize',fs)
xlabel('$r / d_{i}$','Interpreter','latex','FontSize',fs)
ylabel('$u_{e}$ components','Interpreter','latex','FontSize',fs)
xlim([4,8])
hold off

f2145 = figure(2145);
plot(xll,ve_1_1d,'Color',colorblind(1,:), 'LineWidth', lw)
hold on
plot(xll,ve_2_1d,'Color',colorblind(2,:), 'LineWidth', lw)
plot(xll,ve_3_1d,'Color',colorblind(5,:), 'LineWidth', lw)
xline(4,'--k'); xline(8,'--k'); xline(5.7,'--b')
ax = gca;
ax.FontSize = fs;
ax.TickLabelInterpreter ='latex';
legend('$u_{p}$','$u_{a}$','$u_{r}$','Interpreter','latex','FontSize',fs)
xlabel('$r / d_{i}$','Interpreter','latex','FontSize',fs)
ylabel('$u_{e}$ components','Interpreter','latex','FontSize',fs)
xlim([4,8])
hold off

f215 = figure(215);
% Kinetic subplot(1,2,1)
h1=subaxis(1,2,1,'SV',0,'SH',0.05,'MR',0.05,'ML',0.05,'PL',0.05,'PR',0.01);
plot(xll,ve_grad_Ke_1d,'Color',colorblind(1,:), 'LineWidth', lw)
hold on
plot(xll,ve_Div_Pet_1d,'Color',colorblind(2,:), 'LineWidth', lw)
plot(xll,Ke_Div_ve_1d,'Color',colorblind(5,:), 'LineWidth', lw)
plot(xll,-qe_ne_ve_E_Gp_1d,'Color',colorblind(6,:), 'LineWidth', lw)
xline(4,'--k'); xline(8,'--k'); xline(5.7,'--b')
ax = gca; ax.FontSize = fs; ax.TickLabelInterpreter ='latex';
legend('$u_{e} \cdot \nabla \varepsilon^{k}_{e}$','$u_{e} \cdot \nabla \cdot \overline{P}_{e}$',...
    '$\varepsilon^{k}_{e}\nabla \cdot u_{e}$','$-qE\cdot u_{e}$','Interpreter','latex','FontSize',fs,...
    'Location','southeast','NumColumns',2)
xlabel('$r / d_{i}$','Interpreter','latex','FontSize',fs)
%xlim([4,8])
ylabel('kinetic Energy rate terms','Interpreter','latex','FontSize',fs)
hold off
% Thermal
%------------------------------------
%subplot(1,2,2)
h1=subaxis(1,2,2,'SV',0,'SH',0.05,'MR',0.05,'ML',0.05,'PL',0.05,'PR',0.01);
plot(xll,ve_grad_Ue_1d,'Color',colorblind(1,:), 'LineWidth', lw)
hold on
plot(xll,dkQiik_05_e_Gp_1d,'Color',colorblind(2,:), 'LineWidth', lw)
plot(xll,Pet_grad_ve_Gp_1d,'Color',colorblind(5,:), 'LineWidth', lw)
plot(xll,Pet_grad_ve_ii_Gp_1d,'Color',colorblind(6,:), 'LineWidth', lw)
plot(xll,Pet_grad_ve_ij_Gp_1d,'Color',colorblind(8,:), 'LineWidth', lw)
plot(xll,Ue_Div_ve_Gp_1d,'Color',colorblind(9,:), 'LineWidth', lw)
xline(4,'--k'); xline(8,'--k'); xline(5.7,'--b')
ax = gca; ax.FontSize = fs; ax.TickLabelInterpreter ='latex';
legend('$u_{e} \cdot \nabla \varepsilon^{th}_{e}$','$\nabla \cdot h_{e}$',...
    '$\nabla u_{e}:\overline{P}_{e}$','$p\theta_{e}$','$PiD_{e}$',...
    '$\varepsilon^{th}_{e} \nabla \cdot u_{e}$','Interpreter','latex','FontSize',fs,...
    'Location','southeast','NumColumns',2)
xlabel('$r / d_{i}$','Interpreter','latex','FontSize',fs)
%xlim([4,8])
ylabel('Thermal Energy rate terms','Interpreter','latex','FontSize',fs)
hold off

f2155 = figure(2155);
% Kinetic
%subplot(1,2,1)
h1=subaxis(1,2,1,'SV',0,'SH',0.05,'MR',0.05,'ML',0.05,'PL',0.05,'PR',0.01);
plot(xll,ve_grad_Ke_1d,'Color',colorblind(1,:), 'LineWidth', lw)
hold on
plot(xll,ve_Div_Pet_1d,'Color',colorblind(2,:), 'LineWidth', lw)
plot(xll,Ke_Div_ve_1d,'Color',colorblind(5,:), 'LineWidth', lw)
plot(xll,-qe_ne_ve_E_Gp_1d,'Color',colorblind(6,:), 'LineWidth', lw)
xline(4,'--k'); xline(8,'--k'); xline(5.7,'--b')
ax = gca; ax.FontSize = fs; ax.TickLabelInterpreter ='latex';
legend('$u_{e} \cdot \nabla \varepsilon^{k}_{e}$','$u_{e} \cdot \nabla \cdot \overline{P}_{e}$',...
    '$\varepsilon^{k}_{e}\nabla \cdot u_{e}$','$-qE\cdot u_{e}$','Interpreter','latex','FontSize',fs,...
    'Location','southeast','NumColumns',2)
xlabel('$r / d_{i}$','Interpreter','latex','FontSize',fs)
xlim([4,8])
ylabel('kinetic Energy rate terms','Interpreter','latex','FontSize',fs)
hold off
% Thermal
%------------------------------------
%subplot(1,2,2)
h1=subaxis(1,2,2,'SV',0,'SH',0.05,'MR',0.05,'ML',0.05,'PL',0.05,'PR',0.01);
plot(xll,ve_grad_Ue_1d,'Color',colorblind(1,:), 'LineWidth', lw)
hold on
plot(xll,dkQiik_05_e_Gp_1d,'Color',colorblind(2,:), 'LineWidth', lw)
plot(xll,Pet_grad_ve_Gp_1d,'Color',colorblind(5,:), 'LineWidth', lw)
plot(xll,Pet_grad_ve_ii_Gp_1d,'Color',colorblind(6,:), 'LineWidth', lw)
plot(xll,Pet_grad_ve_ij_Gp_1d,'Color',colorblind(8,:), 'LineWidth', lw)
plot(xll,Ue_Div_ve_Gp_1d,'Color',colorblind(9,:), 'LineWidth', lw)
xline(4,'--k'); xline(8,'--k'); xline(5.7,'--b')
ax = gca; ax.FontSize = fs; ax.TickLabelInterpreter ='latex';
legend('$u_{e} \cdot \nabla \varepsilon^{th}_{e}$','$\nabla \cdot h_{e}$',...
    '$\nabla u_{e}:\overline{P}_{e}$','$p\theta_{e}$','$PiD_{e}$',...
    '$\varepsilon^{th}_{e} \nabla \cdot u_{e}$','Interpreter','latex','FontSize',fs,...
    'Location','southeast','NumColumns',2)
xlabel('$r / d_{i}$','Interpreter','latex','FontSize',fs)
xlim([4,8])
ylabel('Thermal Energy rate terms','Interpreter','latex','FontSize',fs)
hold off

f216 = figure(216);
plot(xll,Div_Pes_dot_ve_Gp_1d,'Color',colorblind(1,:), 'LineWidth', lw)
hold on
plot(xll,Pet_grad_ve_Gp_1d,'Color',colorblind(2,:), 'LineWidth', lw)
plot(xll,ve_Div_Pet_1d,'Color',colorblind(5,:), 'LineWidth', lw)
%plot(xll,Pet_grad_ve_Gp_1d + ve_Div_Pet_1d,'g')
xline(4,'--k'); xline(8,'--k'); xline(5.7,'--b')
ax = gca; ax.FontSize = fs; ax.TickLabelInterpreter ='latex';
%legend('dKe/dt','vedivPe','Kedivve','-qEe')
legend('$\nabla \cdot( \overline{P}_{e} \cdot u_{e})$ ','$\overline{P}_{e}:\nabla u_{e}$',...
    '$u_{e} \cdot(\nabla \cdot \overline{P}_{e})$','Interpreter','latex','FontSize',fs)
%plot(xll,ve_2_1d)
xlabel('$r / d_{i}$','Interpreter','latex','FontSize',fs)
ylabel('Kinetic Energy rate terms','Interpreter','latex','FontSize',fs)
%xlim([4,8])
hold off

f2165 = figure(2165);
plot(xll,Div_Pes_dot_ve_Gp_1d,'Color',colorblind(1,:), 'LineWidth', lw)
hold on
plot(xll,Pet_grad_ve_Gp_1d,'Color',colorblind(2,:), 'LineWidth', lw)
plot(xll,ve_Div_Pet_1d,'Color',colorblind(5,:), 'LineWidth', lw)
%plot(xll,Pet_grad_ve_Gp_1d + ve_Div_Pet_1d,'g')
xline(4,'--k'); xline(8,'--k'); xline(5.7,'--b')
ax = gca; ax.FontSize = fs; ax.TickLabelInterpreter ='latex';
%legend('dKe/dt','vedivPe','Kedivve','-qEe')
legend('$\nabla \cdot( \overline{P}_{e} \cdot u_{e})$ ','$\overline{P}_{e}:\nabla u_{e}$',...
    '$u_{e} \cdot(\nabla \cdot \overline{P}_{e})$','Interpreter','latex','FontSize',fs)
%plot(xll,ve_2_1d)
xlabel('$r / d_{i}$','Interpreter','latex','FontSize',fs)
ylabel('Kinetic Energy rate terms','Interpreter','latex','FontSize',fs)
xlim([4,8])
hold off

f217 = figure(217);
plot(xll,Pet_grad_ve_Gp_1d,'Color',colorblind(1,:), 'LineWidth', lw)
hold on
plot(xll,Pet_grad_ve_ii_Gp_1d,'Color',colorblind(2,:), 'LineWidth', lw)
plot(xll,Pet_grad_ve_ij_Gp_1d,'Color',colorblind(5,:), 'LineWidth', lw)
%plot(xll,Pet_grad_ve_Gp_1d + ve_Div_Pet_1d,'g')
xline(4,'--k'); xline(8,'--k'); xline(5.7,'--b')
ax = gca; ax.FontSize = fs; ax.TickLabelInterpreter ='latex';
%legend('dKe/dt','vedivPe','Kedivve','-qEe')
legend('$\overline{P}_{e}:\nabla u_{e}$','$p\theta_{e}$','$PiD_{e}$','Interpreter','latex','FontSize',fs)
%plot(xll,ve_2_1d)
xlabel('$r / d_{i}$','Interpreter','latex','FontSize',fs)
ylabel('Energy rate terms','Interpreter','latex','FontSize',fs)
%xlim([4,8])
hold off

f2175 = figure(2175);
plot(xll,Pet_grad_ve_Gp_1d,'Color',colorblind(1,:), 'LineWidth', lw)
hold on
plot(xll,Pet_grad_ve_ii_Gp_1d,'Color',colorblind(2,:), 'LineWidth', lw)
plot(xll,Pet_grad_ve_ij_Gp_1d,'Color',colorblind(5,:), 'LineWidth', lw)
%plot(xll,Pet_grad_ve_Gp_1d + ve_Div_Pet_1d,'g')
xline(4,'--k'); xline(8,'--k'); xline(5.7,'--b')
ax = gca; ax.FontSize = fs; ax.TickLabelInterpreter ='latex';
%legend('dKe/dt','vedivPe','Kedivve','-qEe')
legend('$\overline{P}_{e}:\nabla u_{e}$','$p\theta_{e}$','$PiD_{e}$','Interpreter','latex','FontSize',fs)
%plot(xll,ve_2_1d)
xlabel('$r / d_{i}$','Interpreter','latex','FontSize',fs)
ylabel('Energy rate terms','Interpreter','latex','FontSize',fs)
xlim([4,8])
hold off

f218 = figure(218);
%dum1=(Pet_grad_ve_Gp_1d-Pet_grad_ve_ii_Gp_1d)./Pet_grad_ve_Gp_1d;
dum1=(Pet_grad_ve_ij_Gp_1d)./Pet_grad_ve_Gp_1d;
semilogy(xll,10+dum1,'Color',colorblind(1,:), 'LineWidth', lw)
hold on
dum2=(Pet_grad_ve_ii_Gp_1d)./Pet_grad_ve_Gp_1d;
semilogy(xll,10+dum2,'Color',colorblind(2,:), 'LineWidth', lw)
xline(4,'--k'); xline(8,'--k'); xline(5.7,'--b')
%legend('dKe/dt','vedivPe','Kedivve','-qEe')
ax = gca; ax.FontSize = fs; ax.TickLabelInterpreter ='latex';
legend('$10+\frac{PiD_{e}}{P_{e}:\nabla u_{e}}$',...
    '$10+\frac{p\theta_{e}}{P_{e}:\nabla u_{e}}$',...
    'Interpreter','latex','FontSize',fs)
%plot(xll,ve_2_1d)
xlabel('$r / d_{i}$','Interpreter','latex','FontSize',fs)
ylabel('Energy rate terms','Interpreter','latex','FontSize',fs)
%xlim([4,8])
hold off

f2185 = figure(2185);
%dum1=(Pet_grad_ve_Gp_1d-Pet_grad_ve_ii_Gp_1d)./Pet_grad_ve_Gp_1d;
dum1=(Pet_grad_ve_ij_Gp_1d)./Pet_grad_ve_Gp_1d;
semilogy(xll,10+dum1,'Color',colorblind(1,:), 'LineWidth', lw)
hold on
dum2=(Pet_grad_ve_ii_Gp_1d)./Pet_grad_ve_Gp_1d;
semilogy(xll,10+dum2,'Color',colorblind(2,:), 'LineWidth', lw)
xline(4,'--k'); xline(8,'--k'); xline(5.7,'--b')
%legend('dKe/dt','vedivPe','Kedivve','-qEe')
ax = gca; ax.FontSize = fs; ax.TickLabelInterpreter ='latex';
legend('$10+\frac{PiD_{e}}{P_{e}:\nabla u_{e}}$',...
    '$10+\frac{p\theta_{e}}{P_{e}:\nabla u_{e}}$',...
    'Interpreter','latex','FontSize',fs)
%plot(xll,ve_2_1d)
xlabel('$r / d_{i}$','Interpreter','latex','FontSize',fs)
ylabel('Energy rate terms','Interpreter','latex','FontSize',fs)
xlim([4,8])
hold off
%}
%--------------------------------------------------------------------------



% 1D Plots for electros smoothed   <------------------------------------
%--------------------------------------------------------------------------
%
cd '/Users/jeffersson_agudelo/Documents/CB104_local_data/Matlab_Scripts_here'

% Components of u electrons
%--------------------------------------------------------------------------
%
f2142 = figure(2142);
plot(yll,smoothdata(ve_1_1d)./VA2c,'Color',colorblind(1,:), 'LineWidth', lw)
hold on
plot(yll,smoothdata(ve_2_1d)./VA2c,'Color',colorblind(2,:), 'LineWidth', lw)
plot(yll,smoothdata(ve_3_1d)./VA2c,'Color',colorblind(5,:), 'LineWidth', lw)
xline(4,'--k'); xline(8,'--k'); xline(5.7,'--b')
ax = gca;
ax.FontSize = fs;
ax.TickLabelInterpreter ='latex';
legend('$u_{p}/V_{A,i}$','$u_{a}/V_{A,i}$','$u_{r}/V_{A,i}$','Interpreter','latex','FontSize',fs)
xlabel('$r / d_{i}$','Interpreter','latex','FontSize',fs)
ylabel('$\mathbf{u}_{e}$ components','Interpreter','latex','FontSize',fs)
xlim([0,10])
ylim([-1.9, 1.9])
sgtitle('$$t = \ $$'+string(round(time_t,2)) + '$$ \ \omega_{i}^{-1}$$','Interpreter','latex','FontSize',20)
hold off
%{
f2145 = figure(2145);
plot(xll,smoothdata(ve_1_1d)./VA2c,'Color',colorblind(1,:), 'LineWidth', lw)
hold on
plot(xll,smoothdata(ve_2_1d)./VA2c,'Color',colorblind(2,:), 'LineWidth', lw)
plot(xll,smoothdata(ve_3_1d)./VA2c,'Color',colorblind(5,:), 'LineWidth', lw)
xline(4,'--k'); xline(8,'--k'); xline(5.7,'--b')
ax = gca;
ax.FontSize = fs;
ax.TickLabelInterpreter ='latex';
legend('$u_{p}/V_{A,i}$','$u_{a}/V_{A,i}$','$u_{r}/V_{A,i}$','Interpreter','latex','FontSize',fs)
xlabel('$r / d_{i}$','Interpreter','latex','FontSize',fs)
ylabel('$\mathbf{u}_{e}$ components','Interpreter','latex','FontSize',fs)
xlim([4,8])
ylim([-1.9, 1.9])
sgtitle('$$t = \ $$'+string(round(time_t,2)) + '$$ \ \omega_{i}^{-1}$$','Interpreter','latex','FontSize',20)
hold off
%}
%--------------------------------------------------------------------------

% Components of u ions
%--------------------------------------------------------------------------
%{
f2146 = figure(2146);
plot(xll,smoothdata(vi_1_1d)./VA2c,'Color',colorblind(1,:), 'LineWidth', lw)
%plot(xll,vi_1_1d./VA2c,'Color',colorblind(1,:), 'LineWidth', lw)
hold on
plot(xll,smoothdata(vi_2_1d)./VA2c,'Color',colorblind(2,:), 'LineWidth', lw)
plot(xll,smoothdata(vi_3_1d)./VA2c,'Color',colorblind(5,:), 'LineWidth', lw)
xline(4,'--k'); xline(8,'--k'); xline(5.7,'--b')
ax = gca;
ax.FontSize = fs;
ax.TickLabelInterpreter ='latex';
legend('$u_{p}/V_{A,i}$','$u_{a}/V_{A,i}$','$u_{r}/V_{A,i}$','Interpreter','latex','FontSize',fs)
xlabel('$r / d_{i}$','Interpreter','latex','FontSize',fs)
ylabel('$\mathbf{u}_{i}$ components','Interpreter','latex','FontSize',fs)
%xlim([4,8])
ylim([-1.9, 1.9])
sgtitle('$$t = \ $$'+string(round(time_t,2)) + '$$ \ \omega_{i}^{-1}$$','Interpreter','latex','FontSize',20)
hold off

f21456 = figure(21456);
plot(xll,smoothdata(vi_1_1d)./VA2c,'Color',colorblind(1,:), 'LineWidth', lw)
hold on
plot(xll,smoothdata(vi_2_1d)./VA2c,'Color',colorblind(2,:), 'LineWidth', lw)
plot(xll,smoothdata(vi_3_1d)./VA2c,'Color',colorblind(5,:), 'LineWidth', lw)
xline(4,'--k'); xline(8,'--k'); xline(5.7,'--b')
ax = gca;
ax.FontSize = fs;
ax.TickLabelInterpreter ='latex';
legend('$u_{p}/V_{A,i}$','$u_{a}/V_{A,i}$','$u_{r}/V_{A,i}$','Interpreter','latex','FontSize',fs)
xlabel('$r / d_{i}$','Interpreter','latex','FontSize',fs)
ylabel('$\mathbf{u}_{i}$ components','Interpreter','latex','FontSize',fs)
xlim([4,8])
ylim([-1.9, 1.9])
sgtitle('$$t = \ $$'+string(round(time_t,2)) + '$$ \ \omega_{i}^{-1}$$','Interpreter','latex','FontSize',20)
hold off
%}
%--------------------------------------------------------------------------

% density fluctuations, magnetic and pressure
%--------------------------------------------------------------------------
%{
f2148 = figure(2148);
plot(xll,smoothdata(sqrt(ve_1_1d.^2 + ve_2_1d.^2 + ve_3_1d.^2))./VA2c,'Color',colorblind(1,:), 'LineWidth', lw)
hold on
plot(xll,smoothdata(sqrt(vi_1_1d.^2 + vi_2_1d.^2 + vi_3_1d.^2))./VA2c,'Color',colorblind(6,:), 'LineWidth', lw)
plot(xll,smoothdata(ne_Gp_1d),'Color',colorblind(2,:), 'LineWidth', lw)
plot(xll,smoothdata(ni_Gp_1d),'Color',colorblind(5,:), 'LineWidth', lw)
plot(xll,smoothdata(pe_e_Gp_1d)./Pi0,'Color',colorblind(3,:), 'LineWidth', lw)
plot(xll,smoothdata(pe_i_Gp_1d)./Pi0,'Color',colorblind(7,:), 'LineWidth', lw)
plot(xll,smoothdata(sqrt(B_1_1d.^2 + B_2_1d.^2 + B_3_1d.^2))./B0,'Color',colorblind(4,:), 'LineWidth', lw)
plot(xll,smoothdata(pe_e_Gp_1d./ne_Gp_1d)./Pi0,'Color',colorblind(8,:), 'LineWidth', lw)
plot(xll,smoothdata(pe_i_Gp_1d./ni_Gp_1d)./Pi0,'Color',colorblind(9,:), 'LineWidth', lw)
xline(4,'--k'); xline(8,'--k'); xline(5.7,'--b')
ax = gca;
ax.FontSize = fs;
ax.TickLabelInterpreter ='latex';
legend('$|u_{e}|/V_{A,i}$',...
    '$|u_{i}|/V_{A,i}$',...
    '$n_{e}/n_{e0}$',...
    '$n_{i}/n_{i0}$',...
    '$p_{e}/P_{i0}$',...
    '$p_{i}/P_{i0}$',...
    '$|B|/B_{0}$',...
    '$p_{e}/n_{e}$',...
    '$p_{i}/n_{i}$',...
    'Interpreter','latex','FontSize',fs,'Location','southeast','NumColumns',2)
xlabel('$r / d_{i}$','Interpreter','latex','FontSize',fs)
ylabel('Several variables','Interpreter','latex','FontSize',fs)
%xlim([4,8])
ylim([-2.9, 2.9])
sgtitle('$$t = \ $$'+string(round(time_t,2)) + '$$ \ \omega_{i}^{-1}$$','Interpreter','latex','FontSize',20)
hold off

%}
%--------------------------------------------------------------------------


% Kinetic and thermal energy terms for electrons
%--------------------------------------------------------------------------
%{
f215 = figure(215);
sgtitle('$$t = \ $$'+string(round(time_t,2)) + '$$ \ \omega_{i}^{-1}$$','Interpreter','latex','FontSize',20)
% Kinetic subplot(1,2,1)
h1=subaxis(1,2,1,'SV',0,'SH',0.05,'MR',0.05,'ML',0.05,'PL',0.05,'PR',0.01);
plot(xll,smoothdata(total_dKedt_1d)./(VA2c*Pi0),'Color',colorblind(10,:), 'LineWidth', lw)
hold on
plot(xll,smoothdata(ve_Div_Pet_1d)./(VA2c*Pi0),'Color',colorblind(2,:), 'LineWidth', lw)
plot(xll,smoothdata(Ke_Div_ve_1d)./(VA2c*Pi0),'Color',colorblind(5,:), 'LineWidth', lw)
plot(xll,smoothdata(-qe_ne_ve_E_Gp_1d)./(VA2c*Pi0),'Color',colorblind(6,:), 'LineWidth', lw)
plot(xll,smoothdata(dum_pXIe_Gp_1d)./(VA2c*Pi0),'Color',colorblind(1,:), 'LineWidth', lw)
xline(4,'--k'); xline(8,'--k'); xline(5.7,'--b')
ax = gca; ax.FontSize = fs; ax.TickLabelInterpreter ='latex';
legend('$\frac{d \varepsilon^{k}_{e} /d t}{\Delta\varepsilon_{0}} $',...
    '$\frac{u_{e} \cdot \nabla \cdot \overline{P}_{e}}{\Delta\varepsilon_{0}}$',...
    '$\frac{\varepsilon^{k}_{e}\nabla \cdot u_{e}}{\Delta\varepsilon_{0}}$',...
    '$\frac{-q_{e}n_{e}E\cdot u_{e}}{\Delta\varepsilon_{0}}$',...
    '$\frac{m_{e}\mathbf{u}_{e} \cdot \mathbf{\Xi}_{e}}{\Delta\varepsilon_{0}}$',...
    'Interpreter','latex','FontSize',fs,...
    'Location','southeast','NumColumns',2)
xlabel('$r / d_{i}$','Interpreter','latex','FontSize',fs)
%xlim([4,8])
%ylim([-0.005, 0.005])
ylim([-3.2, 3.2]) %<------------
ylabel('kinetic Energy rate terms','Interpreter','latex','FontSize',fs)
hold off
% Thermal
%------------------------------------
%subplot(1,2,2)
h1=subaxis(1,2,2,'SV',0,'SH',0.05,'MR',0.05,'ML',0.05,'PL',0.05,'PR',0.01);
plot(xll,smoothdata(total_dUedt_1d),'Color',colorblind(10,:), 'LineWidth', lw)
hold on
plot(xll,smoothdata(dkQiik_05_e_Gp_1d)./(VA2c*Pi0),'Color',colorblind(2,:), 'LineWidth', lw)
plot(xll,smoothdata(Pet_grad_ve_Gp_1d)./(VA2c*Pi0),'Color',colorblind(5,:), 'LineWidth', lw)
plot(xll,smoothdata(Ue_Div_ve_Gp_1d)./(VA2c*Pi0),'Color',colorblind(6,:), 'LineWidth', lw)
plot(xll,smoothdata(dum_pXI2e_Gp_1d)./(VA2c*Pi0),'Color',colorblind(1,:), 'LineWidth', lw)
xline(4,'--k'); xline(8,'--k'); xline(5.7,'--b')
ax = gca; ax.FontSize = fs; ax.TickLabelInterpreter ='latex';
legend('$\frac{d \varepsilon^{th}_{e} / d t }{\Delta\varepsilon_{0}}$',...
    '$\frac{\nabla \cdot h_{e}}{\Delta\varepsilon_{0}}$',...
    '$\frac{\nabla u_{e}:\overline{P}_{e}}{\Delta\varepsilon_{0}}$',...
    '$\frac{\varepsilon^{th}_{e} \nabla \cdot u_{e}}{\Delta\varepsilon_{0}}$',...
    '$\frac{0.5 Tr(m_{e}\overline{\mathbf{\Xi}}_{e})}{\Delta\varepsilon_{0}}$',...
    'Interpreter','latex','FontSize',fs,...
    'Location','southeast','NumColumns',2)
xlabel('$r / d_{i}$','Interpreter','latex','FontSize',fs)
%xlim([4,8])
ylim([-10.2, 10.2]) %<--------------
%ylim([-0.015, 0.015])
ylabel('Thermal Energy rate terms','Interpreter','latex','FontSize',fs)
hold off

f2155 = figure(2155);
sgtitle('$$t = \ $$'+string(round(time_t,2)) + '$$ \ \omega_{i}^{-1}$$','Interpreter','latex','FontSize',20)
% Kinetic subplot(1,2,1)
h1=subaxis(1,2,1,'SV',0,'SH',0.05,'MR',0.05,'ML',0.05,'PL',0.05,'PR',0.01);
plot(xll,smoothdata(total_dKedt_1d)./(VA2c*Pi0),'Color',colorblind(10,:), 'LineWidth', lw)
hold on
plot(xll,smoothdata(ve_Div_Pet_1d)./(VA2c*Pi0),'Color',colorblind(2,:), 'LineWidth', lw)
plot(xll,smoothdata(Ke_Div_ve_1d)./(VA2c*Pi0),'Color',colorblind(5,:), 'LineWidth', lw)
plot(xll,smoothdata(-qe_ne_ve_E_Gp_1d)./(VA2c*Pi0),'Color',colorblind(6,:), 'LineWidth', lw)
plot(xll,smoothdata(dum_pXIe_Gp_1d)./(VA2c*Pi0),'Color',colorblind(1,:), 'LineWidth', lw)
xline(4,'--k'); xline(8,'--k'); xline(5.7,'--b')
ax = gca; ax.FontSize = fs; ax.TickLabelInterpreter ='latex';
legend('$\frac{d \varepsilon^{k}_{e} /d t}{\Delta\varepsilon_{0}} $',...
    '$\frac{u_{e} \cdot \nabla \cdot \overline{P}_{e}}{\Delta\varepsilon_{0}}$',...
    '$\frac{\varepsilon^{k}_{e}\nabla \cdot u_{e}}{\Delta\varepsilon_{0}}$',...
    '$\frac{-q_{e}n_{e}E\cdot u_{e}}{\Delta\varepsilon_{0}}$',...
    '$\frac{m_{e}\mathbf{u}_{e} \cdot \mathbf{\Xi}_{e}}{\Delta\varepsilon_{0}}$',...
    'Interpreter','latex','FontSize',fs,...
    'Location','southeast','NumColumns',2)
xlabel('$r / d_{i}$','Interpreter','latex','FontSize',fs)
xlim([4,8])
%ylim([-0.005, 0.005])
ylim([-3.2, 3.2]) %<------------
ylabel('kinetic Energy rate terms','Interpreter','latex','FontSize',fs)
hold off
% Thermal
%------------------------------------
%subplot(1,2,2)
h1=subaxis(1,2,2,'SV',0,'SH',0.05,'MR',0.05,'ML',0.05,'PL',0.05,'PR',0.01);
plot(xll,smoothdata(total_dUedt_1d),'Color',colorblind(10,:), 'LineWidth', lw)
hold on
plot(xll,smoothdata(dkQiik_05_e_Gp_1d)./(VA2c*Pi0),'Color',colorblind(2,:), 'LineWidth', lw)
plot(xll,smoothdata(Pet_grad_ve_Gp_1d)./(VA2c*Pi0),'Color',colorblind(5,:), 'LineWidth', lw)
plot(xll,smoothdata(Ue_Div_ve_Gp_1d)./(VA2c*Pi0),'Color',colorblind(6,:), 'LineWidth', lw)
plot(xll,smoothdata(dum_pXI2e_Gp_1d)./(VA2c*Pi0),'Color',colorblind(1,:), 'LineWidth', lw)
xline(4,'--k'); xline(8,'--k'); xline(5.7,'--b')
ax = gca; ax.FontSize = fs; ax.TickLabelInterpreter ='latex';
legend('$\frac{d \varepsilon^{th}_{e} / d t }{\Delta\varepsilon_{0}}$',...
    '$\frac{\nabla \cdot h_{e}}{\Delta\varepsilon_{0}}$',...
    '$\frac{\nabla u_{e}:\overline{P}_{e}}{\Delta\varepsilon_{0}}$',...
    '$\frac{\varepsilon^{th}_{e} \nabla \cdot u_{e}}{\Delta\varepsilon_{0}}$',...
    '$\frac{0.5 Tr(m_{e}\overline{\mathbf{\Xi}}_{e})}{\Delta\varepsilon_{0}}$',...
    'Interpreter','latex','FontSize',fs,...
    'Location','southeast','NumColumns',2)
xlabel('$r / d_{i}$','Interpreter','latex','FontSize',fs)
xlim([4,8])
ylim([-10.2, 10.2]) %<--------------
%ylim([-0.015, 0.015])
ylabel('Thermal Energy rate terms','Interpreter','latex','FontSize',fs)
hold off

% f2155 control previous plot
%-----------------------------
%{
f2155 = figure(2155);
sgtitle('$$t = \ $$'+string(round(time_t,2)) + '$$ \ \omega_{i}^{-1}$$','Interpreter','latex','FontSize',20)
% Kinetic subplot(1,2,1)
h1=subaxis(1,2,1,'SV',0,'SH',0.05,'MR',0.05,'ML',0.05,'PL',0.05,'PR',0.01);
plot(xll,smoothdata(dtKe_in_Gp_1d)./(VA2c*Pi0),'Color',colorblind(10,:), 'LineWidth', lw)
hold on
plot(xll,smoothdata(ve_grad_Ke_1d)./(VA2c*Pi0),'Color',colorblind(1,:), 'LineWidth', lw)
plot(xll,smoothdata(ve_Div_Pet_1d)./(VA2c*Pi0),'Color',colorblind(2,:), 'LineWidth', lw)
plot(xll,smoothdata(Ke_Div_ve_1d)./(VA2c*Pi0),'Color',colorblind(5,:), 'LineWidth', lw)
plot(xll,smoothdata(-qe_ne_ve_E_Gp_1d)./(VA2c*Pi0),'Color',colorblind(6,:), 'LineWidth', lw)
xline(4,'--k'); xline(8,'--k'); xline(5.7,'--b')
ax = gca; ax.FontSize = fs; ax.TickLabelInterpreter ='latex';
legend('$\frac{\partial \varepsilon^{k}_{e} /\partial t }{\Delta\varepsilon_{0}}$',...
    '$\frac{u_{e} \cdot \nabla \varepsilon^{k}_{e}}{\Delta\varepsilon_{0}}$',...
    '$\frac{u_{e} \cdot \nabla \cdot \overline{P}_{e}}{\Delta\varepsilon_{0}}$',...
    '$\frac{\varepsilon^{k}_{e}\nabla \cdot u_{e}}{\Delta\varepsilon_{0}}$',...
    '$\frac{-q_{e}n_{e}E\cdot u_{e}}{\Delta\varepsilon_{0}}$',...
    'Interpreter','latex','FontSize',fs,...
    'Location','southeast','NumColumns',2)
xlabel('$r / d_{i}$','Interpreter','latex','FontSize',fs)
xlim([4,8])
%ylim([-0.005, 0.005])
ylim([-3.2, 3.2]) %<------------
ylabel('kinetic Energy rate terms','Interpreter','latex','FontSize',fs)
hold off
% Thermal
%------------------------------------
%subplot(1,2,2)
h1=subaxis(1,2,2,'SV',0,'SH',0.05,'MR',0.05,'ML',0.05,'PL',0.05,'PR',0.01);
plot(xll,smoothdata(dtUe_in_Gp_1d)./(VA2c*Pi0),'Color',colorblind(10,:), 'LineWidth', lw)
hold on
plot(xll,smoothdata(ve_grad_Ue_1d)./(VA2c*Pi0),'Color',colorblind(1,:), 'LineWidth', lw)
plot(xll,smoothdata(dkQiik_05_e_Gp_1d)./(VA2c*Pi0),'Color',colorblind(2,:), 'LineWidth', lw)
plot(xll,smoothdata(Pet_grad_ve_Gp_1d)./(VA2c*Pi0),'Color',colorblind(5,:), 'LineWidth', lw)
plot(xll,smoothdata(Pet_grad_ve_ii_Gp_1d)./(VA2c*Pi0),'Color',colorblind(6,:), 'LineWidth', lw)
plot(xll,smoothdata(Pet_grad_ve_ij_Gp_1d)./(VA2c*Pi0),'Color',colorblind(8,:), 'LineWidth', lw)
plot(xll,smoothdata(Ue_Div_ve_Gp_1d)./(VA2c*Pi0),'Color',colorblind(9,:), 'LineWidth', lw)
xline(4,'--k'); xline(8,'--k'); xline(5.7,'--b')
ax = gca; ax.FontSize = fs; ax.TickLabelInterpreter ='latex';
legend('$\frac{\partial \varepsilon^{th}_{e} / \partial t }{\Delta\varepsilon_{0}}$',...
    '$\frac{u_{e} \cdot \nabla \varepsilon^{th}_{e}}{\Delta\varepsilon_{0}}$',...
    '$\frac{\nabla \cdot h_{e}}{\Delta\varepsilon_{0}}$',...
    '$\frac{\nabla u_{e}:\overline{P}_{e}}{\Delta\varepsilon_{0}}$',...
    '$\frac{p\theta_{e}}{\Delta\varepsilon_{0}}$',...
    '$\frac{PiD_{e}}{\Delta\varepsilon_{0}}$',...
    '$\frac{\varepsilon^{th}_{e} \nabla \cdot u_{e}}{\Delta\varepsilon_{0}}$',...
    'Interpreter','latex','FontSize',fs,...
    'Location','southeast','NumColumns',2)
xlabel('$r / d_{i}$','Interpreter','latex','FontSize',fs)
xlim([4,8])
ylim([-8.2, 8.2]) %<--------------
%ylim([-0.015, 0.015])
ylabel('Thermal Energy rate terms','Interpreter','latex','FontSize',fs)
hold off
%}
%-----------------------------
%}

%--------------------------------------------------------------------------
% Kinetic and thermal energy terms for Ions
%--------------------------------------------------------------------------
%{
f2156 = figure(2156);
sgtitle('$$t = \ $$'+string(round(time_t,2)) + '$$ \ \omega_{i}^{-1}$$','Interpreter','latex','FontSize',20)
% Kinetic subplot(1,2,1)
h1=subaxis(1,2,1,'SV',0,'SH',0.05,'MR',0.05,'ML',0.05,'PL',0.05,'PR',0.01);
plot(xll,smoothdata(total_dKidt_1d)./(VA2c*Pi0),'Color',colorblind(10,:), 'LineWidth', lw)
hold on
plot(xll,smoothdata(vi_Div_Pit_1d)./(VA2c*Pi0),'Color',colorblind(2,:), 'LineWidth', lw)
plot(xll,smoothdata(Ki_Div_vi_1d)./(VA2c*Pi0),'Color',colorblind(5,:), 'LineWidth', lw)
plot(xll,smoothdata(-qi_ni_vi_E_Gp_1d)./(VA2c*Pi0),'Color',colorblind(6,:), 'LineWidth', lw)
plot(xll,smoothdata(dum_pXIi_Gp_1d)./(VA2c*Pi0),'Color',colorblind(1,:), 'LineWidth', lw)
xline(4,'--k'); xline(8,'--k'); xline(5.7,'--b')
ax = gca; ax.FontSize = fs; ax.TickLabelInterpreter ='latex';
legend('$\frac{d \varepsilon^{k}_{i} /d t}{\Delta\varepsilon_{0}} $',...
    '$\frac{u_{i} \cdot \nabla \cdot \overline{P}_{i}}{\Delta\varepsilon_{0}}$',...
    '$\frac{\varepsilon^{k}_{i}\nabla \cdot u_{i}}{\Delta\varepsilon_{0}}$',...
    '$\frac{-q_{i}n_{i}E\cdot u_{i}}{\Delta\varepsilon_{0}}$',...
    '$\frac{m_{i}\mathbf{u}_{i} \cdot \mathbf{\Xi}_{i}}{\Delta\varepsilon_{0}}$',...
    'Interpreter','latex','FontSize',fs,...
    'Location','southeast','NumColumns',2)
xlabel('$r / d_{i}$','Interpreter','latex','FontSize',fs)
%xlim([4,8])
%ylim([-0.005, 0.005])
ylim([-3.2, 3.2]) %<------------
ylabel('kinetic Energy rate terms','Interpreter','latex','FontSize',fs)
hold off
% Thermal
%------------------------------------
%subplot(1,2,2)
h1=subaxis(1,2,2,'SV',0,'SH',0.05,'MR',0.05,'ML',0.05,'PL',0.05,'PR',0.01);
plot(xll,smoothdata(total_dUidt_1d),'Color',colorblind(10,:), 'LineWidth', lw)
hold on
plot(xll,smoothdata(dkQiik_05_i_Gp_1d)./(VA2c*Pi0),'Color',colorblind(2,:), 'LineWidth', lw)
plot(xll,smoothdata(Pit_grad_vi_Gp_1d)./(VA2c*Pi0),'Color',colorblind(5,:), 'LineWidth', lw)
plot(xll,smoothdata(Ui_Div_vi_Gp_1d)./(VA2c*Pi0),'Color',colorblind(6,:), 'LineWidth', lw)
plot(xll,smoothdata(dum_pXI2i_Gp_1d)./(VA2c*Pi0),'Color',colorblind(1,:), 'LineWidth', lw)
xline(4,'--k'); xline(8,'--k'); xline(5.7,'--b')
ax = gca; ax.FontSize = fs; ax.TickLabelInterpreter ='latex';
legend('$\frac{d \varepsilon^{th}_{i} / d t }{\Delta\varepsilon_{0}}$',...
    '$\frac{\nabla \cdot h_{i}}{\Delta\varepsilon_{0}}$',...
    '$\frac{\nabla u_{i}:\overline{P}_{i}}{\Delta\varepsilon_{0}}$',...
    '$\frac{\varepsilon^{th}_{i} \nabla \cdot u_{i}}{\Delta\varepsilon_{0}}$',...
    '$\frac{0.5 Tr(m_{i}\overline{\mathbf{\Xi}}_{i})}{\Delta\varepsilon_{0}}$',...
    'Interpreter','latex','FontSize',fs,...
    'Location','southeast','NumColumns',2)
xlabel('$r / d_{i}$','Interpreter','latex','FontSize',fs)
%xlim([4,8])
ylim([-10.2, 10.2]) %<--------------
%ylim([-0.015, 0.015])
ylabel('Thermal Energy rate terms','Interpreter','latex','FontSize',fs)
hold off

f21556 = figure(21556);
sgtitle('$$t = \ $$'+string(round(time_t,2)) + '$$ \ \omega_{i}^{-1}$$','Interpreter','latex','FontSize',20)
% Kinetic subplot(1,2,1)
h1=subaxis(1,2,1,'SV',0,'SH',0.05,'MR',0.05,'ML',0.05,'PL',0.05,'PR',0.01);
plot(xll,smoothdata(total_dKidt_1d)./(VA2c*Pi0),'Color',colorblind(10,:), 'LineWidth', lw)
hold on
plot(xll,smoothdata(vi_Div_Pit_1d)./(VA2c*Pi0),'Color',colorblind(2,:), 'LineWidth', lw)
plot(xll,smoothdata(Ki_Div_vi_1d)./(VA2c*Pi0),'Color',colorblind(5,:), 'LineWidth', lw)
plot(xll,smoothdata(-qi_ni_vi_E_Gp_1d)./(VA2c*Pi0),'Color',colorblind(6,:), 'LineWidth', lw)
plot(xll,smoothdata(dum_pXIi_Gp_1d)./(VA2c*Pi0),'Color',colorblind(1,:), 'LineWidth', lw)
xline(4,'--k'); xline(8,'--k'); xline(5.7,'--b')
ax = gca; ax.FontSize = fs; ax.TickLabelInterpreter ='latex';
legend('$\frac{d \varepsilon^{k}_{i} /d t}{\Delta\varepsilon_{0}} $',...
    '$\frac{u_{i} \cdot \nabla \cdot \overline{P}_{i}}{\Delta\varepsilon_{0}}$',...
    '$\frac{\varepsilon^{k}_{i}\nabla \cdot u_{i}}{\Delta\varepsilon_{0}}$',...
    '$\frac{-q_{i}n_{i}E\cdot u_{i}}{\Delta\varepsilon_{0}}$',...
    '$\frac{m_{i}\mathbf{u}_{i} \cdot \mathbf{\Xi}_{i}}{\Delta\varepsilon_{0}}$',...
    'Interpreter','latex','FontSize',fs,...
    'Location','southeast','NumColumns',2)
xlabel('$r / d_{i}$','Interpreter','latex','FontSize',fs)
xlim([4,8])
%ylim([-0.005, 0.005])
ylim([-3.2, 3.2]) %<------------
ylabel('kinetic Energy rate terms','Interpreter','latex','FontSize',fs)
hold off
% Thermal
%------------------------------------
%subplot(1,2,2)
h1=subaxis(1,2,2,'SV',0,'SH',0.05,'MR',0.05,'ML',0.05,'PL',0.05,'PR',0.01);
plot(xll,smoothdata(total_dUidt_1d),'Color',colorblind(10,:), 'LineWidth', lw)
hold on
plot(xll,smoothdata(dkQiik_05_i_Gp_1d)./(VA2c*Pi0),'Color',colorblind(2,:), 'LineWidth', lw)
plot(xll,smoothdata(Pit_grad_vi_Gp_1d)./(VA2c*Pi0),'Color',colorblind(5,:), 'LineWidth', lw)
plot(xll,smoothdata(Ui_Div_vi_Gp_1d)./(VA2c*Pi0),'Color',colorblind(6,:), 'LineWidth', lw)
plot(xll,smoothdata(dum_pXI2i_Gp_1d)./(VA2c*Pi0),'Color',colorblind(1,:), 'LineWidth', lw)
xline(4,'--k'); xline(8,'--k'); xline(5.7,'--b')
ax = gca; ax.FontSize = fs; ax.TickLabelInterpreter ='latex';
legend('$\frac{d \varepsilon^{th}_{i} / d t }{\Delta\varepsilon_{0}}$',...
    '$\frac{\nabla \cdot h_{i}}{\Delta\varepsilon_{0}}$',...
    '$\frac{\nabla u_{i}:\overline{P}_{i}}{\Delta\varepsilon_{0}}$',...
    '$\frac{\varepsilon^{th}_{i} \nabla \cdot u_{i}}{\Delta\varepsilon_{0}}$',...
    '$\frac{0.5 Tr(m_{i}\overline{\mathbf{\Xi}}_{i})}{\Delta\varepsilon_{0}}$',...
    'Interpreter','latex','FontSize',fs,...
    'Location','southeast','NumColumns',2)
xlabel('$r / d_{i}$','Interpreter','latex','FontSize',fs)
xlim([4,8])
ylim([-10.2, 10.2]) %<--------------
%ylim([-0.015, 0.015])
ylabel('Thermal Energy rate terms','Interpreter','latex','FontSize',fs)
hold off
%}
%--------------------------------------------------------------------------

% Kinetic and thermal energy both electrons and ions
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
f27552 = figure(27552);
sgtitle('$$t = \ $$'+string(round(time_t,2)) + '$$ \ \omega_{i}^{-1}$$','Interpreter','latex','FontSize',20)
h1=subaxis(2,2,1,'SV',0.02,'SH',0.0,'MR',0.0,'ML',0.0,'PL',0.05,'PR',0.01);
plot(yll,smoothdata(total_dKedt_1d)./(VA2c*Pi0),'Color',colorblind(10,:), 'LineWidth', lw)
hold on
plot(yll,smoothdata(ve_Div_Pet_1d)./(VA2c*Pi0),'Color',colorblind(2,:), 'LineWidth', lw)
plot(yll,smoothdata(Ke_Div_ve_1d)./(VA2c*Pi0),'Color',colorblind(5,:), 'LineWidth', lw)
plot(yll,smoothdata(-qe_ne_ve_E_Gp_1d)./(VA2c*Pi0),'Color',colorblind(6,:), 'LineWidth', lw)
plot(yll,smoothdata(dum_pXIe_Gp_1d)./(VA2c*Pi0),'Color',colorblind(1,:), 'LineWidth', lw)
xline(4,'--k'); xline(8,'--k'); xline(5.7,'--b')
ax = gca; ax.FontSize = fs; ax.TickLabelInterpreter ='latex';
legend('$\frac{d \varepsilon^{k}_{e} /d t}{\Delta\varepsilon_{0}} $',...
    '$\frac{u_{e} \cdot \nabla \cdot \overline{P}_{e}}{\Delta\varepsilon_{0}}$',...
    '$\frac{\varepsilon^{k}_{e}\nabla \cdot u_{e}}{\Delta\varepsilon_{0}}$',...
    '$\frac{-q_{e}n_{e}E\cdot u_{e}}{\Delta\varepsilon_{0}}$',...
    '$\frac{m_{e}\mathbf{u}_{e} \cdot \mathbf{\Xi}_{e}}{\Delta\varepsilon_{0}}$',...
    'Interpreter','latex','FontSize',fs,...
    'Location','northeast','NumColumns',2)
%xlabel('$r / d_{i}$','Interpreter','latex','FontSize',fs)
xlim([4,8])
xticks([]);
ylim([-3.2, 3.2]) %<------------
ylabel('kinetic Energy rate terms','Interpreter','latex','FontSize',fs)
hold off
% Thermal
%------------------------------------
h1=subaxis(2,2,2,'SV',0.02,'SH',0.0,'MR',0.05,'ML',0.0,'PL',0.03,'PR',0.01);
plot(yll,smoothdata(total_dUedt_1d),'Color',colorblind(10,:), 'LineWidth', lw)
hold on
plot(yll,smoothdata(dkQiik_05_e_Gp_1d)./(VA2c*Pi0),'Color',colorblind(2,:), 'LineWidth', lw)
plot(yll,smoothdata(Pet_grad_ve_Gp_1d)./(VA2c*Pi0),'Color',colorblind(5,:), 'LineWidth', lw)
plot(yll,smoothdata(Ue_Div_ve_Gp_1d)./(VA2c*Pi0),'Color',colorblind(6,:), 'LineWidth', lw)
plot(yll,smoothdata(dum_pXI2e_Gp_1d)./(VA2c*Pi0),'Color',colorblind(1,:), 'LineWidth', lw)
xline(4,'--k'); xline(8,'--k'); xline(5.7,'--b')
ax = gca; ax.FontSize = fs; ax.TickLabelInterpreter ='latex';
%%%%%%%%%%%%
ax.YAxisLocation = 'right';
%%%%%%%%%%%%
legend('$\frac{d \varepsilon^{th}_{e} / d t }{\Delta\varepsilon_{0}}$',...
    '$\frac{\nabla \cdot h_{e}}{\Delta\varepsilon_{0}}$',...
    '$\frac{\nabla u_{e}:\overline{P}_{e}}{\Delta\varepsilon_{0}}$',...
    '$\frac{\varepsilon^{th}_{e} \nabla \cdot u_{e}}{\Delta\varepsilon_{0}}$',...
    '$\frac{0.5 Tr(m_{e}\overline{\mathbf{\Xi}}_{e})}{\Delta\varepsilon_{0}}$',...
    'Interpreter','latex','FontSize',fs,...
    'Location','northeast','NumColumns',2)
%xlabel('$r / d_{i}$','Interpreter','latex','FontSize',fs)
xticks([])
xlim([4,8])
ylim([-15.2, 15.2]) %<--------------
ylabel('Thermal Energy rate terms','Interpreter','latex','FontSize',fs)
hold off
%--------------------------------------------------------------------------
% Kinetic subplot(1,2,1)
h1=subaxis(2,2,3,'SV',0.02,'SH',0.0,'MR',0.00,'ML',0.00,'PL',0.05,'PR',0.01);
plot(yll,smoothdata(total_dKidt_1d)./(VA2c*Pi0),'Color',colorblind(10,:), 'LineWidth', lw)
hold on
plot(yll,smoothdata(vi_Div_Pit_1d)./(VA2c*Pi0),'Color',colorblind(2,:), 'LineWidth', lw)
plot(yll,smoothdata(Ki_Div_vi_1d)./(VA2c*Pi0),'Color',colorblind(5,:), 'LineWidth', lw)
plot(yll,smoothdata(-qi_ni_vi_E_Gp_1d)./(VA2c*Pi0),'Color',colorblind(6,:), 'LineWidth', lw)
plot(yll,smoothdata(dum_pXIi_Gp_1d)./(VA2c*Pi0),'Color',colorblind(1,:), 'LineWidth', lw)
xline(4,'--k'); xline(8,'--k'); xline(5.7,'--b')
ax = gca; ax.FontSize = fs; ax.TickLabelInterpreter ='latex';
legend('$\frac{d \varepsilon^{k}_{i} /d t}{\Delta\varepsilon_{0}} $',...
    '$\frac{u_{i} \cdot \nabla \cdot \overline{P}_{i}}{\Delta\varepsilon_{0}}$',...
    '$\frac{\varepsilon^{k}_{i}\nabla \cdot u_{i}}{\Delta\varepsilon_{0}}$',...
    '$\frac{-q_{i}n_{i}E\cdot u_{i}}{\Delta\varepsilon_{0}}$',...
    '$\frac{m_{i}\mathbf{u}_{i} \cdot \mathbf{\Xi}_{i}}{\Delta\varepsilon_{0}}$',...
    'Interpreter','latex','FontSize',fs,...
    'Location','northeast','NumColumns',2)
xlabel('$r / d_{i}$','Interpreter','latex','FontSize',fs)
xlim([4,8])
ylim([-3.2, 3.2]) %<------------
ylabel('kinetic Energy rate terms','Interpreter','latex','FontSize',fs)
hold off
% Thermal
%------------------------------------
h1=subaxis(2,2,4,'SV',0.02,'SH',0.0,'MR',0.05,'ML',0.00,'PL',0.03,'PR',0.01);
%h1=subaxis(2,2,2,'SV',0.01,'SH',0.0,'MR',0.05,'ML',0.0,'PL',0.02,'PR',0.01);
plot(yll,smoothdata(total_dUidt_1d),'Color',colorblind(10,:), 'LineWidth', lw)
hold on
plot(yll,smoothdata(dkQiik_05_i_Gp_1d)./(VA2c*Pi0),'Color',colorblind(2,:), 'LineWidth', lw)
plot(yll,smoothdata(Pit_grad_vi_Gp_1d)./(VA2c*Pi0),'Color',colorblind(5,:), 'LineWidth', lw)
plot(yll,smoothdata(Ui_Div_vi_Gp_1d)./(VA2c*Pi0),'Color',colorblind(6,:), 'LineWidth', lw)
plot(yll,smoothdata(dum_pXI2i_Gp_1d)./(VA2c*Pi0),'Color',colorblind(1,:), 'LineWidth', lw)
xline(4,'--k'); xline(8,'--k'); xline(5.7,'--b')
ax = gca; ax.FontSize = fs; ax.TickLabelInterpreter ='latex';
%%%%%%%%%%%%
ax.YAxisLocation = 'right';
%%%%%%%%%%%%
legend('$\frac{d \varepsilon^{th}_{i} / d t }{\Delta\varepsilon_{0}}$',...
    '$\frac{\nabla \cdot h_{i}}{\Delta\varepsilon_{0}}$',...
    '$\frac{\nabla u_{i}:\overline{P}_{i}}{\Delta\varepsilon_{0}}$',...
    '$\frac{\varepsilon^{th}_{i} \nabla \cdot u_{i}}{\Delta\varepsilon_{0}}$',...
    '$\frac{0.5 Tr(m_{i}\overline{\mathbf{\Xi}}_{i})}{\Delta\varepsilon_{0}}$',...
    'Interpreter','latex','FontSize',fs,...
    'Location','northeast','NumColumns',2)
xlabel('$r / d_{i}$','Interpreter','latex','FontSize',fs)
xlim([4,8])
ylim([-10.2, 10.2]) %<--------------
ylabel('Thermal Energy rate terms','Interpreter','latex','FontSize',fs)
hold off


f275512 = figure(275512);
sgtitle('$$t = \ $$'+string(round(time_t,2)) + '$$ \ \omega_{i}^{-1}$$','Interpreter','latex','FontSize',20)
h1=subaxis(2,2,1,'SV',0.02,'SH',0.0,'MR',0.0,'ML',0.0,'PL',0.05,'PR',0.01);
plot(yll,smoothdata(total_dKedt_1d)./(VA2c*Pi0),'Color',colorblind(10,:), 'LineWidth', lw)
hold on
plot(yll,smoothdata(ve_Div_Pet_1d)./(VA2c*Pi0),'Color',colorblind(2,:), 'LineWidth', lw)
plot(yll,smoothdata(Ke_Div_ve_1d)./(VA2c*Pi0),'Color',colorblind(5,:), 'LineWidth', lw)
plot(yll,smoothdata(-qe_ne_ve_E_Gp_1d)./(VA2c*Pi0),'Color',colorblind(6,:), 'LineWidth', lw)
plot(yll,smoothdata(dum_pXIe_Gp_1d)./(VA2c*Pi0),'Color',colorblind(1,:), 'LineWidth', lw)
xline(4,'--k'); xline(8,'--k'); xline(5.7,'--b')
ax = gca; ax.FontSize = fs; ax.TickLabelInterpreter ='latex';
legend('$\frac{d \varepsilon^{k}_{e} /d t}{\Delta\varepsilon_{0}} $',...
    '$\frac{u_{e} \cdot \nabla \cdot \overline{P}_{e}}{\Delta\varepsilon_{0}}$',...
    '$\frac{\varepsilon^{k}_{e}\nabla \cdot u_{e}}{\Delta\varepsilon_{0}}$',...
    '$\frac{-q_{e}n_{e}E\cdot u_{e}}{\Delta\varepsilon_{0}}$',...
    '$\frac{m_{e}\mathbf{u}_{e} \cdot \mathbf{\Xi}_{e}}{\Delta\varepsilon_{0}}$',...
    'Interpreter','latex','FontSize',fs,...
    'Location','northeast','NumColumns',2)
%xlabel('$r / d_{i}$','Interpreter','latex','FontSize',fs)
xlim([0,10])
xticks([]);
ylim([-3.2, 3.2]) %<------------
ylabel('kinetic Energy rate terms','Interpreter','latex','FontSize',fs)
hold off
% Thermal
%------------------------------------
h1=subaxis(2,2,2,'SV',0.02,'SH',0.0,'MR',0.05,'ML',0.0,'PL',0.03,'PR',0.01);
plot(yll,smoothdata(total_dUedt_1d),'Color',colorblind(10,:), 'LineWidth', lw)
hold on
plot(yll,smoothdata(dkQiik_05_e_Gp_1d)./(VA2c*Pi0),'Color',colorblind(2,:), 'LineWidth', lw)
plot(yll,smoothdata(Pet_grad_ve_Gp_1d)./(VA2c*Pi0),'Color',colorblind(5,:), 'LineWidth', lw)
plot(yll,smoothdata(Ue_Div_ve_Gp_1d)./(VA2c*Pi0),'Color',colorblind(6,:), 'LineWidth', lw)
plot(yll,smoothdata(dum_pXI2e_Gp_1d)./(VA2c*Pi0),'Color',colorblind(1,:), 'LineWidth', lw)
xline(4,'--k'); xline(8,'--k'); xline(5.7,'--b')
ax = gca; ax.FontSize = fs; ax.TickLabelInterpreter ='latex';
%%%%%%%%%%%%
ax.YAxisLocation = 'right';
%%%%%%%%%%%%
legend('$\frac{d \varepsilon^{th}_{e} / d t }{\Delta\varepsilon_{0}}$',...
    '$\frac{\nabla \cdot h_{e}}{\Delta\varepsilon_{0}}$',...
    '$\frac{\nabla u_{e}:\overline{P}_{e}}{\Delta\varepsilon_{0}}$',...
    '$\frac{\varepsilon^{th}_{e} \nabla \cdot u_{e}}{\Delta\varepsilon_{0}}$',...
    '$\frac{0.5 Tr(m_{e}\overline{\mathbf{\Xi}}_{e})}{\Delta\varepsilon_{0}}$',...
    'Interpreter','latex','FontSize',fs,...
    'Location','northeast','NumColumns',2)
%xlabel('$r / d_{i}$','Interpreter','latex','FontSize',fs)
xticks([])
xlim([0,10])
ylim([-15.2, 15.2]) %<--------------
ylabel('Thermal Energy rate terms','Interpreter','latex','FontSize',fs)
hold off
%--------------------------------------------------------------------------
% Kinetic subplot(1,2,1)
h1=subaxis(2,2,3,'SV',0.02,'SH',0.0,'MR',0.00,'ML',0.00,'PL',0.05,'PR',0.01);
plot(yll,smoothdata(total_dKidt_1d)./(VA2c*Pi0),'Color',colorblind(10,:), 'LineWidth', lw)
hold on
plot(yll,smoothdata(vi_Div_Pit_1d)./(VA2c*Pi0),'Color',colorblind(2,:), 'LineWidth', lw)
plot(yll,smoothdata(Ki_Div_vi_1d)./(VA2c*Pi0),'Color',colorblind(5,:), 'LineWidth', lw)
plot(yll,smoothdata(-qi_ni_vi_E_Gp_1d)./(VA2c*Pi0),'Color',colorblind(6,:), 'LineWidth', lw)
plot(yll,smoothdata(dum_pXIi_Gp_1d)./(VA2c*Pi0),'Color',colorblind(1,:), 'LineWidth', lw)
xline(4,'--k'); xline(8,'--k'); xline(5.7,'--b')
ax = gca; ax.FontSize = fs; ax.TickLabelInterpreter ='latex';
legend('$\frac{d \varepsilon^{k}_{i} /d t}{\Delta\varepsilon_{0}} $',...
    '$\frac{u_{i} \cdot \nabla \cdot \overline{P}_{i}}{\Delta\varepsilon_{0}}$',...
    '$\frac{\varepsilon^{k}_{i}\nabla \cdot u_{i}}{\Delta\varepsilon_{0}}$',...
    '$\frac{-q_{i}n_{i}E\cdot u_{i}}{\Delta\varepsilon_{0}}$',...
    '$\frac{m_{i}\mathbf{u}_{i} \cdot \mathbf{\Xi}_{i}}{\Delta\varepsilon_{0}}$',...
    'Interpreter','latex','FontSize',fs,...
    'Location','northeast','NumColumns',2)
xlabel('$r / d_{i}$','Interpreter','latex','FontSize',fs)
xlim([0,10])
ylim([-3.2, 3.2]) %<------------
ylabel('kinetic Energy rate terms','Interpreter','latex','FontSize',fs)
hold off
% Thermal
%------------------------------------
h1=subaxis(2,2,4,'SV',0.02,'SH',0.0,'MR',0.05,'ML',0.00,'PL',0.03,'PR',0.01);
%h1=subaxis(2,2,2,'SV',0.01,'SH',0.0,'MR',0.05,'ML',0.0,'PL',0.02,'PR',0.01);
plot(yll,smoothdata(total_dUidt_1d),'Color',colorblind(10,:), 'LineWidth', lw)
hold on
plot(yll,smoothdata(dkQiik_05_i_Gp_1d)./(VA2c*Pi0),'Color',colorblind(2,:), 'LineWidth', lw)
plot(yll,smoothdata(Pit_grad_vi_Gp_1d)./(VA2c*Pi0),'Color',colorblind(5,:), 'LineWidth', lw)
plot(yll,smoothdata(Ui_Div_vi_Gp_1d)./(VA2c*Pi0),'Color',colorblind(6,:), 'LineWidth', lw)
plot(yll,smoothdata(dum_pXI2i_Gp_1d)./(VA2c*Pi0),'Color',colorblind(1,:), 'LineWidth', lw)
xline(4,'--k'); xline(8,'--k'); xline(5.7,'--b')
ax = gca; ax.FontSize = fs; ax.TickLabelInterpreter ='latex';
%%%%%%%%%%%%
ax.YAxisLocation = 'right';
%%%%%%%%%%%%
legend('$\frac{d \varepsilon^{th}_{i} / d t }{\Delta\varepsilon_{0}}$',...
    '$\frac{\nabla \cdot h_{i}}{\Delta\varepsilon_{0}}$',...
    '$\frac{\nabla u_{i}:\overline{P}_{i}}{\Delta\varepsilon_{0}}$',...
    '$\frac{\varepsilon^{th}_{i} \nabla \cdot u_{i}}{\Delta\varepsilon_{0}}$',...
    '$\frac{0.5 Tr(m_{i}\overline{\mathbf{\Xi}}_{i})}{\Delta\varepsilon_{0}}$',...
    'Interpreter','latex','FontSize',fs,...
    'Location','northeast','NumColumns',2)
xlabel('$r / d_{i}$','Interpreter','latex','FontSize',fs)
xlim([0,10])
ylim([-10.2, 10.2]) %<--------------
ylabel('Thermal Energy rate terms','Interpreter','latex','FontSize',fs)
hold off
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------

%There is no difference. The only difference is the notation.
%--------------------------------------------------------------------------
%{
f216 = figure(216);
sgtitle('$$t = \ $$'+string(round(time_t,2)) + '$$ \ \omega_{i}^{-1}$$','Interpreter','latex','FontSize',20)
plot(xll,smoothdata(Div_Pes_dot_ve_Gp_1d)./(VA2c*Pi0),'Color',colorblind(1,:), 'LineWidth', lw)
hold on
plot(xll,smoothdata(Pet_grad_ve_Gp_1d)./(VA2c*Pi0),'Color',colorblind(2,:), 'LineWidth', lw)
plot(xll,smoothdata(ve_Div_Pet_1d)./(VA2c*Pi0),'Color',colorblind(5,:), 'LineWidth', lw)
%plot(xll,Pet_grad_ve_Gp_1d + ve_Div_Pet_1d,'g')
xline(4,'--k'); xline(8,'--k'); xline(5.7,'--b')
ax = gca; ax.FontSize = fs; ax.TickLabelInterpreter ='latex';
%legend('dKe/dt','vedivPe','Kedivve','-qEe')
legend('$\frac{\nabla \cdot( \overline{P}_{e} \cdot u_{e})}{\Delta\varepsilon_{0}}$ ',...
    '$\frac{\overline{P}_{e}:\nabla u_{e}}{\Delta\varepsilon_{0}}$',...
    '$\frac{u_{e} \cdot(\nabla \cdot \overline{P}_{e})}{\Delta\varepsilon_{0}}$',...
    'Interpreter','latex','FontSize',fs)
%plot(xll,ve_2_1d)
xlabel('$r / d_{i}$','Interpreter','latex','FontSize',fs)
ylabel('Kinetic Energy rate terms','Interpreter','latex','FontSize',fs)
ylim([-10.2, 10.2])
%xlim([4,8])
hold off

f2165 = figure(2165);
sgtitle('$$t = \ $$'+string(round(time_t,2)) + '$$ \ \omega_{i}^{-1}$$','Interpreter','latex','FontSize',20)
plot(xll,smoothdata(Div_Pes_dot_ve_Gp_1d)./(VA2c*Pi0),'Color',colorblind(1,:), 'LineWidth', lw)
hold on
plot(xll,smoothdata(Pet_grad_ve_Gp_1d)./(VA2c*Pi0),'Color',colorblind(2,:), 'LineWidth', lw)
plot(xll,smoothdata(ve_Div_Pet_1d)./(VA2c*Pi0),'Color',colorblind(5,:), 'LineWidth', lw)
%plot(xll,Pet_grad_ve_Gp_1d + ve_Div_Pet_1d,'g')
xline(4,'--k'); xline(8,'--k'); xline(5.7,'--b')
ax = gca; ax.FontSize = fs; ax.TickLabelInterpreter ='latex';
%legend('dKe/dt','vedivPe','Kedivve','-qEe')
legend('$\frac{\nabla \cdot( \overline{P}_{e} \cdot u_{e})}{\Delta\varepsilon_{0}}$ ',...
    '$\frac{\overline{P}_{e}:\nabla u_{e}}{\Delta\varepsilon_{0}}$',...
    '$\frac{u_{e} \cdot(\nabla \cdot \overline{P}_{e})}{\Delta\varepsilon_{0}}$',...
    'Interpreter','latex','FontSize',fs)
%plot(xll,ve_2_1d)
xlabel('$r / d_{i}$','Interpreter','latex','FontSize',fs)
ylabel('Kinetic Energy rate terms','Interpreter','latex','FontSize',fs)
ylim([-10.2, 10.2])
xlim([4,8])
hold off
%}
%--------------------------------------------------------------------------

% For Ions
%--------------------------------------------------------------------------
%{
f2166 = figure(2166);
sgtitle('$$t = \ $$'+string(round(time_t,2)) + '$$ \ \omega_{i}^{-1}$$','Interpreter','latex','FontSize',20)
plot(xll,smoothdata(Div_Pis_dot_vi_Gp_1d)./(VA2c*Pi0),'Color',colorblind(1,:), 'LineWidth', lw)
hold on
plot(xll,smoothdata(Pit_grad_vi_Gp_1d)./(VA2c*Pi0),'Color',colorblind(2,:), 'LineWidth', lw)
plot(xll,smoothdata(vi_Div_Pit_1d)./(VA2c*Pi0),'Color',colorblind(5,:), 'LineWidth', lw)
%plot(xll,Pet_grad_ve_Gp_1d + ve_Div_Pet_1d,'g')
xline(4,'--k'); xline(8,'--k'); xline(5.7,'--b')
ax = gca; ax.FontSize = fs; ax.TickLabelInterpreter ='latex';
%legend('dKe/dt','vedivPe','Kedivve','-qEe')
legend('$\frac{\nabla \cdot( \overline{P}_{i} \cdot u_{i})}{\Delta\varepsilon_{0}}$ ',...
    '$\frac{\overline{P}_{i}:\nabla u_{i}}{\Delta\varepsilon_{0}}$',...
    '$\frac{u_{i} \cdot(\nabla \cdot \overline{P}_{i})}{\Delta\varepsilon_{0}}$',...
    'Interpreter','latex','FontSize',fs)
%plot(xll,ve_2_1d)
xlabel('$r / d_{i}$','Interpreter','latex','FontSize',fs)
ylabel('Kinetic Energy rate terms','Interpreter','latex','FontSize',fs)
ylim([-3.2, 3.2])
%xlim([4,8])
hold off

f21656 = figure(21656);
sgtitle('$$t = \ $$'+string(round(time_t,2)) + '$$ \ \omega_{i}^{-1}$$','Interpreter','latex','FontSize',20)
plot(xll,smoothdata(Div_Pis_dot_vi_Gp_1d)./(VA2c*Pi0),'Color',colorblind(1,:), 'LineWidth', lw)
hold on
plot(xll,smoothdata(Pit_grad_vi_Gp_1d)./(VA2c*Pi0),'Color',colorblind(2,:), 'LineWidth', lw)
plot(xll,smoothdata(vi_Div_Pit_1d)./(VA2c*Pi0),'Color',colorblind(5,:), 'LineWidth', lw)
%plot(xll,Pet_grad_ve_Gp_1d + ve_Div_Pet_1d,'g')
xline(4,'--k'); xline(8,'--k'); xline(5.7,'--b')
ax = gca; ax.FontSize = fs; ax.TickLabelInterpreter ='latex';
%legend('dKe/dt','vedivPe','Kedivve','-qEe')
legend('$\frac{\nabla \cdot( \overline{P}_{i} \cdot u_{i})}{\Delta\varepsilon_{0}}$ ',...
    '$\frac{\overline{P}_{i}:\nabla u_{i}}{\Delta\varepsilon_{0}}$',...
    '$\frac{u_{i} \cdot(\nabla \cdot \overline{P}_{i})}{\Delta\varepsilon_{0}}$',...
    'Interpreter','latex','FontSize',fs)
%plot(xll,ve_2_1d)
xlabel('$r / d_{i}$','Interpreter','latex','FontSize',fs)
ylabel('Kinetic Energy rate terms','Interpreter','latex','FontSize',fs)
ylim([-3.2, 3.2])
xlim([4,8])
hold off
%}
%--------------------------------------------------------------------------

%{
f217 = figure(217);
sgtitle('$$t = \ $$'+string(round(time_t,2)) + '$$ \ \omega_{i}^{-1}$$','Interpreter','latex','FontSize',20)
plot(xll,smoothdata(Pet_grad_ve_Gp_1d)./(VA2c*Pi0),'Color',colorblind(1,:), 'LineWidth', lw)
hold on
plot(xll,smoothdata(p_theta_e_Gp_1d)./(VA2c*Pi0),'Color',colorblind(2,:), 'LineWidth', lw)
plot(xll,smoothdata(PiD_e_Gp_1d)./(VA2c*Pi0),'Color',colorblind(5,:), 'LineWidth', lw)
%plot(xll,Pet_grad_ve_Gp_1d + ve_Div_Pet_1d,'g')
xline(4,'--k'); xline(8,'--k'); xline(5.7,'--b')
ax = gca; ax.FontSize = fs; ax.TickLabelInterpreter ='latex';
%legend('dKe/dt','vedivPe','Kedivve','-qEe')
legend('$\frac{\nabla u_{e}:\overline{P}_{e}}{\Delta\varepsilon_{0}}$',...
    '$\frac{p\theta_{e}}{\Delta\varepsilon_{0}}$',...
    '$\frac{PiD_{e}}{\Delta\varepsilon_{0}}$',...
    'Interpreter','latex','FontSize',fs)
%plot(xll,ve_2_1d)
xlabel('$r / d_{i}$','Interpreter','latex','FontSize',fs)
ylabel('Energy rate terms','Interpreter','latex','FontSize',fs)
xlim([4,8])
ylim([-6.2, 6.2])
hold off

f2178 = figure(2178);
sgtitle('$$t = \ $$'+string(round(time_t,2)) + '$$ \ \omega_{i}^{-1}$$','Interpreter','latex','FontSize',20)
plot(xll,smoothdata(Pit_grad_vi_Gp_1d)./(VA2c*Pi0),'Color',colorblind(1,:), 'LineWidth', lw)
hold on
plot(xll,smoothdata(p_theta_i_Gp_1d)./(VA2c*Pi0),'Color',colorblind(2,:), 'LineWidth', lw)
plot(xll,smoothdata(PiD_i_Gp_1d)./(VA2c*Pi0),'Color',colorblind(5,:), 'LineWidth', lw)
%plot(xll,Pet_grad_ve_Gp_1d + ve_Div_Pet_1d,'g')
xline(4,'--k'); xline(8,'--k'); xline(5.7,'--b')
ax = gca; ax.FontSize = fs; ax.TickLabelInterpreter ='latex';
%legend('dKe/dt','vedivPe','Kedivve','-qEe')
legend('$\frac{\nabla u_{i}:\overline{P}_{i}}{\Delta\varepsilon_{0}}$',...
    '$\frac{p\theta_{i}}{\Delta\varepsilon_{0}}$',...
    '$\frac{PiD_{i}}{\Delta\varepsilon_{0}}$',...
    'Interpreter','latex','FontSize',fs)
%plot(xll,ve_2_1d)
xlabel('$r / d_{i}$','Interpreter','latex','FontSize',fs)
ylabel('Energy rate terms','Interpreter','latex','FontSize',fs)
xlim([4,8])
ylim([-2, 2])
hold off

% This is about my old use of ii and ij instead of the yang definition
%--------------------------------------------------------------------------
%{
f2177 = figure(2177);
sgtitle('$$t = \ $$'+string(round(time_t,2)) + '$$ \ \omega_{i}^{-1}$$','Interpreter','latex','FontSize',20)
plot(xll,smoothdata(Pet_grad_ve_Gp_1d)./(VA2c*Pi0),'Color',colorblind(1,:), 'LineWidth', lw)
hold on
plot(xll,smoothdata(Pet_grad_ve_ii_Gp_1d)./(VA2c*Pi0),'Color',colorblind(2,:), 'LineWidth', lw)
plot(xll,smoothdata(Pet_grad_ve_ij_Gp_1d)./(VA2c*Pi0),'Color',colorblind(5,:), 'LineWidth', lw)
%plot(xll,Pet_grad_ve_Gp_1d + ve_Div_Pet_1d,'g')
xline(4,'--k'); xline(8,'--k'); xline(5.7,'--b')
ax = gca; ax.FontSize = fs; ax.TickLabelInterpreter ='latex';
%legend('dKe/dt','vedivPe','Kedivve','-qEe')
legend('$\frac{\overline{P}_{e}:\nabla u_{e}}{\Delta\varepsilon_{0}}$',...
    '$\frac{p\theta_{e}}{\Delta\varepsilon_{0}}$',...
    '$\frac{PiD_{e}}{\Delta\varepsilon_{0}}$',...
    'Interpreter','latex','FontSize',fs)
%plot(xll,ve_2_1d)
xlabel('$r / d_{i}$','Interpreter','latex','FontSize',fs)
ylabel('Energy rate terms','Interpreter','latex','FontSize',fs)
xlim([4,8])
%ylim([-0.03, 0.03])
hold off
%}
%--------------------------------------------------------------------------

% things
%{
f218 = figure(218);
sgtitle('$$t = \ $$'+string(round(time_t,2)) + '$$ \ \omega_{i}^{-1}$$','Interpreter','latex','FontSize',20)
%dum1=(Pet_grad_ve_Gp_1d-Pet_grad_ve_ii_Gp_1d)./Pet_grad_ve_Gp_1d;
dum1=smoothdata((Pet_grad_ve_ij_Gp_1d)./Pet_grad_ve_Gp_1d);
semilogy(xll,10+dum1,'Color',colorblind(1,:), 'LineWidth', lw)
hold on
dum2=smoothdata((Pet_grad_ve_ii_Gp_1d)./Pet_grad_ve_Gp_1d);
semilogy(xll,10+dum2,'Color',colorblind(2,:), 'LineWidth', lw)
xline(4,'--k'); xline(8,'--k'); xline(5.7,'--b')
%legend('dKe/dt','vedivPe','Kedivve','-qEe')
ax = gca; ax.FontSize = fs; ax.TickLabelInterpreter ='latex';
legend('$10+\frac{PiD_{e}}{P_{e}:\nabla u_{e}}$',...
    '$10+\frac{p\theta_{e}}{P_{e}:\nabla u_{e}}$',...
    'Interpreter','latex','FontSize',fs)
%plot(xll,ve_2_1d)
xlabel('$r / d_{i}$','Interpreter','latex','FontSize',fs)
ylabel('Energy rate terms','Interpreter','latex','FontSize',fs)
%xlim([4,8])
ylim([0.1, 100])
hold off


f2185 = figure(2185);
sgtitle('$$t = \ $$'+string(round(time_t,2)) + '$$ \ \omega_{i}^{-1}$$','Interpreter','latex','FontSize',20)
%dum1=(Pet_grad_ve_Gp_1d-Pet_grad_ve_ii_Gp_1d)./Pet_grad_ve_Gp_1d;
dum1=smoothdata((Pet_grad_ve_ij_Gp_1d)./Pet_grad_ve_Gp_1d);
semilogy(xll,10+dum1,'Color',colorblind(1,:), 'LineWidth', lw)
hold on
dum2=smoothdata((Pet_grad_ve_ii_Gp_1d)./Pet_grad_ve_Gp_1d);
semilogy(xll,10+dum2,'Color',colorblind(2,:), 'LineWidth', lw)
xline(4,'--k'); xline(8,'--k'); xline(5.7,'--b')
%legend('dKe/dt','vedivPe','Kedivve','-qEe')
ax = gca; ax.FontSize = fs; ax.TickLabelInterpreter ='latex';
legend('$10+\frac{PiD_{e}}{P_{e}:\nabla u_{e}}$',...
    '$10+\frac{p\theta_{e}}{P_{e}:\nabla u_{e}}$',...
    'Interpreter','latex','FontSize',fs)
%plot(xll,ve_2_1d)
xlabel('$r / d_{i}$','Interpreter','latex','FontSize',fs)
ylabel('Energy rate terms','Interpreter','latex','FontSize',fs)
xlim([4,8])
ylim([0.1, 100])
hold off
%}

% This is aboutsome olf plotd of the diagonal and off diagonal terms
%--------------------------------------------------------------------------
%{
f218 = figure(218);
dum1=sign(Pet_grad_ve_ii_Gp_1d).*log(Pet_grad_ve_ii_Gp_1d./Pet_grad_ve_Gp_1d);
plot(xll,dum1,'b','LineWidth',lw)
hold on
dum2=sign(Pet_grad_ve_ij_Gp_1d).*log(Pet_grad_ve_ij_Gp_1d./Pet_grad_ve_Gp_1d);
plot(xll,dum2,'m','LineWidth',lw)
xline(4,'--k'); xline(8,'--k'); xline(5.7,'--b')
%legend('dKe/dt','vedivPe','Kedivve','-qEe')
legend('(P_{e}:\nabla u_{e})_{ii}/P_{e}:\nabla u_{e}','(P_{e}:\nabla u_{e})_{ij}/P_{e}:\nabla u_{e}')
%plot(xll,ve_2_1d)
xlabel('r / d_{i}')
ylabel('log(P_{e}:\nabla u_{e}_{ij})/P_{e}:\nabla u_{e})')
%xlim([4,8])
hold off

f219 = figure(219);
dum1=(Pet_grad_ve_ii_Gp_1d./Pet_grad_ve_Gp_1d);
semilogy(xll,dum1,'b','LineWidth',lw)
hold on
dum2=1+(Pet_grad_ve_ij_Gp_1d./Pet_grad_ve_Gp_1d);
semilogy(xll,dum2,'m','LineWidth',lw)
xline(4,'--k'); xline(8,'--k'); xline(5.7,'--b')
%legend('dKe/dt','vedivPe','Kedivve','-qEe')
legend('(P_{e}:\nabla u_{e})_{ii}/P_{e}:\nabla u_{e}','1+(P_{e}:\nabla u_{e})_{ij}/P_{e}:\nabla u_{e}')
%plot(xll,ve_2_1d)
xlabel('r / d_{i}')
ylabel('(P_{e}:\nabla u_{e}_{ij})/P_{e}:\nabla u_{e})')
%xlim([4,8])
hold off
%}
%--------------------------------------------------------------------------

%This is to save the plots
%--------------------------------------------------------------------------
%
cd '/Volumes/PSC_DiRAC_DATA/Energy_Budget_things'
saveas(f214,strcat('s_upuaur_'+ string(N_steps) +'.png'));
saveas(f2145,strcat('s_upuaur_subset_'+ string(N_steps) +'.png'));
saveas(f215,strcat('s_kinetic_thermal_'+ string(N_steps) +'.png'));
saveas(f2155,strcat('s_kinetic_thermal_subset_'+ string(N_steps) +'.png'));
saveas(f216,strcat('s_divPu_PGu_udivP_'+ string(N_steps) +'.png'));
saveas(f2165,strcat('s_divPu_PGu_udivP_subset_'+ string(N_steps) +'.png'));
saveas(f217,strcat('s_PGu_PGuii_PGuij_'+ string(N_steps) +'.png'));
saveas(f2175,strcat('s_PGu_PGuii_PGuij_subset_'+ string(N_steps) +'.png'));
saveas(f218,strcat('s_PGu_PGuij_over_PGu_10_'+ string(N_steps) +'.png'));
saveas(f2185,strcat('s_PGu_PGuij_over_PGu_10_subset'+ string(N_steps) +'.png'));

cd '/Volumes/PSC_DiRAC_DATA/Energy_Budget_things'
saveas(f27551,strcat('electron_ion_energy_terms_'+ string(N_steps) +'.png'));
saveas(f2755,strcat('electron_ion_energy_terms_short_'+ string(N_steps) +'.png'));
cd '/Users/jeffersson_agudelo/Documents/CB104_local_data/Matlab_Scripts_here';


clf(f214); clf(f2145); 
clf(f215); clf(f2155);
clf(f216); clf(f2165);
clf(f217); clf(f2175);
clf(f218); clf(f2185);
%}
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------

cd '/Users/jeffersson_agudelo/Documents/CB104_local_data/Matlab_Scripts_here'

%--------------------------------------------------------------------------