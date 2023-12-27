%% Frequency plot

% Load data
clear all
load('101x101_old.mat')

% Select and plot
for i = 1:length(wi)
    % Frequency vectors
    Wi = wi{i};
    Wr = wr{i};
    
    % Selection of mode 1
    ind1 = find(Wr >3 & Wr < 8 & Wi > -0.5);
    wr1(i)  = Wr(ind1);
    wi1(i)  = Wi(ind1);
    
    % Selection of mode 2
    ind2 = find(Wr >9 & Wr < 18 & Wi > -0.5);
    wr2(i)  = Wr(ind2);
    wi2(i)  = Wi(ind2);
end

%% Represent
subplot(1,2,1)
yyaxis left
plot(I,wr1,'linewidth',1)
hold on
plot(I, 0.16*I+4.3,'--')
grid on
ylabel('Re(\omega_1)')
legend('Code','Experimental')
yyaxis right
plot(I, wi1,'linewidth',1)
ylabel('Im(\omega_1)')
xlabel('I (A)')
subplot(1,2,2)
yyaxis left
plot(I,wr2,'linewidth',1)
hold on
plot(I, 0.14*I+10,'--')
grid on
ylabel('Re(\omega_2)')
legend('Code','Experimental')
yyaxis right
plot(I, wi2,'linewidth',1)
ylabel('Im(\omega_2)')
xlabel('I (A)')

% Write data
dlmwrite('./wn_CFD.txt',[I',wr1',wr2'],',')
dlmwrite('./wi_CFD.txt',[I',wi1',wi2'],',')