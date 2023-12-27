%% Compute files
diary logfile

% Parameters
Ncores = 12;
Ib_ini = 1;
Ib_end = 22;

% Start parallel pool
% parpool('local',Ncores)

% Solve
parfor Ib_fe = Ib_ini:Ib_end
    try
        disp(['I = ',num2str(Ib_fe),' STARTED ----------------------------------------']);
        mainfcn(Ib_fe),
        disp(['I = ',num2str(Ib_fe),' COMPLETED --------------------------------------']);
    catch
        disp(['I = ',num2str(Ib_fe),' FAILED -----------------------------------------']);
    end
end

diary off