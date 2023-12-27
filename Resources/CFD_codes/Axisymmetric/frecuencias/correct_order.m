%% Correct order of EIGS
list_files = [{'91x91.mat'},{'101x101.mat'},{'111x111.mat'}];
for i = 1:length(list_files)
    load(list_files{i})
    for j = 1:23
        if wr{j}(1) > wr{j}(2)
            wrtemp   = wr{j}(1);
            wr{j}(1) = wr{j}(2);
            wr{j}(2) = wrtemp;
            witemp   = wi{j}(1);
            wi{j}(1) = wi{j}(2);
            wi{j}(2) = witemp;
        end
    end
    save(list_files{i},'I','wr','wi')
end