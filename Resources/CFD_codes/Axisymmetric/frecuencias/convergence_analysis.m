%% Convergence table

list_files = [{'91x91.mat'},{'101x101.mat'},{'111x111.mat'}];
for i = 1:length(list_files)
    load(list_files{i})
    Wr(:,i) = wr';
    Wi(:,i) = wi';
end
for I = 1:5:21
    T((I-1)/5+1,:)  = [Wr{I,1}(1), Wr{I,2}(1), Wr{I,3}(1), Wi{I,1}(1), Wi{I,2}(1), Wi{I,3}(1)];
    TP((I-1)/5+1,:) = [(Wr{I,1}(1)-Wr{I,2}(1))/Wr{I,2}(1)*100, 0,(Wr{I,3}(1)-Wr{I,2}(1))/Wr{I,2}(1)*100,...
        (Wi{I,1}(1)-Wi{I,2}(1))/Wi{I,2}(1)*100, 0,(Wi{I,3}(1)-Wi{I,2}(1))/Wi{I,2}(1)*100];
    
    T2((I-1)/5+1,:)  = [Wr{I,1}(2), Wr{I,2}(2), Wr{I,3}(2), Wi{I,1}(2), Wi{I,2}(2), Wi{I,3}(2)];
    TP2((I-1)/5+1,:) = [(Wr{I,1}(2)-Wr{I,2}(2))/Wr{I,2}(2)*100, 0,(Wr{I,3}(2)-Wr{I,2}(2))/Wr{I,2}(2)*100,...
        (Wi{I,1}(2)-Wi{I,2}(2))/Wi{I,2}(2)*100, 0,(Wi{I,3}(2)-Wi{I,2}(2))/Wi{I,2}(2)*100];
end