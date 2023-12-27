%% Expansion of xo
% We pass from a single unknown vector _x0_ to "readable" variables
% just for inspection/drawing purposes
%
nbi = 0;
for j=1:nbl
    bl = list_block{j};
    NVAR = eval(['NV' bl]);
    nunkn = eval(['nt' bl]); %number of gridpoints in block
    for i=1:NVAR
        inic = nbi+(i-1)*nunkn+1;
        ifin = nbi+i*nunkn;
        lv = eval(['list_var_' bl '{' num2str(i) '}']);
        order = [lv '= full(reshape(x0(' num2str(inic) ':' num2str(ifin) '),nr' bl ',nz' bl '));'];
      
        evalc(order);
    end
    nbi = nbi + nunkn*NVAR; 
end
%% Parameters are set in a readeable way

 


  
  
  
  
  
  
  
  
  
  
      % We construct an unique vector of unknowns _x0_
    %
    

  
    
    order='[';
    for j=1:nbl
        bl = list_block{j};
        NVAR = eval(['NV' bl]);
        for i=1:NVAR
            lv = eval(['list_var_' bl '{' num2str(i) '}']);
            order = [order 'reshape(' lv ',nt' bl ',1);'];
        end
    end
    x0=eval([order ']']);
