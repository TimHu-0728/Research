%% Numerical computation of variables and derivatives
  
%% Computation of the time derivative
%
% The generic time integration is a second order backward, 
%
% $$ \dot{y}(t^n) = -\frac{2 \Delta t_1 +\Delta t_2}{\Delta t_1(\Delta t_1
% +\Delta t_2)} y^n + \frac{\Delta t_1 +\Delta t_2}{\Delta t_1 \Delta t_1
% } y^{n-1} + \frac{\Delta t_1}{\Delta t_2 (\Delta t_1
% +\Delta t_2)} y^{n-2}  $$
%
% in other words,  
%
% $$ \dot{y}(t^n) = b_p \,  y^{n-1} + b_m \, y^{n-1} + b_{mm} \, y^{n-2}  $$
% 
list_corner_z = {'1' 'nz'};
list_corner_r = {'1' 'nr'};
list_pref= {'l' 'r'};

%fileID = fopen([path_jacobian 'matrix.m'],'w');
fileID = fopen('matrix.m','w');

fprintf(fileID,'%%% Matrix generation');
fprintf(fileID,'%');
fprintf(fileID,'%% Time derivative');
fprintf(fileID,' \n dt2=dt1+dt; \n bm= -(dt2/dt)/(dt2-dt); \n \n');
fprintf(fileID,' \n  bmm= (dt/dt2)/(dt2-dt); \n bp=-((dt2/dt)^2-1)/((dt2/dt)*(dt-dt2)); \n \n');
fprintf(fileID,' \n xt = bp*x0 + bm*x0m + bmm*x0mm;\n \n');

%% Initializing the variables 
%
%From the entire array _x0_ we have to extract the specific unknown (for
%example the axial velocity of block B, _wB_) compute its derivatives and
%asign them to the symbolic vectors). To do this we have to take care 
%of the order of the derivatives set in the corresponding block.
%

fprintf(fileID,' %% Initializing the variables \n \n');
fprintf(fileID,'nbi = 0;\n'); % #num block inic
for j=1:nbl
    bl = list_block{j};
    NVAR = ['NV' bl];
    NDER = ['ND' bl];
    NT = ['nt' bl]; %number of total unknown in block
    order = ['yf' bl  '= zeros(' NVAR ',' NDER ',' NT ');'];
    fprintf(fileID,'%s\n\n',order);
    for i=1:eval(NVAR)
        fprintf(fileID,'inic = nbi+%s*%s+1; \n',num2str(i-1),NT); 
        fprintf(fileID,'ifin = nbi+%s*%s;\n',num2str(i),NT); 
        
        lv = eval(['list_var_' bl '{' num2str(i) '}']);
        order = [lv '= x0(inic:ifin);'];
        fprintf(fileID,'%s\n',order);
        %%
        % the variable "per se"
        order = [ 'yf' bl '(' num2str(i) ',1,:)=' lv ';'];
        fprintf(fileID,'%s\n',order);
        for k = 2:eval(NDER)-1
            order = ['list_der' '{' num2str(k) '}'];
            deriv = [eval(order) bl];
            order = ['yf' bl, '(' num2str(i) ',' num2str(k) ',: ) = dd' deriv '*' lv ';'];
            fprintf(fileID,'%s\n',order);
        end
        order = ['yf' bl, '(' num2str(i) ',' num2str(eval(NDER)) ',: ) = xt(inic:ifin);'];
         fprintf(fileID,'%s\n \n \n',order); 
    end
    fprintf(fileID,'nbi = nbi + %s*%s; \n \n', NT,NVAR); 
end

%% Calculation of numerical functions and Jacobians
%
%% Evaluation of the numerical Jacobian
%
% Note that:
%
% * i=1 (z=0) correspond to left
% * i=NZ (z=H) correspond to right
% * j=1 (r=0) correspond to bottom
% * j=NR (z=H) correspond to top
% * Out of the extremes correspond to the bulk


fprintf(fileID,'\n%%\n');
fprintf(fileID,'%% Initializing the jacobian submatrix \n');
fprintf(fileID,'%%\n \n');

for ib = 1:nbl
    bli = lower(list_block{ib});
    for jb = 1:nbl
    blj = lower(list_block{jb});
        order1 = [bli blj ' = 0*speye(nt' upper(bli) '*NV' upper(bli) ',nt' upper(blj) '*NV' upper(blj) ');'];
        fprintf(fileID, '%s \n', order1);
    end
end

 for k=1:nbl %for each block
    bl = list_block{k};
    fprintf(fileID,'\n%% Block %s\n \n',bl);
    NVAR = ['NV' bl];
    NDER = ['ND' bl];
    NZ = ['nz' bl];
    NR = ['nr' bl];
    NT = ['nt' bl]; %number of grid points in block
    NVyD = [NDER '*' NVAR]; %extended number of uknowns
    order = ['F' bl bl '= zeros(' NVAR ',' NT ');']; 
    fprintf(fileID,'%s\n',order);
    order = ['DF' bl bl '= zeros(' NVAR ',' NVyD ',' NT ');'];
    fprintf(fileID,'%s\n',order);
    fprintf(fileID,'\n %% bulk \n');
%     % *bulk*
     fprintf(fileID,'for i=2:%s-1 \n',NZ); 
          fprintf(fileID,'for j=2:%s-1 \n',NR); 
              fprintf(fileID,'l=sub2ind([%s,%s],j,i);\n', NR, NZ);
              order = ['x = reshape(yf' bl '(:,:,l)' '''' ',' NVyD ',1);'];
              fprintf(fileID,'%s\n',order);
              equat =['equationF' bl bl ...
                      '(z0' bl '(i)'...
                     ', r0' bl '(j),x,pa);'];
              order = ['[F' bl bl '(:,l),DF' bl bl '(:,:,l)]=' equat];
              fprintf(fileID,'%s\n',order);
           fprintf(fileID,'end\n'); 
              fprintf(fileID,'\n%%\n');
              fprintf(fileID,'%% *top (j=NR)* and *bottom(j=1)* \n');
              fprintf(fileID,'%%\n \n');
              fprintf(fileID,'l=sub2ind([%s,%s],1,i); %% bottom\n',NR,NZ); 
              order = ['x = reshape(yf' bl '(:,:,l)' '''' ',' NVyD ',1);'];
              fprintf(fileID,'%s\n',order);
              equat =['equationF' bl bl ...
                 'b' '(z0' bl '(i), r0' bl '(1),x,pa);'];
              order = ['[F' bl bl '(:,l),DF' bl bl '(:,:,l)]=' equat];
              fprintf(fileID,'%s\n\n',order);
              %
              % The top and bottom row of points
              %
 
              fprintf(fileID,'l=sub2ind([%s,%s],%s,i); %%top \n',NR,NZ,NR); 
              order = [' x = reshape(yf' bl '(:,:,l)' '''' ',' NVyD ',1);'];
              fprintf(fileID,'%s\n',order);
              equat =['equationF' bl bl ...
              't'  '(z0' bl '(i)'...
               ', r0' bl '(' NR '),x,pa);'];
              order = ['[F' bl bl '(:,l),DF' bl bl '(:,:,l)]=' equat];
              fprintf(fileID,'%s\n',order);
      fprintf(fileID,'end\n');

%     %%
%     % *Left (i=1)* and *right(i=NZ)*
%     %      
      fprintf(fileID,'\n%%\n');
      fprintf(fileID,'%% *Left (i=1)* and *right(i=NZ)* \n');
      fprintf(fileID,'%%\n \n');
      fprintf(fileID,'for j=2:%s-1 \n',NR); 
%         % *Left (i=1)*
          fprintf(fileID,'l = sub2ind([%s,%s],j,1); \n',NR,NZ); 
          order = ['x = reshape(yf' bl '(:,:,l)' '''' ',' NVyD ',1);'];
          fprintf(fileID,'%s\n',order);
          equat =['equationF' bl bl ...
             'l' '(z0' bl '(1), r0' bl '(j),x,pa);'];
         order = ['[F' bl bl '(:,l),DF' bl bl '(:,:,l)]=' equat];
         fprintf(fileID,'%s\n \n',order);
%        %%
%        % *Right (i=NZ)*
          fprintf(fileID,'l = sub2ind([%s,%s],j,%s); \n',NR,NZ,NZ); 
        order = ['x = reshape(yf' bl '(:,:,l)' '''' ',' NVyD ',1);'];
        fprintf(fileID,'%s\n',order);

        equat =['equationF' bl bl ...
            'r' '(z0' bl '(' NZ ')'...
               ', r0' bl '(j),x,pa);'];
        order = ['[F' bl bl '(:,l),DF' bl bl '(:,:,l)]=' equat];
        fprintf(fileID,'%s\n',order);
      fprintf(fileID,'end\n');
      
% 
%    Corners
%
%
% preference at corners: corners belong to left and right side by default

      fprintf(fileID,'\n%%corners (Must match with those in connections!) \n\n');
      for ii=1:2
          ii_corner = list_corner_z{ii};
          if(ii>1)
          ii_corner = [ii_corner bl];
          end
          pref = list_pref{ii};
          for jj =1:2
            jj_corner = list_corner_r{jj};
             if(jj>1)
            jj_corner = [jj_corner bl];
             end
            fprintf(fileID,'\n %% corner (z) i=%s and (r) j=%s \n\n',jj_corner,ii_corner);
            fprintf(fileID,'l = sub2ind([%s,%s],%s,%s); \n',NR,NZ,jj_corner,ii_corner); 
            order = ['x = reshape(yf' bl '(:,:,l)' '''' ',' NVyD ',1);'];
            fprintf(fileID,'%s\n',order);

            equat =['equationF' bl bl ...
               pref '(z0' bl '(' ii_corner ')'...
               ', r0' bl '(' jj_corner '),x,pa);'];
               order = ['[F' bl bl '(:,l),DF' bl bl '(:,:,l)]=' equat];
               fprintf(fileID,'%s\n',order);
          end
      end
        
      fprintf(fileID,'\n%%\n');
      fprintf(fileID,'%% Mounting the jacobians \n');
      fprintf(fileID,'%%\n \n');
%     
%      Mounting the jacobians and b
%     
     order = ['x' lower(bl) ' = 0*speye(' NT '*' NVAR ',1);'];
     fprintf(fileID, '%s \n\n', order);
      fprintf(fileID,'for jj=1:%s \n',NVAR); 
          fprintf(fileID,'rowi = (jj-1)*%s + 1;\n',NT);
          fprintf(fileID,'rowf = jj*%s;\n',NT);
          fprintf(fileID,'for kk=1:%s \n',NVAR); 

            fprintf(fileID,'km=(kk-1)*%s+1;\n',NDER);
            fprintf(fileID,'kp=kk*%s;\n',NDER);
            order = ['D=squeeze(DF' bl bl '(jj,km:kp,:));'];
            fprintf(fileID,'%s;\n',order);
            order =[ 'B=spdiags(D(1,:)'',0,nt' bl ', nt' bl ') ...'];
            fprintf(fileID,'%s\n',order);
            for jd=2:eval(NDER)-1
                deriv = ['list_der' '{' num2str(jd) '}'];
                deriv = [eval(deriv) bl];
                order =  ['+ spdiags(D(' num2str(jd) ',:)' '''' ',0, nt' bl ', nt' bl ')* dd' deriv '...'];
            fprintf(fileID,'%s\n',order);
            end
            order =['+ spdiags(D(' num2str(eval(NDER)) ',:)' '''' ',0, nt' bl ', nt' bl ')* bp;'];
            fprintf(fileID,'%s\n',order);
            fprintf(fileID,'coli = (kk-1)*%s + 1;\n',NT);
            fprintf(fileID,'colf = kk*%s;\n',NT);
            order = [lower(bl) lower(bl) '(rowi:rowf,coli:colf) = sparse(B);'];
            fprintf(fileID,'%s\n',order);            
          fprintf(fileID,'end\n');
          order = ['x' lower(bl) '(rowi:rowf) = - F' bl bl '(jj,:)'';'];
          fprintf(fileID,'%s\n',order);
      fprintf(fileID,'end\n');
%     
      fprintf(fileID,'\n%%\n');
      fprintf(fileID,'%% Connections \n');
      fprintf(fileID,'%%\n \n');

     order = ['list_con_' bl];
     [basura, ncon] = size(eval(order));
%     
     for nc=1:ncon
         tco = eval(['list_con_' bl '{' num2str(nc) '}']);
%        %% Geometrical information to the connected block
         % tco(1) this the block to connect with!
         fprintf(fileID,'%% Connection of %s with %s \n', bl, tco(1));
         fprintf(fileID,'%%\n \n');
         NVAR_c = ['NV' tco(1)];
         NDER_c = ['ND' tco(1)];
         NZ_c = ['nz' tco(1)];
         NR_c = ['nr' tco(1)];
         NT_c = ['nt' tco(1)]; %number of grid points in block
         NVyD_c = [NDER_c '*' NVAR_c]; %extended number of uknowns
%
         order = ['DF' bl tco(1) '= zeros(' NVAR ',' NVyD_c ',' NT_c ');'];
         fprintf(fileID,'%s\n',order);
%         %%
%         % tco(2) indicates where to connect
%         %
         switch tco(2);
             case {'t', 'b'}
                 if(tco(2) == 't')
                     j = NR;
                     jc = '1';
                 else
                     j = '1';
                     jc = NR_c;
                 end
%                 
                 order = ['F' bl tco(1) '= zeros(' NVAR ',' NZ ');'];
                 fprintf(fileID,'%s\n',order);

                 fprintf(fileID,'jl = [];\n');
                 fprintf(fileID,'jl_c = [];\n');
                 fprintf(fileID,'il = [];\n');

                 fprintf(fileID,'for i=2:%s-1\n',NZ);
                     fprintf(fileID,'ic = i;\n');
                     fprintf(fileID,'l = sub2ind([%s,%s],%s,i);\n', NR,NZ,j);
                     fprintf(fileID,'l_c = sub2ind([%s,%s],%s,ic);\n', NR_c,NZ_c,jc);
                     fprintf(fileID,'jl = [jl l];\n');
                     fprintf(fileID,'jl_c = [jl_c l_c];\n');
                     order = [' x = reshape(yf' tco(1) '(:,:,l_c)' '''' ',' NVyD_c ',1);'];
                     fprintf(fileID,'%s\n',order);
                     equat =['equationF' bl tco ...
                         '(z0' tco(1) '(ic)'...
                         ', r0' tco(1) '(' jc '),x,pa);'];
                     order = ['[F' bl tco(1) '(:,i),DF' bl tco(1) '(:,:,l_c)]=' equat];
                     fprintf(fileID,'%s\n',order);                                          
                     fprintf(fileID,'il = [il i];\n');
                     fprintf(fileID,'end\n');

          for ii=1:2
             ii_corner = list_corner_z{ii};
             if(ii>1)
                 ii_corner = [ii_corner tco(1)];
             end
             fprintf(fileID,'\n %% corners CHECK COMPABILILITY WITH CORNERS ABOVE  \n\n');            
             fprintf(fileID,'\n %% corner (z) i=%s and (r) j=%s \n\n',ii_corner,j);

             fprintf(fileID,'%% l = sub2ind([%s,%s],%s,%s);\n', NR,NZ,j,ii_corner);
             fprintf(fileID,'%% l_c = sub2ind([%s,%s],%s,%s);\n\n', NR_c,NZ_c,jc,ii_corner);
             fprintf(fileID,'%% jl = [l jl];\n');
             fprintf(fileID,'%% jl_c = [l_c jl_c];\n');
             fprintf(fileID,'%% il = [il i];\n');

             order = [' x = reshape(yf' tco(1) '(:,:,l_c)' '''' ',' NVyD_c ',1);'];
             fprintf(fileID,'%% %s\n',order);
             equat =['equationF' bl tco ...
                     '(z0' tco(1) '(' ii_corner ')'...
                     ', r0' tco(1) '(' jc '),x,pa);'];
             order = ['[F' bl tco(1) '(:,' ii_corner '),DF' bl tco(1) '(:,:,l_c)]=' equat];
             fprintf(fileID,'%% %s\n',order);
          end                   
             case {'l','r'}
                 if( tco(2) == 'l')
                     i = '1';
                     ic = NZ_c;
                 else
                     i = NZ;
                     ic = '1';
                 end

                 order = ['F' bl tco(1) '= zeros(' NVAR ',' NR ');'];
                 fprintf(fileID,'%s\n',order);
                 
                 fprintf(fileID,'jl = [];\n');
                 fprintf(fileID,'jl_c = [];\n');
                 fprintf(fileID,'il = [];\n');
                 
                 fprintf(fileID,'for j=2:%s-1\n',NR);
                     fprintf(fileID,'jc = j;\n');
                     fprintf(fileID,'l = sub2ind([%s,%s],j,%s);\n',NR,NZ,i);
                     fprintf(fileID,'l_c = sub2ind([%s,%s],j,%s);\n',NR_c,NZ_c,ic);
                     fprintf(fileID,'jl = [jl l];\n');
                     fprintf(fileID,'jl_c = [jl_c l_c];\n');
                     order = ['x = reshape(yf' tco(1) '(:,:,l_c)' '''' ',' NVyD_c ',1);'];
                     fprintf(fileID,'%s\n',order);
                     equat =['equationF' bl tco ...
                         '(z0' tco(1) '(' ic ')'...
                         ', r0' tco(1) '(jc),x,pa);'];
                     order = ['[F' bl tco(1) '(:,j),DF' bl tco(1) '(:,:,l_c)]=' equat];
                     fprintf(fileID,'%s\n',order);
                     fprintf(fileID,'il = [il i];\n');
                  fprintf(fileID,'end\n\n');
                  
         for jj=1:2
             jj_corner = list_corner_r{jj};
             if(ii>1)
                 jj_corner = [jj_corner bl];
             end
             fprintf(fileID,'\n %% corners CHECK COMPABILILITY WITH CORNERS ABOVE  \n\n');            
             fprintf(fileID,'\n %% corner (z) i=%s and (r) j=%s \n\n',i,jj_corner);

             fprintf(fileID,'%% l = sub2ind([%s,%s],%s,%s);\n', NR,NZ,jj_corner,i);
             fprintf(fileID,'%% l_c = sub2ind([%s,%s],%s,%s);\n', NR_c,NZ_c,jj_corner,ic);
             fprintf(fileID,'%% jl = [jl l];\n');
             fprintf(fileID,'%% jl_c = [jl_c l_c];\n');
             fprintf(fileID,'%% il = [il i];\n');   
             
             order = [' x = reshape(yf' tco(1) '(:,:,l_c)' '''' ',' NVyD_c ',1);'];
             fprintf(fileID,'%% %s\n',order);
             equat =['equationF' bl tco ...
                     '(z0' tco(1) '(' ii_corner ')'...
                     ', r0' tco(1) '(' jc '),x,pa);'];
             order = ['[F' bl tco(1) '(:,i),DF' bl tco(1) '(:,:,l_c)]=' equat];
             fprintf(fileID,'%% %s\n',order);
          end                                     
         end
         
      fprintf(fileID,'\n%%\n');
      fprintf(fileID,'%% Constructing %s%s \n', lower(bl), lower(tco(1)));
      fprintf(fileID,'%%\n \n');


         fprintf(fileID,'for jj=1:%s\n',NVAR);
             fprintf(fileID,'list = jl +(jj-1)*%s;\n',NT);
             fprintf(fileID,'for kk=1:%s\n',NVAR_c);
                 fprintf(fileID,'km=(kk-1)*%s+1;\n',NDER_c);
                 fprintf(fileID,'kp=kk*%s;\n',NDER_c);
                 order = ['D=squeeze(DF' bl tco(1) '(jj,km:kp,:));'];
                 fprintf(fileID,'%s\n',order);
                 order =[ 'B = spdiags(D(1,:)'',0,nt' tco(1) ', nt' tco(1) ') ...'];
                 fprintf(fileID,'%s\n',order);
                 for jd=2:eval(NDER_c)-1
                     deriv = ['list_der' '{' num2str(jd) '}'];
                     deriv = [eval(deriv) tco(1)];
                     order =['+ spdiags(D(' num2str(jd) ',:)' '''' ',0, nt' tco(1) ', nt' tco(1) ')* dd' deriv '...'];
                     fprintf(fileID,'%s\n',order);
                 end
                 order =[ '+ spdiags(D(' num2str(eval(NDER_c)) ',:)' '''' ',0, nt' tco(1) ', nt' tco(1) ')* bp;'];
                 fprintf(fileID,'%s\n',order);
                 fprintf(fileID,'aux = sparse(B);\n');
                 fprintf(fileID,'coli = (kk-1)*%s + 1;\n',NT_c);
                 fprintf(fileID,'colf = kk*%s;\n',NT_c);
                 order = [lower(bl) lower(tco(1)) '(list,coli:colf) = aux(jl_c,:);'];
                 fprintf(fileID,'%s\n',order);
             fprintf(fileID,'end\n');
             order = ['x' lower(bl) '(list) = x' lower(bl) '(list) - F' bl tco(1) '(jj,il)'';'];
             fprintf(fileID,'%s\n',order);
             
         fprintf(fileID,'end\n \n'); 
     end
 end


fprintf(fileID, '%% \n');
fprintf(fileID, '%% Mounting of vector b and off diagonal jacobian\n');
fprintf(fileID, '%% \n\n');

order ='[';
order1 = '[';
for ib = 1:nbl
    bli = lower(list_block{ib});
    order1 = [ order1 ' x' bli ';'];
    for jb = 1:nbl
    blj = lower(list_block{jb});
    order = [ order ' ' bli blj];
    end
    order = [ order ' ; ' ];
end
fprintf(fileID, '\na = %s ];\n', order);
fprintf(fileID, 'b = %s ];\n', order1);
fclose(fileID);
