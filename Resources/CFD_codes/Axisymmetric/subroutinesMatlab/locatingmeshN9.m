la=zeros(nzb*nrb,8);

for j=1:nzb
        jmm=j-2;
        jm=j-1;
        jp=j+1;
        jpp=j+2;

    for i=1:nrb
        imm=i-2;
        im=i-1;
        ip=i+1;
        ipp=i+2;
        l=sub2ind([nrb,nzb,1],i,j,1);
        
        
    if ((i>1)&&(i<nrb) &&( j>1) && (j<nzb))   
        la(l,1)=sub2ind([nrb,nzb,1],i,j,1);
        la(l,2)=sub2ind([nrb,nzb,1],i,jm,1);
        la(l,3)=sub2ind([nrb,nzb,1],i,jp,1);
        
        la(l,4)=sub2ind([nrb,nzb,1],im,jm,1);
        la(l,5)=sub2ind([nrb,nzb,1],im,j,1);
        la(l,6)=sub2ind([nrb,nzb,1],im,jp,1);
       
        la(l,7)=sub2ind([nrb,nzb,1],ip,jm,1);
        la(l,8)=sub2ind([nrb,nzb,1],ip,j,1);
        la(l,9)=sub2ind([nrb,nzb,1],ip,jp,1);
       
        
        
        end
    if (i==1)&&(j>1)&&(j<nzb)
     
         
        
        
       
        la(l,1)=sub2ind([nrb,nzb,1],i,j,1);
        la(l,2)=sub2ind([nrb,nzb,1],i,jm,1);
        la(l,3)=sub2ind([nrb,nzb,1],i,jp,1);
        
        la(l,4)=sub2ind([nrb,nzb,1],ip,j,1);
        la(l,5)=sub2ind([nrb,nzb,1],ip,jm,1);
        la(l,6)=sub2ind([nrb,nzb,1],ip,jp,1);
        
        la(l,7)=sub2ind([nrb,nzb,1],ipp,j,1);
        la(l,8)=sub2ind([nrb,nzb,1],ipp,jm,1);
        la(l,9)=sub2ind([nrb,nzb,1],ipp,jp,1);
        
      
        la(l,8)=sub2ind([nrb,nzb,1],ipp+1,j,1);
        
    end   
    
      if ((i==nrb)&&(j>1)&&(j<nzb))
 
          la(l,1)=sub2ind([nrb,nzb,1],i,j,1);
          la(l,2)=sub2ind([nrb,nzb,1],i,jm,1);
          la(l,3)=sub2ind([nrb,nzb,1],i,jp,1);
         
          la(l,4)=sub2ind([nrb,nzb,1],im,j,1);
          la(l,5)=sub2ind([nrb,nzb,1],im,jm,1);
          la(l,6)=sub2ind([nrb,nzb,1],im,jp,1);
        
        la(l,7)=sub2ind([nrb,nzb,1],imm,j,1);
        la(l,8)=sub2ind([nrb,nzb,1],imm,jm,1);
        la(l,9)=sub2ind([nrb,nzb,1],imm,jp,1);
       % la(l,8)=sub2ind([nrb,nzb,1],imm-1,j,1);
       end 
     if (j==1)&&(i>1)&&(i<nrb)
        la(l,1)=sub2ind([nrb,nzb,1],i,j,1);
        la(l,2)=sub2ind([nrb,nzb,1],im,j,1);
        la(l,3)=sub2ind([nrb,nzb,1],ip,j,1);
        
        la(l,4)=sub2ind([nrb,nzb,1],i,jp,1);
        la(l,5)=sub2ind([nrb,nzb,1],im,jp,1);
        la(l,6)=sub2ind([nrb,nzb,1],ip,jp,1);
       
        la(l,7)=sub2ind([nrb,nzb,1],i,jpp,1);
        la(l,8)=sub2ind([nrb,nzb,1],im,jpp,1);
        la(l,9)=sub2ind([nrb,nzb,1],ip,jpp,1);
       
     end
    
        if (j==nzb)&&(i>1)&&(i<nrb)
        
        la(l,1)=sub2ind([nrb,nzb,1],i,j,1);    
        la(l,2)=sub2ind([nrb,nzb,1],im,j,1);
        la(l,3)=sub2ind([nrb,nzb,1],ip,j,1);
        
        la(l,4)=sub2ind([nrb,nzb,1],i,jm,1);
        la(l,5)=sub2ind([nrb,nzb,1],im,jm,1);
        la(l,6)=sub2ind([nrb,nzb,1],ip,jm,1);
       
        la(l,7)=sub2ind([nrb,nzb,1],i,jmm,1);
        la(l,8)=sub2ind([nrb,nzb,1],im,jmm,1);
        la(l,9)=sub2ind([nrb,nzb,1],ip,jmm,1);
        end  
   
       if ((i==1)&& (j==1))
         la(l,1)=sub2ind([nrb,nzb,1],i,j,1);   
         la(l,2)=sub2ind([nrb,nzb,1],i,jp,1);
         la(l,3)=sub2ind([nrb,nzb,1],i,jpp,1);
        
        la(l,4)=sub2ind([nrb,nzb,1],ip,j,1);
        la(l,5)=sub2ind([nrb,nzb,1],ip,jp,1);
        la(l,6)=sub2ind([nrb,nzb,1],ip,jpp,1); 
        
        la(l,7)=sub2ind([nrb,nzb,1],ipp,j,1);  
        la(l,8)=sub2ind([nrb,nzb,1],ipp,jp,1);  
        la(l,9)=sub2ind([nrb,nzb,1],ipp,jpp,1); 
       end
    
           if ((i==1)&& (j==nzb))
         la(l,1)=sub2ind([nrb,nzb,1],i,j,1);   
         la(l,2)=sub2ind([nrb,nzb,1],i,jm,1);
         la(l,3)=sub2ind([nrb,nzb,1],i,jmm,1);
        
        la(l,4)=sub2ind([nrb,nzb,1],ip,j,1);
        la(l,5)=sub2ind([nrb,nzb,1],ip,jm,1);
        la(l,6)=sub2ind([nrb,nzb,1],ip,jmm,1);
        
        la(l,7)=sub2ind([nrb,nzb,1],ipp,j,1);  
        la(l,8)=sub2ind([nrb,nzb,1],ipp,jm,1);  
        la(l,9)=sub2ind([nrb,nzb,1],ipp,jmm,1); 
           end
     
     if ((i==nrb)&& (j==1))
         la(l,1)=sub2ind([nrb,nzb,1],i,j,1);   
         la(l,2)=sub2ind([nrb,nzb,1],i,jp,1);
         la(l,3)=sub2ind([nrb,nzb,1],i,jpp,1);
         
        la(l,4)=sub2ind([nrb,nzb,1],im,j,1);
        la(l,5)=sub2ind([nrb,nzb,1],im,jp,1);
        la(l,6)=sub2ind([nrb,nzb,1],im,jpp,1);
        
       la(l,7)=sub2ind([nrb,nzb,1],imm,j,1); 
       la(l,8)=sub2ind([nrb,nzb,1],imm,jp,1); 
       la(l,9)=sub2ind([nrb,nzb,1],imm,jpp,1);   
       end
       
           if ((i==nrb)&& (j==nzb))
         la(l,1)=sub2ind([nrb,nzb,1],i,j,1);   
         la(l,2)=sub2ind([nrb,nzb,1],i,jm,1);
         la(l,3)=sub2ind([nrb,nzb,1],i,jmm,1);
        la(l,4)=sub2ind([nrb,nzb,1],im,j,1);
        la(l,5)=sub2ind([nrb,nzb,1],im,jm,1);
        la(l,6)=sub2ind([nrb,nzb,1],im,jmm,1);   
        la(l,7)=sub2ind([nrb,nzb,1],imm,j,1);
        la(l,8)=sub2ind([nrb,nzb,1],imm,jm,1);
        la(l,9)=sub2ind([nrb,nzb,1],imm,jmm,1);
           end
       end
       
       
       
       
       
end