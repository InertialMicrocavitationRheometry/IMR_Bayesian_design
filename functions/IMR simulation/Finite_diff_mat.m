
function [Diff_Matrix ] = Finite_diff_mat(Nodes,order,Tm_check)
% Creates finite diffrence matrices
% Nodes: Number of nodes  
% order: order of differentiation 
% Tm_check: 0 not for external temp , 1 used for ext. temp 


if Tm_check == 0 
    
    N = Nodes-1; 
    deltaY = 1/N;
    K = 1:1:N+1;
    yk = (K-1)*deltaY; 
    
elseif Tm_check == 1 
    
    N = Nodes-1; 
    deltaY = -2/N;
    K = 1:1:N+1;
    yk = 1+(K-1)*deltaY; 
    
end 

Diff_Matrix = zeros(Nodes,1);

if order == 1
    
    
    
        for counter=2:(N)  %in between

            Diff_Matrix(counter,counter+1) = 1/2 ;
            Diff_Matrix(counter,counter-1) = -1/2 ;
            
        end
         
        if Tm_check == 0 
            Diff_Matrix(end,end) = 3/2; 
            Diff_Matrix(end,end-1) = -2;
            Diff_Matrix(end,end-2) = 1/2;

        end
  
        Diff_Matrix = Diff_Matrix / deltaY ; 
    
elseif order == 2 
    
      
        for counter=2:(N)  %in between

            if Tm_check == 0 
                
                Diff_Matrix(counter,counter+1) = 1+deltaY/yk(counter) ;
                Diff_Matrix(counter,counter) = -2 ;
                Diff_Matrix(counter,counter-1) = 1-deltaY/yk(counter) ;
                
            elseif Tm_check == 1
                
                Diff_Matrix(counter,counter+1) = 1 ;
                Diff_Matrix(counter,counter) = -2 ;
                Diff_Matrix(counter,counter-1) = 1 ;
            end
        end
        
            if Tm_check == 0     
                
                Diff_Matrix(1,1)= -6; Diff_Matrix(1,2) = 6;
                
            end
            
            
            Diff_Matrix = Diff_Matrix / (deltaY^2) ;
 end
    
    
end

