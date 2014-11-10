function delete=equation_line2(x_input,y_input,x_1,y_1,x_2,y_2,direction)

delete=0.*x_input+1;
%delete=1; % yes

if x_1==x_2
    % direction, 1-delete left cell, 2-right
    for i=1:1:length(delete)    
        if direction==2        
            if x_input(i)<x_1            
                delete(i)=2;        
            end            
        end
        
        if direction==1        
            if x_input(i)>=x_1            
                delete(i)=2;        
            end            
        end
    end    
else

% direction, 1-above, 2-below

a=(y_2-y_1)/(x_2-x_1);
b=y_1-a*x_1;

y_line=a*x_input+b;

for i=1:1:length(delete)
    if direction==1    
        if y_input(i)<=y_line        
            delete(i)=2;    
        end        
    elseif direction==2
        if y_input(i)>=y_line
            delete(i)=2;
        end
    end
end

end