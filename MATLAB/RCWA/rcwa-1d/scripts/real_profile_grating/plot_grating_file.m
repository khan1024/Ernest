figure;

for i=1:1:number_of_layers
    for j=2:2:size_result_matrix_1(2)
        if isnan(result_matrix(i,j+2))==1
            
        else
            x_begin=result_matrix(i,j);
            x_end=result_matrix(i,j+2);
            x_width=x_end-x_begin;
            y_begin=y_max-result_matrix(i,1)*i;
            y_height=result_matrix(i,1);
            barva=[abs(cos(3.5*real(result_matrix(i,j+1)))),abs(cos(5.26*real(result_matrix(i,j+1)))),abs(sin(35*real(result_matrix(i,j+1))))];
            if x_width>0
                rectangle('Position',[x_begin,y_begin,x_width,y_height],'FaceColor',barva)
            end
        end
    end    
end
hold on
for i=1:1:number_of_input_file
    plot(data(:,1,i),data(:,2,i),'LineWidth',2,'Color','r')
end
xlabel('\Lambda [\mu{}m]','FontSize',17);
ylabel('d [\mu{}m]','FontSize',17);
