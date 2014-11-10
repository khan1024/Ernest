% Author: Pavel Kwiecien, pavel.kwiecien@seznam.cz
% Czech Technical University in Prage, Optical Physics Group, Czech Republic

% pruchod celou strukturou

% Lambda_x, Lambda_y must be same in all grating

for i=1:1:number_of_layers
    if i==1    
        px_=eval(['coordinate_x_' num2str(i)]);
        py_=eval(['coordinate_y_' num2str(i)]);
        Lambda_x_1=px_(end)-px_(1);
        Lambda_y_1=py_(end)-py_(1);
    else        
        px_=eval(['coordinate_x_' num2str(i)]);    
        py_=eval(['coordinate_y_' num2str(i)]);    
        Lambda_x_=px_(end)-px_(1);    
        Lambda_y_=py_(end)-py_(1);    
    
        if Lambda_x_1==Lambda_x_ && Lambda_y_1==Lambda_y_
            %OK
        else
            error('periods are not identical, check your input file')
        end
    
    end
end
clear Lambda_x Lambda_x_ Lambda_y Lambda_y_

matice_stejnosti=zeros(number_of_layers,2);

%zaloha_px=px;
%zaloha_py=py;
%zaloha_n_mat=n_mat;

for i=1:1:number_of_layers
    if i==1 % prvni vrsva je vzdy uplne nova
        matice_stejnosti(i,1)=i;
        matice_stejnosti(i,2)=i;
        
    else
        matice_stejnosti(i,1)=i;
        
        index_konec=0;
        for j=1:1:i
            if j==i
                matice_stejnosti(i,2)=j;
            else
                px_i=eval(['coordinate_x_' num2str(i)]);
                px_j=eval(['coordinate_x_' num2str(j)]);
                py_i=eval(['coordinate_y_' num2str(i)]);
                py_j=eval(['coordinate_y_' num2str(j)]);
                n_mat_i=eval(['r_index_' num2str(i)]);
                n_mat_j=eval(['r_index_' num2str(j)]);
                
                if length(px_i)==length(px_j) && length(py_i)==length(py_j)
                    pomocny_soucet=sum(abs(px_j-px_i))+sum(abs(py_j-py_i))+sum(sum(abs(n_mat_j-n_mat_i)));
                else
                    pomocny_soucet=1; % they are different
                end
                
                %pomocny_soucet=sum(abs(px(:,:,j)-px(:,:,i)))+...
                %    sum(abs(py(:,:,j)-py(:,:,i)))+sum(sum(abs(n_mat(:,:,j)-n_mat(:,:,i))));
                
                if pomocny_soucet<1E-12
                    % jsou stejne
                    matice_stejnosti(i,2)=j;
                    index_konec=1; % ukoncit vnitrni for-cyklus
                else
                    % jsou jine
                end
                if index_konec==1
                    break
                end
            end
        end

    end
end

% jak nastavit novy vypocet

sort_matice_stejnosti=sort(matice_stejnosti(:,2));

index_struktura=1;
for i=1:1:number_of_layers
    if i==1
        new_structure(index_struktura,2)=sort_matice_stejnosti(1);
        new_structure(index_struktura,1)=index_struktura;
        index_struktura=index_struktura+1;
    else
        if sort_matice_stejnosti(i)==sort_matice_stejnosti(i-1)
        else
            new_structure(index_struktura,2)=sort_matice_stejnosti(i);
            new_structure(index_struktura,1)=index_struktura;
            index_struktura=index_struktura+1;
        end
    end
end

%matice_stejnosti_final[čislo vrstvy;číslo opakující předešle (současné) se vrstvy; číslo v nové struktuře]

for i=1:1:number_of_layers
    matice_stejnosti_final(i,1)=matice_stejnosti(i,1);
    matice_stejnosti_final(i,2)=matice_stejnosti(i,2);
    
    pozice_v_matici=find(new_structure(:,2)==matice_stejnosti(i,2));
    matice_stejnosti_final(i,3)=new_structure(pozice_v_matici,1);
end

% vytvoreni nove struktury
original_number_of_layers=number_of_layers;

number_of_layers=length(new_structure(:,1));
homogenous_layer=zeros(length(new_structure(:,1)),1);
current_hom_n=zeros(length(new_structure(:,1)),1);

exist_nvm=exist('nvm_method','var');

if exist_nvm==1
    if length(nvm_method)~=original_number_of_layers
        nvm_method_o(1:original_number_of_layers)=0;
        nvm_method_o(1:length(nvm_method))=nvm_method_o(1:length(nvm_method))+nvm_method;
    else
        nvm_method_o=nvm_method;
    end
else
    nvm_method_o(1:original_number_of_layers)=0;
end


for i=1:1:number_of_layers
    nvm_method_(i)=nvm_method_o(new_structure(i,2));
    
    temp_px=eval(['coordinate_x_' num2str(new_structure(i,2))]);
    eval(['new_px_' num2str(i), '=temp_px;']);
    
    temp_py=eval(['coordinate_y_' num2str(new_structure(i,2))]);
    eval(['new_py_' num2str(i), '=temp_py;']);
    
    temp_n_mat=eval(['r_index_' num2str(new_structure(i,2))]);
    eval(['new_n_mat_' num2str(i), '=temp_n_mat;']);
    
    %px(:,:,i)=px(:,:,new_structure(i,2));
    %py(:,:,i)=py(:,:,new_structure(i,2));
    %n_mat(:,:,i)=n_mat(:,:,new_structure(i,2));
    
    if abs((max(max(temp_n_mat)))-min(min(temp_n_mat)))<1E-12 % je homogenni
        homogenous_layer(i)=1; % je homogenni
        current_hom_n(i)=temp_n_mat(1,1);
    else
        homogenous_layer(i)=2; % neni homogenni
    end
end
clear temp_px temp_py temp_n_mat