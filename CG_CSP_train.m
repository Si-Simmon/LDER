function [Q,Q_rest] = CG_CSP_train(data1,data2,m)                 %data:N¡ÁS¡ÁT
trial_1 = size(data1,3);                                                      
trial_2 = size(data2,3);
R_1 = 0; R_2 = 0;
for trial1_i = 1:trial_1;
    data_mean_1 = mean(data1(:,:,trial1_i),2)*ones(1,size(data1,2));          
    data1(:,:,trial1_i) = data1(:,:,trial1_i) - data_mean_1;
    R_1 = R_1+(data1(:,:,trial1_i)*data1(:,:,trial1_i)'/trace(data1(:,:,trial1_i)*data1(:,:,trial1_i)'));
end
R_1 = R_1/trial_1;                   
for trial2_i = 1:trial_2;
    data_mean_2 = mean(data2(:,:,trial2_i),2)*ones(1,size(data2,2));         
    data2(:,:,trial2_i) = data2(:,:,trial2_i) - data_mean_2;
    R_2 = R_2+(data2(:,:,trial2_i)*data2(:,:,trial2_i)'/trace(data2(:,:,trial2_i)*data2(:,:,trial2_i)'));
end
R_2 = R_2/trial_2;

R_sum = R_1+R_2;                        

[Uc,eigv] = eig(R_sum);

[eigv,ind] = sort(diag(eigv),'descend');

Uc = Uc(:,ind);

 W = sqrt(inv(diag(eigv))) * Uc';
 
 S_1 = W * R_1 * W';
 S_2 = W * R_2 * W';  
    if any(any(isnan(S_1))) || any(any(isnan(S_2)))
        disp('S_1 or S_2 equals to NaN£¬Jump out of the function£¡');   
        Q = zeros(2*m,size(data1,1));
        Q_rest = zeros(size(data1,1)-2*m,size(data1,1));
        return;       
    end
  
 [B_1,Psi{1}] = eig(S_1);
 [B_2,Psi{2}] = eig(S_2);


 
 [Psi{1},ind] = sort(diag(Psi{1})); B_1 = B_1(:,ind);
 [Psi{2},ind] = sort(diag(Psi{2})); B_2 = B_2(:,ind);
 B_rest = B_1(:,(m+1):1:end-m);
 Q_rest = B_rest'*W;
B_1 = B_1(:,[1:m end-m+1:end]);
 Q = B_1'*W;


