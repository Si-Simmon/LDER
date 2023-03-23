function [label_generation_new] = CVgeneration(alphabet,prob,number)
%% Cross Validation generation  %%
% date: 2023.03.08

% eg. alphabet = [1 2 3 4 5]; prob = [0.7 0.075 0.075 0.075 0.075];

label_generation = randsrc(number,1,[alphabet;prob]);
label_generation_new = label_generation; 

end



