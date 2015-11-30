function [predicted_obs, true_obs] = HMM_predict(data,subset_size,num_states,num_discrete_obs,num_samples)
%data=data
%subset=first (subset) days of data we use for training
assert(length(data) > subset_size);

[B,b_1,b_inf]=HMM_calculate_params(data(1:subset_size),num_states,num_discrete_obs,num_samples);

% create a table of length one tenth of the original dataset
% since we are looking at returns over 10 day intervals

num_days=10;

%now everything is going to have discrete observations

training_length = floor(subset_size/num_days);
prediction_length=floor(length(data)/num_days);
predicted_obs=zeros(1,prediction_length)';

%timestep = (1:prediction_length)';

%add our training data to our prediction table
true_obs = HMM_discretize(aggregate(data,num_days),num_discrete_obs);
train_obs = HMM_discretize(aggregate(data(1:subset_size),num_days),num_discrete_obs);
predicted_obs(1:training_length)=train_obs;

%b represents internal states
b=cell(1,prediction_length)';
b{1}=b_1;
for i=1:(training_length-1)
    xi = predicted_obs(i);
    b{i+1} = (B{xi}* b{i})/((b_inf)'*B{xi}*b{i});
end

%use this training data to predict the rest
for i = training_length:(prediction_length-1)
    %internal state update
    xi = predicted_obs(i);
    b{i+1} = (B{xi}* b{i})/((b_inf)'*B{xi}*b{i});
    
    %find the highest probability next observation state
    %calculate conditional probabilities
    cond_probs = cellfun(@(y) (b_inf)'*y*b{i} / sum(cellfun(@(x) (b_inf)'*x*b{i},B)),B);
    [~,predicted_observation]=max(cond_probs);
    predicted_obs(i+1)=predicted_observation;
    %commented out because this is vectorized above
    %for j=1:num_discrete_obs
        %cond_probs{j}=( (b_inf)'*B{j}*b{i} ) / sum(cellfun(@(x) (b_inf)'*x*b{i},B));
    %end
end


%finish this
    

