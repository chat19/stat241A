function [predicted_levels] = HMM_predict(data,train_len,num_states,num_discrete_obs,num_samples)

assert(length(data) > train_len);

[B,b_1,b_inf]=HMM_calculate_params(data(1:train_len),num_states,num_discrete_obs,num_samples)

% create a table of length one tenth of the original dataset
% since we are looking at returns over 10 day intervals
prediction_length=floor(length(data)/10);
predicted_levels=ones(prediction_length);

predicted_levels(1:floor(train_len/10))=data(1:floor(train_len/10));

for i = (floor(train_len/10)+1):prediction_length
    %calculate conditional probabilities

%finish this
    

