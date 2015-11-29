function [ joint_prob ] = HMM_calculate_joint(data,num_states,num_discrete_obs,num_samples)

observed_sequence=HMM_discretize(aggregate(data,10),num_discrete_obs);
[B,b_1,b_inf]=HMM_calculate_params(data,num_states,num_discrete_obs,num_samples);

B_vals = B(observed_sequence);

% p is product of B_xt for t=1 to N
% where N is the number of observations

if numel(B_vals) >= 1
    p = B_vals{1};
    for i=2:numel(B_vals)
        p = p*B_vals{i};
    end
end

joint_prob=(b_inf)'*p*b_1;
