function [ ranked_obs ] = HMM_discretize(observations,num_obs_states)
%rank -- takes observed values and num_states as an array and integer
% and gives each element a rank corresponding to the quantile

%do this to transform our observations into those suitable for HMM training
levels=quantile(observations,num_obs_states-1);
ranked_obs = ones(1,length(observations))';
%initialize ranks to 1

for i=1:length(observations)
    %find corresponding rank for each index in observations
    for j=1:num_obs_states-1
        if (observations(i)>=levels(j))
            ranked_obs(i) = ranked_obs(i) + 1;
        end
    end
    
end

end
