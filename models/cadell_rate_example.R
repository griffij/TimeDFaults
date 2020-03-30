##Example of calculating the present-day rate of the Cadell fault, 
##using Bayes Theorem


llh_active = 1-pexp(32000,1/7500)
llh_quiescent = 1-pexp(32000,1/1500000)

# Note we assume inform priors weights of 0.5 on active and quiescent states

p_active = llh_active*0.5/((llh_active + llh_quiescent)*0.5)
p_quiescent = llh_quiescent*0.5/((llh_active + llh_quiescent)*0.5)

sum_probs = p_active + p_quiescent
print(p_active)
print(p_quiescent)
print(sum_probs)
print('This should equal 1.0')

# Calculate probability in next 2500 years
p_2500 = p_active*pexp(2500,1/7500) + p_quiescent*pexp(2500,1/1500000)
print('Probability of earthquake in next 2500 years:')
print(p_2500)
print((1/p_2500))
print('Equivalent to annual rate of:')
annual_rate = p_2500/2500.
print(annual_rate)
print((1/annual_rate))