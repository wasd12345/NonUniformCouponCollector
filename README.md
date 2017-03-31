# NonUniformCouponCollector
Coupon Collector Problem for arbitrary probabilities and arbitrary quotas

Has 3 Functions:


NonUniformCouponCollectorQuotas__AnalyticalMean
-> For arbitrary coupon probabilities and arbitrary coupon quotas, calculates exact value of expected #draws to complete coupon collection
-> Slow for large #coupons and/or large quotas


UniformCouponCollectorQuotas__AnalyticalMean
-> For uniform coupon probabilities and uniform coupon quotas (all coupon quotas same but not necessarily 1), calculates exact value of expected #draws to complete coupon collection


NonUniformCouponCollectorQuotas__Simulation
-> For arbitrary coupon probabilities and arbitrary coupon quotas, do a simulation with arbitrarily many trials to get an empirical distribution of #draws to complete coupon collection. In particular, looks at summary statistics of mean, median, variance.
-> Analytical values for the mean, calculated using the 2 functions above, match the simulated values.
