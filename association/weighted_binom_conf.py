import numpy as np
import scipy.stats

def weighted_binom_conf(weights, successes, confidence):
    r'''
    return a weighted binomial confidence interval
    using the Wilson method described here
    https://en.wikipedia.org/wiki/Binomial_proportion_confidence_interval#Wilson_score_interval
    except rederiving for a weighted sum \sum_i w_i X_i
    where the X_i are iid binomial variables

    Derivation starts with
    \sigma = \sqrt(p(1-p)}
    and
    \frac{ \sum_i w_i (X_i - p) } { \sigma \sqrt{ \sum_i w_i^2 } } \sim z
    where z \sim N(0, 1)
    substituting in \sigam, solving for p using the quadratic formula
    and getting
    p = \frac 1 {1 + c^2/t^2}(\hat p + \frac c^2 {2 t^2}) +/-
      \frac{c/t}{1 + c^2/t^2} \sqrt{\hat p (1 - \hat p) + c^2/(4t^2)}
    where t = \sum_i w_i (e.g. n for noneighted problem)
    \hat p = \sum_i w_i X_i / t
    c = z \sqrt{ \sum_i w_i^2 }
    z = F_{N(0,1)}^{-1}(1 - \alpha)
    where
    \alpha is the desired confidence interval

    returns (mean prob, lower ci bound, upper ci bound)

    spot testing on 2021/08/04 shows this code is the same as
    statsmodels.stats.proportion.proportion_confint(method='wilson')
    when all the weights are 1,
    that setting the weigths of a bunch of elements to near zero
    gives about the same effect as them not being there (but slightly
    smaller intervals)
    and setting the weights of a bunch of elements to 0.5 gives about
    the same effect as having half as many of those observations, but
    slightly smaller intervals (in the case where that class of observations
    is already overrepresented)
    '''
    assert weights.shape == successes.shape
    assert len(weights.shape) == 1
    t = np.sum(weights)
    phat = np.dot(weights, successes)/t
    z = scipy.stats.norm.ppf(1 - confidence/2)
    c = z * np.sqrt(np.dot(weights, weights))

    divisor = 2 + 2*c**2/t**2
    center = (2*phat + c**2/t**2)/divisor
    interval_size = c/t*np.sqrt(4*phat*(1-phat) + c**2/t**2)/divisor

    return (phat, center - interval_size, center + interval_size)


