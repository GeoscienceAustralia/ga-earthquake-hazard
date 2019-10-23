def weighted_avg_and_std(values, weights):
    """
    Return the weighted average and standard deviation.

    values, weights -- Numpy ndarrays with the same shape.
    """
    from numpy import average, sqrt

    wtaverage = average(values, weights=weights)
    variance = average((values-wtaverage)**2, weights=weights)  # Fast and numerically precise
    return (wtaverage, sqrt(variance))

# from excel at the moment..
# values from looking up source models and coresponding mags manually
values = [4.03E-04,
2.70E-04,
2.40E-04,
1.39E-04,
1.40E-04,
1.31E-04,
1.39E-04,
4.70E-05,
7.75E-05
]

weights = [0.204171946,
0.209746606,
0.282914027,
0.051538462,
0.051538462,
0.052751131,
0.091859729,
0.036076923,
0.019402715
]

mean, stdev = weighted_avg_and_std(values, weights)
print(mean, stdev)
print(1/mean, 1/stdev)

