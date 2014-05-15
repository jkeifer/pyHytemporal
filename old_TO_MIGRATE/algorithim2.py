__author__ = 'phoetrymaster'

import numpy
from scipy import optimize
from scipy import interpolate

bestguess = -50

#Measured values, pseudo-Argentina
valsf = {-175: -0.2, -159: -0.23, -143: -0.24, -127: -0.25, -111: -0.23, -95: -0.26, -79: -0.28, -63: -0.22,
         -47: -0.24, -31: -0.12, -15: 0.14, 1: 0.35, 17: 0.45, 33: 0.67, 49: 0.8, 65: 0.73, 81: 0.65, 97: 0.51,
         113: 0.27, 129: 0.01, 145: -0.1, 161: -0.19, 177: -0.21}

#Predicted to match, x3
#valsh = {-367: -0.2, -351: -0.2, -335: -0.2, -319: -0.2, -303: -0.2, -287: -0.2, -271: -0.2, -255: -0.2, -239: -0.1,
#         -223: 0.1, -207: 0.32, -191: 0.4, -175: 0.7, -159: 0.81, -143: 0.7, -127: 0.6, -111: 0.5, -95: 0.3, -79: 0,
#         -63: -0.1, -47: -0.2, -31: -0.2, -15: -0.2, 1: -0.2, 17: -0.2, 33: -0.2, 49: -0.2, 65: -0.2, 81: -0.2,
#         97: -0.2, 113: -0.2, 129: -0.1, 145: 0.1, 161: 0.32, 177: 0.4, 193: 0.7, 209: 0.81, 225: 0.7, 241: 0.6,
#         257: 0.5, 273: 0.3, 289: 0, 305: -0.1, 321: -0.2, 337: -0.2, 353: -0.2, 369: -0.2, 385: -0.2, 401: -0.2,
#         417: -0.2, 433: -0.2, 449: -0.2, 465: -0.2, 481: -0.2, 497: -0.1, 513: 0.1, 529: 0.32, 545: 0.4, 561: 0.7,
#         577: 0.81, 593: 0.7, 609: 0.6, 625: 0.5, 641: 0.3, 657: 0, 673: -0.1, 689: -0.2, 705: -0.2, 721: -0.2}

#Predicted to match, pseudo-Kansas, x1
valsh = {1: -0.2, 17: -0.2, 33: -0.2, 49: -0.2, 65: -0.2, 81: -0.2, 97: -0.2, 113: -0.2, 129: -0.1, 145: 0.1, 161: 0.32,
         177: 0.4, 193: 0.7, 209: 0.81, 225: 0.7, 241: 0.6, 257: 0.5, 273: 0.3, 289: 0, 305: -0.1, 321: -0.2, 337: -0.2,
         353: -0.2}


#Predicted to unmatch, psedo-Kansas
#valsh = {-367: -0.2, -351: -0.2, -335: -0.2, -319: -0.2, -303: 0, -287: 0.1, -271: 0.3, -255: 0.7, -239: 0.75,
#         -223: 0.79, -207: 0.72, -191: 0.68, -175: 0.64, -159: 0.57, -143: 0.35, -127: 0.29, -111: 0.1, -95: 0,
#         -79: -0.1, -63: -0.2, -47: -0.2, -31: -0.2, -15: -0.2, 1: -0.2, 17: -0.2, 33: -0.2, 49: -0.2, 65: 0, 81: 0.1,
#         97: 0.3, 113: 0.7, 129: 0.75, 145: 0.79, 161: 0.72, 177: 0.68, 193: 0.64, 209: 0.57, 225: 0.35, 241: 0.29,
#         257: 0.1, 273: 0, 289: -0.1, 305: -0.2, 321: -0.2, 337: -0.2, 353: -0.2, 369: -0.2, 385: -0.2, 401: -0.2,
#         417: -0.2, 433: 0, 449: 0.1, 465: 0.3, 481: 0.7, 497: 0.75, 513: 0.79, 529: 0.72, 545: 0.68, 561: 0.64,
#         577: 0.57, 593: 0.35, 609: 0.29, 625: 0.1, 641: 0, 657: -0.1, 673: -0.2, 689: -0.2, 705: -0.2, 721: -0.2}


##Meausred values, pseudo-Kansas
#valsf = {1: -0.25, 17: -0.23, 33: -0.26, 49: -0.28, 65: -0.22, 81: -0.24, 97: -0.12, 113: 0.14, 129: 0.35, 145: 0.45,
#         161: 0.67, 177: 0.8, 193: 0.73, 209: 0.65, 225: 0.51, 241: 0.27, 257: 0.01, 273: -0.1, 289: -0.19,
#         305: -0.21, 321: -0.2, 337: -0.23, 353: -0.24}
#
##Measured to unmatch, pseudo-Kansas
#valsf = {1: -0.2, 17: -0.2, 33: -0.2, 49: -0.2, 65: 0, 81: 0.1,
#         97: 0.3, 113: 0.7, 129: 0.75, 145: 0.79, 161: 0.72, 177: 0.68, 193: 0.64, 209: 0.57, 225: 0.35, 241: 0.29,
#         257: 0.1, 273: 0, 289: -0.1, 305: -0.2, 321: -0.2, 337: -0.2, 353: -0.2}

def get_sort_dates_values(vals, threshhold=-3000):
    """Gets the DOY dates (the keys) in a list from dictonary vals and sorts those, placing them in chronological order
    (list x0). Then the function iterates over these values and gets the corresponding values, thresholding values if
    they are lower than an optional threshhold value (-3000 default = NoData in MODIS imagery), then appending them to
    the list y. x and y are then returned."""

    x = vals.keys()
    x.sort()
    y = []

    for i in x:
        if vals[i] < threshhold:
            y.append(threshhold)
        else:
            y.append(vals[i])

    print x
    print y

    return x, y


def find_fit(valsf, valsh, threshhold):
    x0, y0 = get_sort_dates_values(valsf, threshhold=threshhold)
    x1, y1 = get_sort_dates_values(valsh)

    tck = interpolate.splrep(x1, y1)

    fun = lambda x: ((1 / 22.8125 * numpy.sum(
        (valsf[i] - (x[0] * interpolate.splev((x[1] * (i + x[2])), tck))) ** 2 for i in x0)) ** (
                         1. / 2))

    bnds = ((0.3, 1.5), (0.3, 1.5), (-180, 365))

    res = optimize.minimize(fun, (1, 1, bestguess), method='TNC', bounds=bnds)

    return res.fun, res.x

res, transforms = find_fit(valsf, valsh, threshhold=-.2)

print res
print 1.0 / transforms[0]
print 1.0 / transforms[1]
print transforms[2]
