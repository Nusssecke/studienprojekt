using MathLink

# Import expansion coefficients
test = W"E^(1 / (6 * (1-u)))"
weval(test;u =0.5)