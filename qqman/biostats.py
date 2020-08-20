import numpy as np
import numbers

def ppoints(n, a=None):
	""" numpy analogue or `R`'s `ppoints` function
		see details at https://www.rdocumentation.org/packages/stats/versions/3.6.2/topics/ppoints
		https://docs.tibco.com/pub/enterprise-runtime-for-R/5.0.0/doc/html/Language_Reference/stats/ppoints.html
		:param n: array type or number"""
	
	if isinstance(n, numbers.Number):
		n = np.float(n)
	else:
		n = np.float(len(n))
	if a == None:
		a = .375 if n<=10 else .5
		
	return (np.arange(n) + 1 - a)/(n + 1 - 2*a)