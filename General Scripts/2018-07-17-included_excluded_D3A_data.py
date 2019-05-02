import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde
from time import sleep as sl

infile = open('/Users/hirokisuyama/Desktop/Book1.csv', 'r').read().splitlines()

length = []
cov = []
for line in infile:

    col = line.split(',')
    if col[0] != '':
        cov.append(col[0].strip('\ufeff'))
    if col[1] != '':
        length.append(col[1])

####Generate data
x = np.array(cov).astype(np.float)
y = np.array(length).astype(np.float)

#define the axis, data and colours
fig, ax = plt.subplots(figsize=(13,5))
ax.scatter(x[1745:], y[1745:], s=5, color='green') ##excluded data
ax.scatter(x[:1744], y[:1744], s=5, color='red') ##kept data


#customise the graph axis and add lengends etc
plt.xlabel('Genome Coverage')
plt.ylabel('Genome Size (million b.bp)')

plt.show()
