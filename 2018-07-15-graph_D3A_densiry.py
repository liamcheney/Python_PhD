import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde

infile = open('/Users/liamcheneyy/Desktop/Workbook1.csv','r').read().splitlines()

length = []
para = []
for line in infile:
    col = line.split(',')
    length.append(col[1].lstrip('\ufeff').lstrip('"'))
    para.append(col[0].lstrip('\ufeff').lstrip('"'))

# Generate data
x = np.array(length).astype(np.float)
y = np.array(para).astype(np.float)

# Calculate the point density
xy = np.vstack([x,y])
z = gaussian_kde(xy)(xy)


idx = z.argsort()
x, y, z = x[idx], y[idx], z[idx]

#define the axis, data and colours
fig, ax = plt.subplots(figsize=(13,5))
sc = ax.scatter(x, y, c=z, s=100, edgecolor='', vmin=0, vmax=1)
# sc1 = ax.scatter(x[6:50], y[6:50], s=100, color='red')
# sc2 = ax.scatter(x[51:100], y[51:100], s=100, color='green')

#customise the graph axis and add lengends etc
plt.xlabel('Number of Genuine Paralogs and Falsely Split Core Genes Each Genome')
plt.ylabel('Genome Size (million b.bp)')
cbar = plt.colorbar(sc)
cbar.set_ticks([0, 0.25, 0.5, 0.75, 1])


plt.show()
