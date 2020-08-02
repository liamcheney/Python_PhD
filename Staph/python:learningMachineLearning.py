import argparse
from time import sleep as sl
import matplotlib
import matplotlib.pyplot as plt
import sklearn
import math


class LennardJonesian(object):
    """
    This is a class object that helps to define the functional form of the Lennard Jones
    potential, and allows you to calculate the value of the potential given a distances
    between the two particles.

    To construct a default object, simply call:
    lj = LennardJonesian()
    with the default parameters of a=1 and b=1 (see the __init__ function)

    If you would like to change these parameters, you can initialize a new object, such as
    lj2 = LennardJonesian(a=2,b=1.5)

    """

    def __init__(self, a=1, b=1):
        self.a = a
        self.b = b

    def potential(self, x):
        """
        This function calculate the value of the Lennard Jones potential at a given separation x.
        """
        return -2 * self.a / math.pow(x, 6) + self.b / math.pow(x, 12)

    def plot_analytic_potential(self, n_samples=500, x_start=0.64, x_end=1.64, c='k', label=None):
        """
        This function plots the curve for Lennard Jones potential.

        Args:
          xs: a list of x values.
        """
        xs = [x_start + i * (x_end - x_start) / n_samples for i in range(n_samples)]
        plt.plot(xs, [self.potential(i) for i in xs], '-', c=c, label=label, lw=2.5, alpha=0.5)
        plt.plot(xs, [0 for _ in xs], 'y:')
        plt.xlabel('$x$', fontsize=20)
        plt.ylabel('$y$', fontsize=20)

    def fake_samples(self, n_samples=500, x_start=0.64, x_end=1.64):
        """
        For the demonstration purpose only, we will generate a set of random samples (x,y)
        that will roughly follow the relationship given by the LJ potential. A random noise
        is added to the y value to make it deviate from the true LJ potential. The idea
        here is to mimic an experimentally/computationally generated database, on which
        ML will be performed. In reality, the database generating process would be
        the most laborious and slow part of the process, which cannot be demonstrated
        in a live event.

        Args:
          n_samples: number of samples to be generated, default is 500.
        """
        import random

        # this draws a random float in the range given by [x_start,x_end]
        x_sample = [x_start + random.randrange(0, 100000) / 100000 * (x_end - x_start) for _ in range(n_samples)]
        x_sample = list(sorted(x_sample))

        # this calculate the corresponding y at the x, plus a random noise
        y_sample = [self.potential(x) * (random.randrange(-50, 50) / 100 + 1) for x in x_sample]
        return x_sample, y_sample