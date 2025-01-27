### --------------------------------------- ###
#-#                 ZERNIKE                 #-#
### --------------------------------------- ###

"""
Python package for the evaluation of Zernike polynomials

Date: Jan 2025
Author: Itay Davidovitch - Zeiss SMT
Email: dv.itay@gmail.com
Version: 1.0
Description: this package implements several methods to compute and plot
Zernike polynomials. It also provides a method to calculate the Zernike transform
coefficients of some function on the unit disk.
"""


import numpy as np
import pandas as pd
from scipy.special import factorial
from scipy.interpolate import griddata
import matplotlib.pyplot as plt


class ZernikeCalculator:
    
    def __init__(self):
        self.loaded_csv = None # Last loaded raw csv file
        self.map = None        # Last processed total map
        self.resolution = 1024 # Resolution of grid on unit circle (1024 by default)


    def radial_polynomial(self, n, m, r):
        R = np.zeros_like(r)
        for k in range((n - abs(m)) // 2 + 1):
            R += r**(n - 2*k) * (-1)**k * factorial(n - k) / \
                 (factorial(k) * factorial((n + abs(m)) // 2 - k) * factorial((n - abs(m)) // 2 - k))
        return R


    def zernike_polynomial(self, n, m, r, theta):
        """
        Calculate the full Zernike polynomial.

        Args:
            n (int): Radial degree.
            m (int): Azimuthal frequency.
            rho (numpy.ndarray): Radial coordinate.
            theta (numpy.ndarray): Angular coordinate.

        Returns:
            numpy.ndarray: Zernike polynomial values.
        """
        radial = self.radial_polynomial(n, m, r)
        if m >= 0:
            return radial * np.cos(m * theta)
        else:
            return radial * np.sin(-m * theta)


    def plot_polynomials(self, max_zern=21):
        """
        Plot the Zernike polynomials up to max_zern.

        Itay: I might later make this a pyramid shape for visual clarity.
        """

        # Get grid
        r = np.linspace(0, 1, 100)
        theta = np.linspace(0, 2 * np.pi, 100)
        R, Theta = np.meshgrid(r, theta)
        X, Y = R * np.cos(Theta), R * np.sin(Theta)

        # Set up the plot
        number_of_rows = (max_zern + 5) // 6  # Determine the number of rows for 6 plots per row
        fig, axes = plt.subplots(number_of_rows, 6, figsize=(20, 4 * number_of_rows), subplot_kw={'projection': 'polar'})
        axes = axes.flatten()  # Flatten the axes array for easier indexing

        zern_counter = 0
        n = 0

        while zern_counter < max_zern:
            for m in range(-n, n + 1, 2):
                if zern_counter >= max_zern:
                    break
                Z = self.zernike_polynomial(n, m, R, Theta)
                ax = axes[zern_counter]
                ax.contourf(Theta, R, Z, 100, cmap='RdBu')
                ax.set_yticklabels([])
                ax.set_xticklabels([])
                ax.set_title(rf"$Z_{{{n}}}^{{{m}}}$", fontsize=18)

                zern_counter += 1
            n += 1
        
        # Hide unused subplots
        for ax in axes[zern_counter:]:
            ax.axis('off')


# Zernike Transform
    def load_csv(self, csv_file_path, resolution=1024):
        
        try:
            df = pd.read_csv(csv_file_path)
        except Exception as e:
            raise ValueError(f"Error reading CSV file: {e}")

        required_columns = ['X', 'Y', 'Z']
        if not all(col in df.columns for col in required_columns):
            raise ValueError(f"CSV file must contain the following columns: {required_columns}")

        X = df['X'].values.astype(np.float32)
        Y = df['Y'].values.astype(np.float32)
        Z = df['Z'].values.astype(np.float32)

        # Data processing
        Z = (Z - np.mean(Z)) * 1000
        # Scale X,Y to fit unit disk
        X = 2 * (X - np.min(X)) / (np.max(X) - np.min(X)) - 1
        Y = 2 * (Y - np.min(Y)) / (np.max(Y) - np.min(Y)) - 1
        # Define X,Y meshgrid with given resolution
        X_grid, Y_grid = np.mgrid[min(X):max(X):resolution*1j, min(Y):max(Y):resolution*1j]
        Z_grid = griddata((X, Y), Z, (X_grid, Y_grid), method='linear')

        # Save map, filepath and grid in module
        self.map = Z_grid  # Store the last processed total map
        self.loaded_csv = csv_file_path  # Store the path of the loaded CSV
        self.resolution = resolution  # Store the resolution of the grid

        return

    # def interpolate_to_polar(self, r, theta):
    #     """
    #     Interpolate the map to polar coordinates (r, theta).

    #     Args:
    #         r (numpy.ndarray): Radial coordinate.
    #         theta (numpy.ndarray): Angular coordinate.

    #     Returns:
    #         numpy.ndarray: Interpolated map.
    #     """
    #     values = np.linspace(-1, 1, self.resolution)
    #     return griddata((values, values), self.map, (r * np.cos(theta), r * np.sin(theta)), method='linear', fill_value=0)


    # def integrand(self, n, m, r, theta):
    #     return self.interpolate_to_polar(r, theta) * self.zernike_polynomial(n, m, r, theta) * r



    def zernike_transform(self, n_max=6):
        """
        Calculate the Zernike coefficients for a given function corresponding to wafer heights.

        Args:
            func (callable): Function to project onto the Zernike polynomials, defined on (r, theta).
            n_max (int): Maximum radial degree of the Zernike polynomials.

        Returns:
            dict: Dictionary of Zernike coefficients {(n, m): coefficient}.
        """

        x, y = np.linspace(-1, 1, self.resolution), np.linspace(-1, 1, self.resolution)
        dx, dy = x[1] - x[0], y[1] - y[0]

        # Create meshgrid
        X, Y = np.meshgrid(x, y)
        R, Theta = np.sqrt(X**2 + Y**2), np.arctan2(Y, X)

        # Filter points outside the unit disk and exclude NaN values
        mask = (R <= 1) & (~np.isnan(self.map))

        coefficients = {}

        for n in range(n_max+1):
            for m in range(-n, n + 1, 2):
                # Compute the Zernike polynomial on the cartesian grid
                zernike = self.zernike_polynomial(n, m, R, Theta)
                # Compute the integral
                coeff = np.sum(zernike[mask] * self.map[mask]) * dx * dy
                coefficients[(n, m)] = coeff * (2*n + 2) / np.pi  # Scale by the normalization factor
                if m==0:
                    coefficients[(n, m)] /= 2 # Neumann factor 

        return coefficients
    

if __name__ == "__main__":

    pass