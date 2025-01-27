# Zernike Polynomial Calculator

## Overview
**Zernike** is a Python package for the evaluation of Zernike polynomials. It includes utilities for computing and plotting Zernike polynomials, as well as calculating Zernike transform coefficients for data on the unit disk. 

This package is particularly useful in fields such as optics, image processing, and wavefront analysis.

## Features
- Compute radial and full Zernike polynomials.
- Plot Zernike polynomials up to a specified degree.
- Perform Zernike transforms to calculate coefficients of a function on the unit disk.
- Load and process CSV data containing \(X, Y, Z\) values for Zernike analysis.

---

## Getting Started

### Prerequisites
Make sure the following Python libraries are installed:
- `numpy`
- `pandas`
- `scipy`
- `matplotlib`

Install them using pip if not already installed:
```bash
pip install -r requirements.txt
```


---

## Usage

### Importing the Module
```python
from zernike import ZernikeCalculator
```

### Example Usage
1. **Initialize the Zernike Calculator**:
   ```python
   zc = ZernikeCalculator()
   ```

2. **Compute and Plot Zernike Polynomials**:
   ```python
   zc.plot_polynomials(max_zern=21)
   ```

3. **Load a CSV File**:
   Ensure your CSV file contains columns `X`, `Y`, and `Z`. Load it as follows:
   ```python
   zc.load_csv("data.csv", resolution=1024)
   ```

4. **Perform Zernike Transform**:
   Compute Zernike coefficients up to a maximum radial degree:
   ```python
   coefficients = zc.zernike_transform(n_max=6)
   print(coefficients)
   ```

---

## Input Format
The input CSV file should have the following structure:

| X     | Y     | Z      |
|-------|-------|--------|
| -0.5  | 0.2   | 0.003  |
| 0.1   | -0.3  | -0.002 |
| ...   | ...   | ...    |

- `X` and `Y` must represent normalized coordinates in the range \([-1, 1]\).
- `Z` is the corresponding value at each \((X, Y)\) point.

---

## Methods

### `radial_polynomial(n, m, r)`
Computes the radial part of the Zernike polynomial for given \(n\), \(m\), and \(r\).

### `zernike_polynomial(n, m, r, theta)`
Computes the full Zernike polynomial for given \(n\), \(m\), \(r\), and \(\theta\).

### `plot_polynomials(max_zern=21)`
Plots Zernike polynomials up to the specified maximum index.

### `load_csv(csv_file_path, resolution=1024)`
Loads data from a CSV file and processes it for Zernike analysis.

### `zernike_transform(n_max=6)`
Calculates Zernike coefficients for the current loaded data.

---

## License
This project is licensed under the [MIT License](LICENSE).

---

## Author
- **Itay Davidovitch** - Zeiss SMT  
- **Email**: [dv.itay@gmail.com](mailto:dv.itay@gmail.com)  

Feel free to reach out for questions or suggestions!
