import numpy as np

def ebisuzaki(series):
    """
    Based on discussion in:
    Ebisuzaki, W., 1997:
    A Method to Estimate the Statistical Significance of a Correlation When the Data Are Serially Correlated.
    J. Climate, 10, 2147–2153. doi:10.1175/1520-0442(1997)010<2147:AMTETS>2.0.CO;2

    INPUT:
        series: 1D NumPy array (or list-like) of evenly spaced time or spatial series
        
    OUTPUT:
        new_series: a surrogate time series constructed by randomly altering the phase of the FFT
                    while preserving the original power spectrum, i.e. surrogate time series retain the same
                    autocorrelation structure and variance, but not the specific distribution of values.
    """

    # Convert input to NumPy array
    series = np.asarray(series)

    # Compute the discrete Fourier transform (DFT) of the input series
    fft_vals = np.fft.fft(series)

    # Extract magnitudes and phases (we will keep magnitudes, randomize phases)
    magnitudes = np.abs(fft_vals)
    phases = np.angle(fft_vals)

    N = len(series)
    Nf = (N - 1) // 2  # Number of positive frequency bins (excluding DC)

    # Generate Nf random phases between 0 and 2π
    random_phases = np.random.uniform(0, 2 * np.pi, Nf)
    phase_shifts = np.exp(1j * random_phases)  # Convert to complex unit vectors

    # Apply the random phase shift to the positive frequency components
    fft_vals[1:Nf+1] *= phase_shifts

    # Reconstruct negative frequencies using conjugate symmetry
    fft_vals[-Nf:] = np.conj(fft_vals[1:Nf+1][::-1])

    # Inverse FFT returns a complex array; take the real part
    new_series = np.fft.ifft(fft_vals).real

    return new_series
