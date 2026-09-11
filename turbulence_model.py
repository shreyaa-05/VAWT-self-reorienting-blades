"""
turbulence_model.py
-------------------
A lightweight, reproducible synthetic turbulence field for engineering
experiments around the VAWT model.

Project connection:
    "Exploiting Turbulence for Furthering Engineering"

The field is deliberately a SURROGATE, not a CFD turbulence model. It creates
zero-mean, spatially correlated longitudinal wind fluctuations as a sum of
Fourier modes. That makes it useful for controlled numerical experiments:
we can vary turbulence intensity and correlation structure, then measure how
an energy-conversion system responds.
"""

from dataclasses import dataclass
import numpy as np


@dataclass
class SyntheticTurbulence:
    """Reproducible correlated longitudinal turbulence for a rotor revolution."""
    intensity: float = 0.10       # RMS fluctuation / mean wind speed
    n_modes: int = 8              # number of azimuthal Fourier modes
    spectral_exponent: float = 1.0 # higher -> more energy in large scales
    seed: int = 7
    strength: float = 1.0         # controlled turbulence shaping multiplier

    def __post_init__(self):
        if not 0.0 <= self.intensity <= 0.5:
            raise ValueError("intensity must be between 0 and 0.5")
        if self.n_modes < 1:
            raise ValueError("n_modes must be >= 1")
        if self.strength < 0.0:
            raise ValueError("strength must be >= 0")

        rng = np.random.default_rng(self.seed)
        self._phase = rng.uniform(0.0, 2.0 * np.pi, self.n_modes)
        k = np.arange(1, self.n_modes + 1, dtype=float)
        self._weights = 1.0 / np.power(k, 0.5 * self.spectral_exponent)

        # Normalize the unscaled field to unit RMS, so `intensity` has a clear
        # interpretation independent of the number of Fourier modes.
        sample_theta = np.linspace(0.0, 2.0 * np.pi, 8192, endpoint=False)
        raw = self._raw(sample_theta)
        self._rms = float(np.sqrt(np.mean(raw ** 2)))

    def _raw(self, theta):
        theta = np.asarray(theta, dtype=float)
        k = np.arange(1, self.n_modes + 1, dtype=float)
        return np.sum(
            self._weights[:, None] *
            np.sin(k[:, None] * theta.reshape(1, -1) + self._phase[:, None]),
            axis=0,
        ) if theta.ndim else float(np.sum(self._weights * np.sin(k * theta + self._phase)))

    def normalized_fluctuation(self, theta, blade_phase=0.0):
        """Return a zero-mean, approximately unit-RMS correlated fluctuation."""
        raw = self._raw(np.asarray(theta) + blade_phase)
        return raw / self._rms

    def local_speed(self, mean_speed, theta, blade_phase=0.0):
        """Apply the synthetic longitudinal fluctuation to a local mean speed."""
        fluct = self.normalized_fluctuation(theta, blade_phase=blade_phase)
        speed = mean_speed * (1.0 + self.strength * self.intensity * fluct)
        return np.maximum(speed, 0.05 * mean_speed)

    def statistics(self, n=720):
        """Return basic statistics of the generated normalized turbulence field."""
        theta = np.linspace(0.0, 2.0 * np.pi, n, endpoint=False)
        u = self.normalized_fluctuation(theta)
        return {
            "mean": float(np.mean(u)),
            "rms": float(np.sqrt(np.mean(u ** 2))),
            "min": float(np.min(u)),
            "max": float(np.max(u)),
        }


if __name__ == "__main__":
    t = SyntheticTurbulence(intensity=0.10)
    print("Synthetic turbulence statistics:", t.statistics())
