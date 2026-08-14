"""Humidity conversions shared by the reanalysis sources and the radiosondes.

Kept source-agnostic on purpose: CARRA-2 needs RH -> q on ingest (it publishes
no specific humidity), the MOSAiC soundings need the same conversion for their
measured RH, and the profile comparison needs q -> RH so that every source is
expressed on one saturation convention rather than each product's own.

All functions use only true ufuncs, so numpy arrays and xarray DataArrays both
work and DataArray inputs keep their dims.
"""

from __future__ import annotations

import numpy as np

#: Ratio of the gas constants for dry air and water vapour.
EPSILON = 0.621981


def saturation_vapour_pressure(t_k, over: str = "water"):
    """Saturation vapour pressure [Pa] from temperature [K].

    Alduchov & Eskridge (1996) improved Magnus coefficients, accurate to the
    Arctic winter temperatures this pipeline works at (better than 0.4 % down
    to -80 C over ice).

    ``over``: 'water', 'ice', or 'mixed' — the IFS-style blend using ice below
    -23 C, water above 0 C and a quadratic ramp between. Radiosondes and CARRA
    both report RH over *water* even far below freezing, which is why 'water'
    is the default: mixing conventions would introduce a bias of tens of
    percent in RH at Arctic winter temperatures.
    """
    tc = t_k - 273.15
    e_water = 610.94 * np.exp(17.625 * tc / (tc + 243.04))
    e_ice = 611.21 * np.exp(22.587 * tc / (tc + 273.86))
    if over == "water":
        return e_water
    if over == "ice":
        return e_ice
    if over == "mixed":
        ramp = (t_k - 250.16) / (273.16 - 250.16)
        alpha = np.minimum(np.maximum(ramp, 0.0), 1.0) ** 2
        return alpha * e_water + (1.0 - alpha) * e_ice
    raise ValueError(f"unknown saturation reference {over!r} "
                     "(expected 'water', 'ice' or 'mixed')")


def specific_humidity_from_rh(rh_percent, t_k, p_pa, over: str = "water"):
    """Specific humidity [kg/kg] from relative humidity [%], T [K], p [Pa].

    Relative humidity is taken as the vapour-pressure ratio e/e_sat (the GRIB
    and ECMWF convention), so e = RH * e_sat and
    q = eps*e / (p - (1 - eps)*e).
    """
    e = np.maximum(rh_percent, 0.0) / 100.0 * saturation_vapour_pressure(t_k, over)
    e = np.minimum(e, 0.99 * p_pa)
    return EPSILON * e / (p_pa - (1.0 - EPSILON) * e)


def relative_humidity_from_q(q, t_k, p_pa, over: str = "water"):
    """Relative humidity [%] from specific humidity [kg/kg], T [K], p [Pa].

    The exact inverse of ``specific_humidity_from_rh``: inverting
    q = eps*e / (p - (1-eps)*e) for the vapour pressure gives
    e = q*p / (eps + (1-eps)*q).
    """
    q = np.maximum(q, 0.0)
    e = q * p_pa / (EPSILON + (1.0 - EPSILON) * q)
    return 100.0 * e / saturation_vapour_pressure(t_k, over)
