"""Reference surface LW fluxes for the radiative-transfer stages, per source.

The two sources store surface radiation very differently:

- ERA5: `strd`/`str` in the sfc file, accumulated J m-2 over the hour ENDING
  at valid_time (divide by 3600 for the 1-h mean W m-2); no explicit upward
  flux, so LWup = strd - str.
- MERRA-2: collection M2T1NXRAD (the `rad` file), 1-h time-averaged W m-2
  stamped at the CENTER of the window (HH:30). LWGAB is the surface-ABSORBED
  downwelling flux (= EMIS * LWdn) and LWGEM the emitted flux, so
  LWdn = LWGAB / EMIS and LWup = LWGEM + (1 - EMIS) * LWdn (emitted +
  reflected). For an instantaneous analysis at HH the two windows stamped
  (HH-1):30 and HH:30 bracket the instant; their mean is used (or whichever
  one exists at day edges).

Both references are 1-h means, never instantaneous — the stage-7c state-time
caveat applies to either source.
"""

from __future__ import annotations

import datetime as dt

import numpy as np

from .config import rad_path
from .io_era5 import open_era5


def load_surface_lw(cfg: dict, date: dt.date, when, sfc=None):
    """Reference surface LW fluxes around analysis instant `when`.

    Returns (lwdn, lwup, note): 2-D [W m-2] arrays on the analysis grid and a
    provenance string. For ERA5 pass the already-open sfc dataset selected at
    `when`; for MERRA-2 the day's rad file is opened here.
    """
    if cfg["source"] == "era5":
        lwdn = sfc["strd"].values / 3600.0
        lwup = lwdn - sfc["str"].values / 3600.0
        return lwdn, lwup, ("ERA5 strd/str 1-h accumulations ending at the "
                            "snapshot, /3600 to W m-2; LWup = strd - str")

    rad = open_era5(rad_path(cfg, date))
    when = np.datetime64(when)
    stamps = [when - np.timedelta64(30, "m"), when + np.timedelta64(30, "m")]
    have = [s for s in stamps if s in rad["valid_time"].values]
    if not have:
        raise KeyError(
            f"no M2T1NXRAD stamps at {stamps[0]} / {stamps[1]} in "
            f"{rad_path(cfg, date)} — re-run merra2_download.py for this day")
    sel = rad.sel(valid_time=have)
    emis = sel["emis"].values
    lwdn_t = sel["lwgab"].values / emis
    lwup_t = sel["lwgem"].values + (1.0 - emis) * lwdn_t
    note = ("MERRA-2 M2T1NXRAD 1-h means stamped "
            + ", ".join(str(s)[11:16] for s in have)
            + " (bracketing the instant); LWdn = LWGAB/EMIS, "
              "LWup = LWGEM + (1-EMIS)*LWdn")
    rad.close()
    return lwdn_t.mean(axis=0), lwup_t.mean(axis=0), note
