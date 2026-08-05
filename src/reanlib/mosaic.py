"""MOSAiC level-2 radiosonde profiles (Maturilli et al. 2021, CC-BY-4.0).

Monthly tab files from the PANGAEA series doi:10.1594/PANGAEA.928656 are
auto-downloaded to data/mosaic/soundings/ on first use. Callers must handle
a None return (month not in SOUNDING_SETS, download failure, or no sounding
in the requested window).
"""

from __future__ import annotations

from .config import REPO_ROOT

SOUNDING_DIR = REPO_ROOT / "data" / "mosaic" / "soundings"
SOUNDING_SETS = {
    (2020, 1): ("PS122_2_radiosonde_202001.tab",
                "https://doi.pangaea.de/10.1594/PANGAEA.928659?format=textfile"),
}
_CACHE: dict = {}


def load_sounding_month(year: int, month: int):
    """Cached DataFrame of every level-2 sounding row in the month, or None."""
    import pandas as pd

    key = (year, month)
    if key in _CACHE:
        return _CACHE[key]
    if key not in SOUNDING_SETS:
        print(f"note: no MOSAiC level-2 sounding file registered for "
              f"{year}-{month:02d} (reanlib/mosaic.py SOUNDING_SETS)")
        _CACHE[key] = None
        return None
    name, url = SOUNDING_SETS[key]
    path = SOUNDING_DIR / name
    if not path.exists():
        SOUNDING_DIR.mkdir(parents=True, exist_ok=True)
        print(f"downloading MOSAiC radiosonde level-2 data (~60 MB) -> {name} ...")
        try:
            import urllib.request
            urllib.request.urlretrieve(url, path)
        except Exception as e:
            print(f"  download failed ({e})")
            _CACHE[key] = None
            return None
    with open(path) as f:
        skip = next(i for i, line in enumerate(f) if line.startswith("*/")) + 1
    df = pd.read_csv(path, sep="\t", skiprows=skip,
                     usecols=["Date/Time", "Altitude [m]", "PPPP [hPa]",
                              "TTT [°C]"])
    df["Date/Time"] = pd.to_datetime(df["Date/Time"], format="ISO8601")
    _CACHE[key] = df
    return df


def sounding_profile(when, window_min: float = 90.0, z_max_m: "float | None" = None):
    """Rows of the sounding nearest `when`, altitude-sorted, or None.

    Launches are ~6 h apart, so a +-90 min window isolates one ascent.
    """
    import pandas as pd

    ts = pd.Timestamp(when)
    df = load_sounding_month(ts.year, ts.month)
    if df is None:
        return None
    m = (df["Date/Time"] - ts).abs() <= pd.Timedelta(minutes=window_min)
    if z_max_m is not None:
        m &= df["Altitude [m]"] <= z_max_m
    sel = df[m].sort_values("Altitude [m]")
    return sel if len(sel) else None
