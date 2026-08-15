"""Horizontal-grid abstraction shared by the analysis stages.

Two grid families occur in the pipeline:

**regular** (ERA5, MERRA-2)
    dims ``latitude`` x ``longitude`` with 1-D coordinates of the same names.
    Relative cell area is proportional to cos(latitude).

**curvilinear** (CARRA-2)
    dims ``y`` x ``x`` on a north polar-stereographic projection, carrying
    2-D ``latitude(y, x)`` / ``longitude(y, x)`` coordinates. Relative cell
    area is proportional to (1 + sin(latitude))**2 -- see ``area_weights``.

Every stage that reduces, plots or samples a horizontal field goes through
these helpers instead of hard-coding ``("latitude", "longitude")``, so the
same code runs on either family.
"""

from __future__ import annotations

import numpy as np
import xarray as xr

EARTH_RADIUS_KM = 6371.0

#: CF grid-mapping variable written into normalized curvilinear files.
GRID_MAPPING_VAR = "polar_stereographic"


def hdims(obj: xr.Dataset | xr.DataArray) -> tuple[str, str]:
    """The object's two horizontal dimension names, outer (y-like) first.

    Raises ValueError if neither known layout is present.
    """
    dims = set(obj.dims)
    if {"latitude", "longitude"} <= dims:
        return ("latitude", "longitude")
    if {"y", "x"} <= dims:
        return ("y", "x")
    raise ValueError(
        f"no recognized horizontal dims in {sorted(dims)}; expected "
        "('latitude', 'longitude') or ('y', 'x')"
    )


def is_curvilinear(obj: xr.Dataset | xr.DataArray) -> bool:
    """True for projected grids with 2-D latitude/longitude coordinates."""
    return hdims(obj) == ("y", "x")


def hshape(obj: xr.Dataset | xr.DataArray) -> tuple[int, int]:
    """Size of the horizontal grid, in ``hdims`` order."""
    return tuple(obj.sizes[d] for d in hdims(obj))  # type: ignore[return-value]


def latlon2d(obj: xr.Dataset | xr.DataArray) -> tuple[np.ndarray, np.ndarray]:
    """2-D (latitude, longitude) arrays over the horizontal grid.

    Regular grids are broadcast from their 1-D coordinates; curvilinear grids
    already store 2-D coordinates and are returned in ``hdims`` order.
    """
    if is_curvilinear(obj):
        lat = obj["latitude"].transpose("y", "x").values
        lon = obj["longitude"].transpose("y", "x").values
        return np.asarray(lat, dtype=float), np.asarray(lon, dtype=float)
    lat = np.asarray(obj["latitude"].values, dtype=float)
    lon = np.asarray(obj["longitude"].values, dtype=float)
    return np.meshgrid(lat, lon, indexing="ij")


def domain_mask(obj: xr.Dataset | xr.DataArray) -> np.ndarray:
    """Boolean 2-D mask of the cells that belong to the requested area.

    A regular lat/lon request covers exactly the box asked for, so every cell
    counts. CARRA-2 is subset in projection space, so its delivered rectangle
    over-covers the latitude band and ``io_carra2.clip_to_area`` records which
    cells are genuinely inside as a ``domain_mask`` variable.
    """
    if "domain_mask" in getattr(obj, "variables", {}):
        return np.asarray(obj["domain_mask"].values, dtype=bool)
    return np.ones(hshape(obj), dtype=bool)


def box_mask(obj: xr.Dataset | xr.DataArray, area) -> np.ndarray:
    """Boolean 2-D mask of the cells inside a CDS-style ``[N, W, S, E]`` box.

    Works on either grid family, since it tests the (2-D) latitude/longitude
    of every cell. The longitude test is wrap-safe: E and W are compared as an
    eastward span from W, and a zero span (e.g. -180..180) means the full
    circle. Combine with ``domain_mask``/``area_weights`` for analysis-time
    subsetting to a domain smaller than the one downloaded.
    """
    lat_n, lon_w, lat_s, lon_e = (float(v) for v in area)
    lat2d, lon2d = latlon2d(obj)
    inside = (lat2d >= lat_s) & (lat2d <= lat_n)
    span = (lon_e - lon_w) % 360.0
    if span > 0.0:
        inside &= (lon2d - lon_w) % 360.0 <= span
    return inside


def area_tag(area) -> str:
    """Short filename tag for an analysis area, e.g. ``85-90N``."""
    lat_n, lon_w, lat_s, lon_e = (float(v) for v in area)
    tag = f"{lat_s:g}-{lat_n:g}N"
    if (lon_e - lon_w) % 360.0 != 0.0:
        tag += f"_{lon_w:g}-{lon_e:g}E"
    return tag


def area_weights(obj: xr.Dataset | xr.DataArray) -> np.ndarray:
    """Relative spherical cell area over the horizontal grid, as a 2-D array.

    Cells outside ``domain_mask`` get zero weight, so area statistics ignore
    the corners a projected bounding box adds. Only ratios matter (callers
    normalize by the weight sum), so the constants are dropped:

    * regular lat/lon -- cell area is proportional to cos(latitude).
    * north polar stereographic -- the point scale factor is
      k = (1 + sin(phi_c)) / (1 + sin(phi)) for standard parallel phi_c
      (Snyder 21-4, spherical), and true area is dx*dy / k**2, so the relative
      weight is (1 + sin(phi))**2. The standard parallel enters only as a
      constant factor and cancels, which is why this needs no projection
      metadata.
    """
    lat2d, _ = latlon2d(obj)
    if is_curvilinear(obj):
        w = (1.0 + np.sin(np.deg2rad(lat2d))) ** 2
    else:
        w = np.cos(np.deg2rad(lat2d))
    return np.where(domain_mask(obj), w, 0.0)


def weighted_mean(field: np.ndarray, weights: np.ndarray) -> float:
    """Area-weighted mean of a 2-D field, ignoring non-finite cells."""
    ok = np.isfinite(field)
    wsum = weights[ok].sum()
    if wsum == 0:
        return float("nan")
    return float((field[ok] * weights[ok]).sum() / wsum)


def horizontal_coords(obj: xr.Dataset | xr.DataArray) -> dict:
    """Coordinate mapping to attach to a dataset built on this grid.

    For regular grids that is the two 1-D coordinates; for curvilinear grids
    the projection coordinates plus the 2-D latitude/longitude fields.
    """
    dims = hdims(obj)
    coords = {d: obj[d] for d in dims if d in obj.coords}
    if is_curvilinear(obj):
        for name in ("latitude", "longitude"):
            coords[name] = (dims, latlon2d(obj)[0 if name == "latitude" else 1])
    return coords


def cell_spacing_km(obj: xr.Dataset | xr.DataArray) -> float:
    """Median great-circle distance between horizontally adjacent cells [km]."""
    lat2d, lon2d = latlon2d(obj)
    step = great_circle_km(lat2d[:, :-1], lon2d[:, :-1],
                           lat2d[:, 1:], lon2d[:, 1:])
    return float(np.nanmedian(step))


def same_grid(a: xr.Dataset, b: xr.Dataset, tol_fraction: float = 0.01) -> bool:
    """True when two datasets share a horizontal grid to within a fraction
    of a cell.

    Exact equality is too strict for CARRA-2: the CDS encodes latitude and
    longitude independently in each GRIB message, so the pressure-level and
    single-level deliveries of one grid disagree in the last few decimals —
    a couple of metres against 2500 m spacing, and up to 0.3 deg of longitude
    beside the pole where longitude is ill-conditioned. Comparing physical
    displacement rather than coordinate values sidesteps both.
    """
    if hdims(a) != hdims(b):
        return False
    lat_a, lon_a = latlon2d(a)
    lat_b, lon_b = latlon2d(b)
    if lat_a.shape != lat_b.shape:
        return False
    displacement = float(np.nanmax(great_circle_km(lat_a, lon_a, lat_b, lon_b)))
    return displacement <= tol_fraction * cell_spacing_km(a)


def grid_template(obj: xr.Dataset) -> xr.Dataset:
    """A small dataset carrying only this grid's description.

    Holds the horizontal coordinates plus, where present, the CF grid mapping
    and the domain mask — enough for ``hdims``, ``area_weights``,
    ``horizontal_coords`` and the plotting helpers, without keeping a month of
    fields alive while the daily files stream past. Selecting the grid
    variables alone is not enough: ERA5 and MERRA-2 have none, and the
    resulting dataset would carry no dims at all.
    """
    keep = [v for v in (GRID_MAPPING_VAR, "domain_mask") if v in obj.variables]
    template = obj[keep] if keep else xr.Dataset()
    return template.assign_coords(horizontal_coords(obj)).load()


def great_circle_km(lat1, lon1, lat2, lon2) -> np.ndarray:
    """Great-circle distance [km]; correct near the pole, unlike a flat
    lat/lon hypotenuse."""
    p1, p2 = np.deg2rad(lat1), np.deg2rad(lat2)
    dlam = np.deg2rad(np.asarray(lon2) - np.asarray(lon1))
    hav = (np.sin((p2 - p1) / 2) ** 2
           + np.cos(p1) * np.cos(p2) * np.sin(dlam / 2) ** 2)
    return 2 * EARTH_RADIUS_KM * np.arcsin(np.sqrt(np.clip(hav, 0, 1)))


def _unit_vectors(lat: np.ndarray, lon: np.ndarray) -> np.ndarray:
    """(N, 3) unit vectors on the sphere for KD-tree nearest-neighbour search.

    Chord distance in 3-D is monotonic in great-circle distance, so the
    3-D nearest neighbour is the true nearest neighbour -- which a lat/lon
    KD-tree would get wrong near the pole and across the date line.
    """
    phi, lam = np.deg2rad(lat.ravel()), np.deg2rad(lon.ravel())
    cos_phi = np.cos(phi)
    return np.column_stack([cos_phi * np.cos(lam), cos_phi * np.sin(lam),
                            np.sin(phi)])


class GridIndex:
    """Nearest-grid-cell lookup for either grid family.

    Build once per opened file and reuse across many queries -- the KD-tree
    over a 2.5 km CARRA-2 domain costs about a second to construct.

        idx = GridIndex(ds)
        (iy, ix), dist_km = idx.query(85.0, 120.0)
        column = ds.isel(**idx.isel(iy, ix))
    """

    def __init__(self, obj: xr.Dataset | xr.DataArray):
        self.dims = hdims(obj)
        self.curvilinear = is_curvilinear(obj)
        self.lat2d, self.lon2d = latlon2d(obj)
        self.shape = self.lat2d.shape
        self.inside = domain_mask(obj)
        self._tree = None
        self._flat: np.ndarray | None = None
        if not self.curvilinear:
            self._lat1d = np.asarray(obj["latitude"].values, dtype=float)
            self._lon1d = np.asarray(obj["longitude"].values, dtype=float)

    @property
    def tree(self):
        if self._tree is None:
            from scipy.spatial import cKDTree

            # only cells inside the requested area are candidates: the corners
            # of a projected bounding box hold no data, and a query near the
            # domain edge would otherwise snap to one and return silent NaNs
            self._flat = np.flatnonzero(self.inside.ravel())
            lat = self.lat2d.ravel()[self._flat]
            lon = self.lon2d.ravel()[self._flat]
            self._tree = cKDTree(_unit_vectors(lat, lon))
        return self._tree

    def query(self, lat: float, lon: float) -> tuple[tuple[int, int], float]:
        """Nearest in-domain cell to (lat, lon): ``(iy, ix)`` and distance [km].

        The distance is the caller's guard against a point outside the domain:
        the nearest cell always exists, so a far-away match shows up as a large
        distance rather than an error.
        """
        if self.curvilinear:
            vec = _unit_vectors(np.array([lat]), np.array([lon]))
            k = int(self.tree.query(vec)[1][0])
            iy, ix = np.unravel_index(self._flat[k], self.shape)
        else:
            iy = int(np.abs(self._lat1d - lat).argmin())
            dlon = (self._lon1d - lon + 180.0) % 360.0 - 180.0
            ix = int(np.abs(dlon).argmin())
        dist = float(great_circle_km(self.lat2d[iy, ix], self.lon2d[iy, ix],
                                     lat, lon))
        return (int(iy), int(ix)), dist

    def isel(self, iy: int, ix: int) -> dict[str, int]:
        """``dict`` of index selectors for ``Dataset.isel``."""
        return {self.dims[0]: iy, self.dims[1]: ix}

    def latlon(self, iy: int, ix: int) -> tuple[float, float]:
        """Cell-centre latitude/longitude of an index pair."""
        return float(self.lat2d[iy, ix]), float(self.lon2d[iy, ix])


def projection_crs(obj: xr.Dataset | xr.DataArray):
    """Cartopy CRS of a curvilinear grid's projection, or None.

    Reads the CF ``polar_stereographic`` grid-mapping variable written by
    ``io_carra2.normalize_carra2``. Returning None means the caller should
    fall back to plotting against 2-D latitude/longitude.
    """
    if not is_curvilinear(obj) or GRID_MAPPING_VAR not in getattr(obj, "variables", {}):
        return None
    import cartopy.crs as ccrs

    a = obj[GRID_MAPPING_VAR].attrs
    # The globe matters: x/y were computed on the file's sphere, and letting
    # cartopy default to WGS84 instead misplaces the data by up to ~19 km at
    # the edge of an 85-90N domain (zero at the pole). Declaring the sphere
    # makes the CRS self-consistent with the stored coordinates — a
    # round-trip through it reproduces the delivered latitude/longitude to
    # 0.2 m. Cartopy still applies a proper datum transformation when drawing
    # on WGS84-referenced axes, which is the CF-correct reading of the
    # grid mapping this file declares.
    radius = float(a.get("earth_radius", 6371229.0))
    globe = ccrs.Globe(semimajor_axis=radius, semiminor_axis=radius, ellipse=None)
    return ccrs.NorthPolarStereo(
        central_longitude=float(a.get("straight_vertical_longitude_from_pole", 0.0)),
        true_scale_latitude=float(a.get("standard_parallel", 90.0)),
        globe=globe,
    )
