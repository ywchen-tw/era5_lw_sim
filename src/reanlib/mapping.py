"""Shared polar-stereographic map panel used by the plotting scripts."""

from __future__ import annotations

import matplotlib.path as mpath
import matplotlib.pyplot as plt
import numpy as np

C_MASKED = "#d9d9d9"


def add_cyclic(field: np.ndarray, lon: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Close the -180/180 seam by repeating the first longitude column."""
    return (np.concatenate([field, field[:, :1]], axis=1),
            np.append(lon, lon[0] + 360.0))


def polar_panel(ax, lat: np.ndarray, lon: np.ndarray, field: np.ndarray, *,
                kind: str = "seq", cbar_label: str = "", south_lat: float | None = None,
                vmin: float | None = None, vmax: float | None = None,
                cmap=None, coast_color: str | None = None,
                cbar_label_size: float = 12, cbar_tick_size: float = 10):
    """North-polar-stereographic pcolormesh panel with circular boundary.

    kind 'seq': sequential colormap (viridis), vmin 0, vmax = 99th percentile
    unless given. kind 'div': diverging (RdBu_r), symmetric about 0. Masked
    (NaN) cells show the gray axes background. Returns the mesh.
    """
    import cartopy.crs as ccrs

    if south_lat is None:
        south_lat = float(np.min(lat))
    field_c, lon_c = add_cyclic(np.asarray(field), np.asarray(lon))

    if kind == "seq":
        cmap = cmap or plt.get_cmap("viridis")
        if vmax is None:
            vmax = max(float(np.nanpercentile(field_c, 99)), 1.0)
        vmin = 0.0 if vmin is None else vmin
        extend = "max"
    else:
        cmap = cmap or plt.get_cmap("RdBu_r")
        if vmax is None:
            vmax = max(float(np.nanpercentile(np.abs(field_c), 99)), 1.0)
        vmin = -vmax if vmin is None else vmin
        extend = "both"
    # masked cells: transparent NaNs over a gray axes background (cartopy
    # requires the colormap's 'bad' to stay transparent for wrapped meshes)
    ax.set_facecolor(C_MASKED)

    mesh = ax.pcolormesh(lon_c, lat, np.ma.masked_invalid(field_c),
                         transform=ccrs.PlateCarree(), cmap=cmap,
                         vmin=vmin, vmax=vmax, shading="nearest")
    # white coastlines on the dark sequential colormap, dark on the light one
    if coast_color is None:
        coast_color = "#ffffff" if kind == "seq" else "#444444"
    ax.coastlines(resolution="50m", lw=0.7, color=coast_color)
    ax.gridlines(draw_labels=False, lw=0.4, color="#bbbbbb",
                 xlocs=range(-180, 181, 45), ylocs=[80, 85])

    ax.set_extent([-180, 180, south_lat, 90], crs=ccrs.PlateCarree())
    theta = np.linspace(0, 2 * np.pi, 200)
    circle = mpath.Path(np.column_stack([np.sin(theta), np.cos(theta)]) * 0.5 + 0.5)
    ax.set_boundary(circle, transform=ax.transAxes)

    # manual meridian labels: cartopy's auto labels collide at the circular edge
    style = dict(transform=ccrs.PlateCarree(), fontsize=8, color="#666666",
                 ha="center", va="center", clip_on=False)
    label_lat = south_lat - 0.6
    ax.text(0, label_lat, "0°", **style)
    ax.text(180, label_lat, "180°", **style)
    ax.text(90, label_lat, "90°E", rotation=-90, **style)
    ax.text(-90, label_lat, "90°W", rotation=90, **style)
    ax.text(-135, 85, "85°N", transform=ccrs.PlateCarree(), fontsize=8,
            color="#666666", ha="center", va="center", alpha=0.9)

    cb = plt.colorbar(mesh, ax=ax, orientation="horizontal", pad=0.05,
                      fraction=0.045, extend=extend)
    cb.set_label(cbar_label, fontsize=cbar_label_size)
    cb.ax.tick_params(labelsize=cbar_tick_size)
    return mesh
