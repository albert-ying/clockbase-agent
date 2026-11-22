#!/usr/bin/env python3
"""Colour roles of Figure 3, one palette per entry of PALETTES, selected by the environment variable FIG3_PALETTE.

B is the Figure 3 palette and the one used when FIG3_PALETTE is unset (8 November 2026). 'july' holds the colours of
the July generators exactly; A and C are the other brighter candidates of the same review. Only colour changes
between palettes: geometry, type, marks, hatching and data are the same in all four. Fixed meanings in every palette:
lower clock output takes the cool hue and higher clock output the warm hue; the cool side of the volcano colour map
runs through the lower-bar hue and its warm side through the higher-bar hue; hatching and marker shape still give the
modality; July's alpha values are kept.

Roles
  lower, higher              bar fill for lower and higher clock outputs (b-g), key swatches
  lower_edge, higher_edge    bar outline and hatch strokes of those bars (hatching takes the edge colour)
  met_fill, met_edge         condition cell: condition met (b-e and key)
  met_text                   the reason text 'all three met' (bold, 6 pt)
  unknown_fill, unknown_edge condition cell: not stated or not agreed
  reason_text                reason text listing the unmet conditions
  ink                        text and the crossed cell of 'not met'
  diverging                  volcano point colour on Cohen's d and its colour bar: a matplotlib name, or colour
                             nodes [(position 0-1, hex)] for a linear segmented map
  density, density_alpha     density shading under the volcano cloud
  category                   volcano label text and category key, per intervention category
  names                      words used for the bar colours in the evidence README
Neutral roles (end-marker face and edge, non-significant points, point outlines, FDR line, key edges) are the same
in every palette and are listed in NEUTRAL.

Derived shades: a value commented as 'fill + x% black' (or white, or a named hex) is the fill mixed with that colour
by fraction x in sRGB, rounded to 8 bits. A category text colour taken from a named set is darkened only as far as
WCAG 4.5:1 on white requires.

Floors met by A, B and C (July is reported, not changed). `python fig3_palette_2026_11_08.py --check` prints every
measurement and lists the floors a palette misses.
  - lower vs higher fill: CIEDE2000 of at least 15 (nature-figure DELTA_E_FLOOR) and at least 8 under simulated
    deuteranopia and protanopia (CVD_DELTA_E_FLOOR; Machado et al. 2009, severity 1)
  - volcano map ends: CIEDE2000 of at least 8 under deuteranopia and under protanopia
  - bar outline and hatch strokes: at least 3:1 against their fill (WCAG 2.1 non-text contrast), so hatching keeps
    the modality readable at print size
  - category label text: at least 4.5:1 on white, and at least 15 CIEDE2000 from every colour of the volcano map
    (its nodes and 81 samples along it) and from the lower and higher fills, so a label hue is not read as a direction
  - condition met: 'all three met' text at least 4.5:1 on white, and the met fill at least 15 CIEDE2000 from the lower
    fill, so a met cell does not read as part of the bar beside it
The luminance contrast and the CIE L* gap of the lower and higher fills (what survives a greyscale print) are reported
beside them.
"""
import os
import sys

NEUTRAL = {
    "marker_edge": "#222", "marker_face": "#f7f7f7",          # bar-end modality marker (b-g)
    "key_marker_edge": "#111", "key_marker_face": "#eee",     # modality marker in the b-e key
    "nonsignificant": "lightgray", "point_outline": "black",   # volcano layers 1 and 2
    "fdr_line": "#666", "category_key_edge": "#333", "modality_key_face": "#999",
    "notmet_fill": "white",
}

PALETTES = {
    # July colours: Okabe-Ito blue and vermilion, ClockBase slate, RdBu_r on Cohen's d over rocket_r density.
    "july": {
        "lower": "#0072B2", "higher": "#D55E00",
        "lower_edge": "#111", "higher_edge": "#111",
        "met_fill": "#33414E", "met_edge": "#33414E", "met_text": "#0072B2",
        "unknown_fill": "#D9D9D9", "unknown_edge": "#B5B5B5", "reason_text": "#555555", "ink": "#222222",
        "diverging": "RdBu_r",
        "density": "rocket_r", "density_alpha": 0.35,
        "category": {"Disease": "#DC143C", "Drug": "#4B0082", "Genetic": "#228B22", "Environment": "#FF8C00",
                     "Other": "#808080"},
        "names": {"lower": "July blue", "higher": "July orange"},
    },
    # A: NPG (Nature Publishing Group) colours. Lower = NPG cyan with NPG navy strokes, higher = NPG red; condition met
    # = NPG navy; the map runs NPG navy, cyan, near-white, salmon, red, deep red. Category text keeps the July hue
    # analogues but stays off the map: Disease raspberry and Drug violet (NPG deep red and navy are the map ends),
    # Genetic NPG green (darkened 18% black for 4.5:1), Environment NPG brown.
    "A": {
        "lower": "#4DBBD5", "higher": "#E64B35",
        "lower_edge": "#3C5488", "higher_edge": "#501A13",          # higher edge: #E64B35 + 65% black
        "met_fill": "#3C5488", "met_edge": "#3C5488", "met_text": "#3C5488",
        "unknown_fill": "#D9D9D9", "unknown_edge": "#B5B5B5", "reason_text": "#555555", "ink": "#222222",
        "diverging": [(0.0, "#3C5488"), (1 / 3, "#4DBBD5"), (0.44, "#A6DDEA"),   # #4DBBD5 + 50% white
                      (0.5, "#F7F7F7"), (0.56, "#F39B7F"), (2 / 3, "#E64B35"), (1.0, "#DC0000")],
        "density": [(0.0, "#FFFFFF"), (1.0, "#5A5A5A")], "density_alpha": 0.35,
        "category": {"Disease": "#A01A58", "Drug": "#7B3294", "Genetic": "#00836F", "Environment": "#7E6148",
                     "Other": "#808080"},
        "names": {"lower": "NPG blue", "higher": "NPG red"},
    },
    # B, the Figure 3 palette: the July hue roles, brighter. Lower = azure, higher = NPG red; condition met = a deep
    # navy (NPG navy #3C5488 lies within 15 CIEDE2000 of the lower fill; see --check). The map runs deep azure,
    # the lower fill, pale blue, near-white, pale red, the higher fill and deep red, over a grey density as in A and C.
    # Disease raspberry and Environment ochre keep the category text off the map's reds and oranges.
    "B": {
        "lower": "#1B6FD1", "higher": "#E64B35",
        "lower_edge": "#092749", "higher_edge": "#501A13",          # fills + 65% black
        "met_fill": "#1F3A68", "met_edge": "#1F3A68", "met_text": "#1F3A68",
        "unknown_fill": "#D9D9D9", "unknown_edge": "#B5B5B5", "reason_text": "#555555", "ink": "#222222",
        "diverging": [(0.0, "#0E4C92"), (1 / 3, "#1B6FD1"), (0.44, "#A4C5ED"),   # #1B6FD1 + 60% white
                      (0.5, "#F7F7F7"), (0.56, "#F6B9A8"), (2 / 3, "#E64B35"), (1.0, "#9E2423")],
        "density": [(0.0, "#FFFFFF"), (1.0, "#5A5A5A")], "density_alpha": 0.35,
        "category": {"Disease": "#A01A58", "Drug": "#5B1A9E", "Genetic": "#0E6E3A", "Environment": "#B25F00",
                     "Other": "#808080"},
        "names": {"lower": "blue", "higher": "red"},
    },
    # C: teal and red. Lower = a light bright teal, higher = a darker vermilion red, so the two fills differ in
    # lightness as well as hue (teal and coral of equal lightness merge under protanopia and in greyscale). The map
    # runs petrol, teal, pale teal, near-white, pale red, red, dark red; the petrol cool end adds the blue that keeps
    # the two ends apart under protanopia. Condition met = deep teal. Teal and red take the direction hues, so the
    # category text uses raspberry, violet, navy (NPG) and NPG brown.
    "C": {
        "lower": "#1FB59A", "higher": "#D6402F",
        "lower_edge": "#0B3F36", "higher_edge": "#4B1610",          # fills + 65% black
        "met_fill": "#00695C", "met_edge": "#00695C", "met_text": "#00695C",
        "unknown_fill": "#D9D9D9", "unknown_edge": "#B5B5B5", "reason_text": "#555555", "ink": "#222222",
        # cool end: petrol #004D5C; warm end: higher + 60% #67001F (July's RdBu_r warm end); light nodes: fills + 60%
        # white
        "diverging": [(0.0, "#004D5C"), (1 / 3, "#1FB59A"), (0.44, "#A5E1D7"),
                      (0.5, "#F7F7F7"), (0.56, "#EFB3AC"), (2 / 3, "#D6402F"), (1.0, "#931A25")],
        "density": [(0.0, "#FFFFFF"), (1.0, "#5A5A5A")], "density_alpha": 0.35,
        "category": {"Disease": "#A01A58", "Drug": "#6A3D9A", "Genetic": "#3C5488", "Environment": "#7E6148",
                     "Other": "#808080"},
        "names": {"lower": "Teal", "higher": "red"},
    },
}

DEFAULT = "july"  # the Figure 3 palette (the original colours)


def name():
    """The palette chosen by FIG3_PALETTE (DEFAULT when unset)."""
    chosen = os.environ.get("FIG3_PALETTE", DEFAULT) or DEFAULT
    if chosen not in PALETTES:
        raise SystemExit(f"FIG3_PALETTE={chosen!r} is not one of {sorted(PALETTES)}")
    return chosen


def active():
    """Colour roles of the chosen palette, with the neutral roles added."""
    return {**NEUTRAL, **PALETTES[name()]}


def _cmap(spec, label):
    """A matplotlib name is passed through unchanged (July); colour nodes become a linear segmented map."""
    if isinstance(spec, str):
        return spec
    from matplotlib.colors import LinearSegmentedColormap
    return LinearSegmentedColormap.from_list(label, spec, N=256)


def diverging(pal=None):
    return _cmap((pal or active())["diverging"], f"fig3_diverging_{name()}")


def density(pal=None):
    return _cmap((pal or active())["density"], f"fig3_density_{name()}")


EDGE_CONTRAST_FLOOR = 3.0  # WCAG 2.1 non-text contrast


def check(names=None):
    """Measurements behind the floors in the module docstring, per palette. Colour distance is CIEDE2000; colour
    deficiency is simulated with the Machado et al. (2009) matrices at severity 1, as in nature-figure's
    simulate_cvd; contrast is the WCAG ratio. 'below_floor' lists every floor the palette misses."""
    sys.path.insert(0, "pipeline_assets/figure_style_2026_09_10")
    import matplotlib
    import numpy as np
    from matplotlib.colors import to_rgb
    from nature_figure_utils import (CVD_DELTA_E_FLOOR, DELTA_E_FLOOR, LABEL_CONTRAST_FLOOR, _srgb_to_lab, delta_e2000,
                                     simulate_cvd, wcag_contrast)
    kinds = ("deuteranopia", "protanopia")

    def lightness(c):
        return float(_srgb_to_lab(np.array(to_rgb(c)))[0])

    def sep(a, b):
        out = {"normal": delta_e2000(a, b)}
        out.update({k: delta_e2000(simulate_cvd(a, k), simulate_cvd(b, k)) for k in kinds})
        return out

    rows = []
    for key in names or PALETTES:
        p = {**NEUTRAL, **PALETTES[key]}
        cmap = matplotlib.colormaps[p["diverging"]] if isinstance(p["diverging"], str) else _cmap(p["diverging"], key)
        ramp = [matplotlib.colors.to_hex(cmap(v)) for v in np.linspace(0.0, 1.0, 81)]
        fills = sep(p["lower"], p["higher"])
        ends = {k: delta_e2000(simulate_cvd(ramp[0], k), simulate_cvd(ramp[-1], k)) for k in kinds}
        edges = {"lower": wcag_contrast(p["lower_edge"], p["lower"]), "higher": wcag_contrast(p["higher_edge"], p["higher"])}
        cats = {c: v for c, v in p["category"].items() if c != "Other"}
        text = {f"category {c}": wcag_contrast(v) for c, v in cats.items()}
        text["met_text"] = wcag_contrast(p["met_text"])
        off_map = {c: min(delta_e2000(v, m) for m in ramp) for c, v in cats.items()}
        nodes = [n for _, n in p["diverging"]] if not isinstance(p["diverging"], str) else ramp
        off_nodes = {c: min(delta_e2000(v, m) for m in nodes) for c, v in cats.items()}
        off_fills = {c: min(delta_e2000(v, p["lower"]), delta_e2000(v, p["higher"])) for c, v in cats.items()}
        met_lower = delta_e2000(p["met_fill"], p["lower"])
        lstar = {"lower": lightness(p["lower"]), "higher": lightness(p["higher"])}
        below = []
        if fills["normal"] < DELTA_E_FLOOR:
            below.append(f"lower vs higher fill dE {fills['normal']:.1f} < {DELTA_E_FLOOR:g}")
        below += [f"lower vs higher fill dE {k} {fills[k]:.1f} < {CVD_DELTA_E_FLOOR:g}" for k in kinds
                  if fills[k] < CVD_DELTA_E_FLOOR]
        below += [f"map ends dE {k} {v:.1f} < {CVD_DELTA_E_FLOOR:g}" for k, v in ends.items() if v < CVD_DELTA_E_FLOOR]
        below += [f"{k} edge vs fill {v:.2f}:1 < {EDGE_CONTRAST_FLOOR:g}:1" for k, v in edges.items()
                  if v < EDGE_CONTRAST_FLOOR]
        below += [f"{k} on white {v:.2f}:1 < {LABEL_CONTRAST_FLOOR:g}:1" for k, v in text.items()
                  if v < LABEL_CONTRAST_FLOOR]
        below += [f"category {c} dE {v:.1f} from the map < {DELTA_E_FLOOR:g}" for c, v in off_map.items()
                  if v < DELTA_E_FLOOR]
        below += [f"category {c} dE {v:.1f} from a map node < {DELTA_E_FLOOR:g}" for c, v in off_nodes.items()
                  if v < DELTA_E_FLOOR]
        below += [f"category {c} dE {v:.1f} from a bar fill < {DELTA_E_FLOOR:g}" for c, v in off_fills.items()
                  if v < DELTA_E_FLOOR]
        if met_lower < DELTA_E_FLOOR:
            below.append(f"met fill vs lower fill dE {met_lower:.1f} < {DELTA_E_FLOOR:g}")
        rows.append({"palette": key, "lower": p["lower"], "higher": p["higher"],
                     "dE_lower_higher": {k: round(v, 1) for k, v in fills.items()},
                     "luminance_contrast_lower_higher": round(wcag_contrast(p["lower"], p["higher"]), 2),
                     "lightness_Lstar_fills": {**{k: round(v, 1) for k, v in lstar.items()},
                                               "gap": round(abs(lstar["lower"] - lstar["higher"]), 1)},
                     "map_ends": ramp[0] + " " + ramp[-1],
                     "dE_map_ends": {k: round(v, 1) for k, v in ends.items()},
                     "edge_vs_fill_contrast": {k: round(v, 2) for k, v in edges.items()},
                     "met_fill": p["met_fill"], "dE_met_vs_lower": round(met_lower, 1),
                     "contrast_on_white": {k: round(v, 2) for k, v in text.items()},
                     "category_dE_from_map": {k: round(v, 1) for k, v in off_map.items()},
                     "category_dE_from_map_nodes": {k: round(v, 1) for k, v in off_nodes.items()},
                     "category_dE_from_fills": {k: round(v, 1) for k, v in off_fills.items()},
                     "below_floor": below})
    return rows


if __name__ == "__main__":
    if "--check" in sys.argv:
        import json
        print(json.dumps(check(), indent=1))
    else:
        print(name(), active())
