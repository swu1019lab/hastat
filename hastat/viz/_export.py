"""Preserve text in vector exports without changing the caller's settings."""
import matplotlib as mpl


def save_figure(figure, *args, **kwargs):
    """Embed TrueType PDF fonts and retain SVG text for downstream editing."""
    with mpl.rc_context({'pdf.fonttype': 42, 'svg.fonttype': 'none'}):
        return figure.savefig(*args, **kwargs)
