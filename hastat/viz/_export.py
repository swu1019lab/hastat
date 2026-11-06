"""Shared plot typography and editable vector exports."""
from contextlib import contextmanager
import matplotlib as mpl
from matplotlib.font_manager import FontProperties, findfont


@contextmanager
def plot_fonts():
    """Use Arial when available and restore the caller's settings afterwards."""
    try:
        findfont(FontProperties(family='Arial'), fallback_to_default=False)
        family = 'Arial'
    except ValueError:
        family = 'DejaVu Sans'
    settings = {
        'font.family': [family],
        'mathtext.fontset': 'custom',
        'mathtext.rm': family,
        'mathtext.it': family + ':italic',
        'mathtext.bf': family + ':bold',
        'mathtext.sf': family,
    }
    if 'mathtext.bfit' in mpl.rcParams:
        settings['mathtext.bfit'] = family + ':bold:italic'
    with mpl.rc_context(settings):
        yield


def save_figure(figure, *args, **kwargs):
    """Embed TrueType PDF fonts and retain SVG text for downstream editing."""
    with plot_fonts(), mpl.rc_context({'pdf.fonttype': 42, 'svg.fonttype': 'none'}):
        return figure.savefig(*args, **kwargs)
