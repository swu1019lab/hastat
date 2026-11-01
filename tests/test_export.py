"""Check exported text objects, font embedding and scoped settings."""
import io
import unittest
import xml.etree.ElementTree as ET

import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt

from hastat.viz._export import save_figure


class ExportValidation(unittest.TestCase):
    def setUp(self):
        self.fig, self.ax = plt.subplots()
        self.ax.text(.5, .5, 'Hap1 (n=25) π', fontsize=12)
        self.addCleanup(plt.close, self.fig)

    def test_pdf_embeds_truetype_instead_of_type3(self):
        stream = io.BytesIO()
        with mpl.rc_context({'pdf.fonttype': 3}):
            save_figure(self.fig, stream, format='pdf')
            self.assertEqual(mpl.rcParams['pdf.fonttype'], 3)
        content = stream.getvalue()
        self.assertIn(b'/CIDFontType2', content)
        self.assertIn(b'/FontFile2', content)
        self.assertIn(b'/ToUnicode', content)
        self.assertNotIn(b'/Type3', content)

    def test_svg_retains_literal_text(self):
        stream = io.BytesIO()
        with mpl.rc_context({'svg.fonttype': 'path'}):
            save_figure(self.fig, stream, format='svg')
            self.assertEqual(mpl.rcParams['svg.fonttype'], 'path')
        root = ET.fromstring(stream.getvalue())
        labels = [''.join(t.itertext()) for t in root.findall('.//{http://www.w3.org/2000/svg}text')]
        self.assertIn('Hap1 (n=25) π', labels)

    def test_settings_restored_after_save_error(self):
        with mpl.rc_context({'pdf.fonttype': 3, 'svg.fonttype': 'path'}):
            with self.assertRaises(ValueError):
                save_figure(self.fig, io.BytesIO(), format='invalid-format')
            self.assertEqual(mpl.rcParams['pdf.fonttype'], 3)
            self.assertEqual(mpl.rcParams['svg.fonttype'], 'path')

    def test_png_is_unchanged(self):
        direct, wrapped = io.BytesIO(), io.BytesIO()
        self.fig.savefig(direct, format='png')
        save_figure(self.fig, wrapped, format='png')
        self.assertEqual(direct.getvalue(), wrapped.getvalue())


if __name__ == '__main__':
    unittest.main()
