# EXTENDED DATA FIG. 6 - merge the artwork onto the text layer.
#
# Clipping is lossless: the source page's content stream is wrapped in a clip path, so nothing is rasterised
# and the result stays fully editable vector art.
#
# Panels a and b are the top and bottom halves of the single block written by R/02_panels_ab.R (532 x 547 pt),
# clipped apart so each can be placed on its own. Panels c, d and e-g are clipped out of the bench artwork
# named by ED6_BENCH - a 612 x 792 pt page whose generating code is not in this repository.
#
# Usage: python python/compose.py
import os
from pypdf import PdfReader, PdfWriter
from pypdf.generic import DecodedStreamObject

REPO  = os.environ.get("ED6_REPO", os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
OUT   = os.environ.get("ED6_OUT", os.path.join(REPO, "figures"))
BASE  = os.path.join(OUT, "ED6 base.pdf")                            # from R/03_assemble_base.R
AB    = os.path.join(OUT, "ED6 panels ab - promoter occupancy.pdf")  # from R/02_panels_ab.R
BENCH = os.environ["ED6_BENCH"] if "ED6_BENCH" in os.environ else os.path.join(OUT, "bench_artwork.pdf")
FINAL = os.path.join(OUT, "Extended Data Fig. 6.pdf")

ITEMS = [
    # source, clip box in the source page, slot on the 612 x 960 pt output page
    (AB,    (0, 338, 532, 547),   (30, 710, 582, 922)),   # a  occupancy over the full ReMap universe
    (AB,    (0, 44, 532, 335),    (30, 410, 582, 700)),   # b  chromatin marks + state-TF binding sites
    (BENCH, (36, 272, 268, 470),  (30, 220, 295, 392)),   # c  sgESR1 knockdown efficiency
    (BENCH, (345, 535, 553, 709), (315, 220, 582, 392)),  # d  tamoxifen tumour growth
    (BENCH, (327, 282, 557, 468), (30, 74, 214, 208)),    # e  SomaScan, ASPC source, recovered scale
    (BENCH, (65, 57, 282, 239),   (214, 74, 398, 208)),   # f  SomaScan, NEPC source, recovered scale
    (BENCH, (333, 57, 570, 239),  (398, 74, 582, 208)),   # g  SomaScan, PCa source, CLR
]

def clipped(path, box):
    w = PdfWriter(clone_from=path); pg = w.pages[0]
    x0, y0, x1, y1 = box
    s = DecodedStreamObject()
    s.set_data(b"q %f %f %f %f re W n\n" % (x0, y0, x1 - x0, y1 - y0) + pg.get_contents().get_data() + b"\nQ")
    pg.replace_contents(s)
    return pg

def place(target, pg, box, slot):
    x0, y0, x1, y1 = box; a, b, c, d = slot
    sw, sh = c - a, d - b
    s = min(sw / (x1 - x0), sh / (y1 - y0))
    target.merge_transformed_page(pg, (s, 0, 0, s,
                                       a + (sw - s * (x1 - x0)) / 2 - s * x0,
                                       b + (sh - s * (y1 - y0)) / 2 - s * y0))

w = PdfWriter(); w.append(PdfReader(BASE)); t = w.pages[0]
for src, box, slot in ITEMS:
    place(t, clipped(src, box), box, slot)
with open(FINAL, "wb") as fh: w.write(fh)
print("written", FINAL)
