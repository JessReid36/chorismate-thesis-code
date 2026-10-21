"""build_charges.py - the optimal-charges page from template.html (beside this script) and <work_dir>/data.json.
Usage: python3 build_charges.py <work_dir> <plotly.min.js>
Writes <work_dir>/optimal_charges.html (Plotly from the jsDelivr CDN) and <work_dir>/chorismate_optimal_charges.html
(Plotly inlined, so it opens offline). Same construction and the same pinned library as site_effect_map/build.py."""
import hashlib, os, sys
WORK, LIB = sys.argv[1], sys.argv[2]
HERE = os.path.dirname(os.path.abspath(__file__))
PLOTLY_SHA = "6d21266ce1bd7d9e5ab4e115989c70c20de0382fd973a8f26ab58619eba4d603"   # plotly.js-dist-min 2.35.2
t = open(os.path.join(HERE, "template.html")).read()
data = open(os.path.join(WORK, "data.json")).read()
assert t.count("__DATA__") == 1
page = t.replace("__DATA__", data)
open(os.path.join(WORK, "optimal_charges.html"), "w").write(page)
cdn = '<script src="https://cdn.jsdelivr.net/npm/plotly.js-dist-min@2.35.2/plotly.min.js"></script>'
assert page.count(cdn) == 1
assert hashlib.sha256(open(LIB, "rb").read()).hexdigest() == PLOTLY_SHA, LIB + " is not plotly.js-dist-min 2.35.2's plotly.min.js"
lib = open(LIB).read()
assert "</script" not in lib
head, body = page.split('<div class="wrap">', 1)
standalone = ('<!doctype html>\n<html lang="en">\n<head>\n<meta charset="utf-8">\n'
              '<meta name="viewport" content="width=device-width, initial-scale=1">\n'
              + head + '</head>\n<body>\n<div class="wrap">'
              + body.replace(cdn, "<script>/* plotly.js-dist-min 2.35.2, MIT licence, inlined so this file works offline */\n" + lib + "\n</script>")
              + "\n</body>\n</html>\n")
open(os.path.join(WORK, "chorismate_optimal_charges.html"), "w").write(standalone)
print("cdn page %d bytes, standalone %d bytes" % (len(page), len(standalone)))
