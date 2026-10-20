"""build.py - the site-effect map pages from template.html (beside this script) and <work_dir>/data.json.
Usage: python3 build.py <work_dir> <plotly.min.js>
Writes <work_dir>/site_effect_map.html (Plotly from the jsDelivr CDN) and <work_dir>/chorismate_site_effect_map.html
(Plotly inlined, so it opens offline). The inlined library must be plotly.js-dist-min 2.35.2 (npm; MIT licence)."""
import hashlib, os, sys
WORK, LIB = sys.argv[1], sys.argv[2]
HERE = os.path.dirname(os.path.abspath(__file__))
PLOTLY_SHA = "6d21266ce1bd7d9e5ab4e115989c70c20de0382fd973a8f26ab58619eba4d603"   # plotly.js-dist-min 2.35.2, plotly.min.js
t = open(os.path.join(HERE, "template.html")).read(); data = open(os.path.join(WORK, "data.json")).read()
page = t.replace("__DATA__", data)
assert "__RESCOMMIT__" not in page and "__SUBDIFF__" not in page, "fill the residue note's placeholders from the real x89 data first"
open(os.path.join(WORK, "site_effect_map.html"), "w").write(page)                       # artifact (Plotly from the CDN)
cdn = '<script src="https://cdn.jsdelivr.net/npm/plotly.js-dist-min@2.35.2/plotly.min.js"></script>'
assert page.count(cdn) == 1
assert hashlib.sha256(open(LIB, "rb").read()).hexdigest() == PLOTLY_SHA, LIB + " is not plotly.js-dist-min 2.35.2's plotly.min.js"
lib = open(LIB).read()
assert "</script" not in lib
head, body = page.split('<div class="wrap">', 1)
standalone = ("<!doctype html>\n<html lang=\"en\">\n<head>\n<meta charset=\"utf-8\">\n<meta name=\"viewport\" content=\"width=device-width, initial-scale=1\">\n"
              + head + "</head>\n<body>\n<div class=\"wrap\">" + body.replace(cdn, "<script>/* plotly.js-dist-min 2.35.2, MIT licence, inlined so this file works offline */\n" + lib + "\n</script>") + "\n</body>\n</html>\n")
open(os.path.join(WORK, "chorismate_site_effect_map.html"), "w").write(standalone)
print(len(page), len(standalone))
