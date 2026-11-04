# cs_citecheck.py's in-text fallback has never fired

4 November 2026. Hand-written. The measurement in section 3 is reproducible with the snippet given.

## 1. The bug

`present()` decides whether a cited DOI has a paper behind it:

```python
def present(doi, lib, libsq):
    suf = squash(doi.split("/", 1)[1])
    if len(suf) < 5:
        return None
    for n in lib:                      # (a) filename match
        if suf in squash(n):
            return n
    for n, t in libsq.items():         # (b) fallback: DOI in the paper's text
        if t and doi in t:
            return n + " (in text)"
    return None
```

`libsq` is built in `main()` as `{n: squash(t) for n, t in lib.items()}` - the SQUASHED text, with
everything outside `[a-z0-9]` removed. Branch (b) tests a RAW DOI, which always contains `.` and
`/`, against that squashed text. A raw DOI can never be a substring of it. **Branch (b) is dead
code and has never returned a match.**

So every presence verdict the tool has ever produced was decided by branch (a), the filename, alone
- which is the behaviour the report explicitly describes as the DEGRADED mode for when pdftotext is
missing:

> "pdftotext was NOT available, so quotation checking is disabled and the DOI check is
> filename-only."

pdftotext has been available. The DOI check was filename-only anyway.

## 2. The fix

One line. Squash the DOI on the same terms as the text it is compared against:

```python
        if t and squash(doi) in t:
```

## 3. What it changes, measured

Against the library at 79 PDFs and the 69 distinct DOIs cited across 1494 committed documents:

    present under the current present()   34
    present under the fixed present()     40
    newly found, were reported ABSENT      6

The six, with the file that holds each:

    10.1002/wcms.70019                          WIREs ... Neese - ORCA 6.0.pdf        (in text)
    10.1016/s1574-1400(09)00502-7               S1574-1400%2809%2900502-7.pdf          (in text)
    10.1021/acs.jctc.4c01698                    ct4c01698.pdf                          (in text)
    10.1039/a909486e                            v09-092.pdf                           (in text)
    10.1146/annurev.physchem.59.032607.093618   v09-092.pdf                           (in text)
    10.33011/livecoms.1.1.5067                  nihms-994026.pdf                       (in text)

All six are papers the library already held. All six were reported as claims that "cannot currently
meet" the project's sourcing rule. They met it.

Reproduce with:

    cd phase2.2 && python3 - <<'EOF'
    import sys,re,subprocess,shutil,tempfile,hashlib,unicodedata
    from pathlib import Path
    exec(open('cs_citecheck.py').read().split('def main(')[0])
    refs=Path('../../chorismate-refs'); res=Path('../../chorismate-thesis-results')
    lib,_,n=extract(refs,tempfile.mkdtemp()); libsq={k:squash(v) for k,v in lib.items()}
    C=dois(docs(res))
    f=lambda doi:(lambda s:(next((k for k in lib if s in squash(k)),None) or
        next((k+' (in text)' for k,t in libsq.items() if t and squash(doi) in t),None)))(
        squash(doi.split('/',1)[1]))
    cur=[d for d in sorted(C) if present(d,lib,libsq)]; fix=[d for d in sorted(C) if f(d)]
    print(len(cur),len(fix)); [print(d,'->',f(d)) for d in fix if d not in cur]
    EOF

## 4. The second blind spot, which the fix does NOT close

Branch (a) matches on filename, so a paper filed under anything other than its own DOI suffix is
invisible to it. Two confirmed in this library:

    nihms-1716744.pdf   is Hennefarth & Alexandrova, Curr Opin Struct Biol 2022, 72, 1-8,
                        doi 10.1016/j.sbi.2021.06.006   (HHS author manuscript)
    au3c00536.pdf       carries its own doi 10.1021/jacsau.3c00536; squash("jacsau3c00536") is not
                        a substring of squash("au3c00536.pdf")

With the fix of section 2 both become matchable via their text, because both stamp their own DOI on
page 1. Without it, neither can ever be found. Author manuscripts from PMC and PubMed Central are
the general case and the library has at least three (`nihms-*.pdf`).

`name_refs.py`, added alongside, reads a PDF's own page-1 DOI and reports the filename that branch
(a) will match, so new papers are filed under a name the tool can see. It also reports when the
library already holds the paper under another name, which is how the `nihms-1716744.pdf` case above
was found.

## 5. Consequence for the record

`CITATION_INTEGRITY_20261026.md` over-reports absence. Its "37 absent" includes at least six papers
the library held, and its counts are in any case stale: it scanned 77 PDFs and the library now holds
79, and it lists `10.1021/ja049376t` as absent when `ja049376t.pdf` is present.

Nothing in it over-reports PRESENCE, so no claim that the report called sourced was unsourced. The
error is in the safe direction, but it has been pointing at work as unverifiable when it was not.

## 6. One cosmetic issue, not fixed here

`main()` writes a hardcoded output name:

    (out / "CITATION_INTEGRITY_20261026.md").write_text(...)

so a re-run on any later date produces a file named for 26 October. Left alone deliberately -
changing it may break whatever reads that path. Rename on copy, or change the literal to a date
argument, as preferred.
