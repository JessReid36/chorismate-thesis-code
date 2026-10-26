#!/usr/bin/env python3
"""Citation integrity check for the chorismate thesis record.

Two questions, both answered mechanically against the committed reference library:

  (1) PRESENCE.  Every DOI cited anywhere in the results repo - is the paper actually in
      chorismate-refs? The project's standing rule is that a claim must come from the full paper
      in that library, not from an abstract, a search result or memory. A cited DOI with no paper
      behind it marks a claim that cannot meet that rule.

  (2) QUOTATION.  Every phrase the documents present as a quotation from a paper - does it
      actually occur in the library?

WHY (2) IS NOT A PLAIN grep, AND WHY THIS SCRIPT EXISTS AT ALL. On 26 October 2026 a plain
substring search over the extracted text of the Dittner dissertation reported that a quotation did
not occur in it. STAGE6_ERRATUM7_20261026.txt records what followed: the quotation was real, the
search failed because the PDF sets "specifically" with the fi ligature U+FB01, and a committed
document went on to accuse an earlier session of fabricating it. Three further quotations were
reported unverifiable the same evening for two more reasons of the same kind. Extracted PDF text
defeats naive matching in at least four ways, all of them present in this library:

    ligatures          speci<U+FB01>cally, e<U+FB00>ect, di<U+FB03>cult
    OCR substitution   Claeyssens renders "mol-1" as "mol21" and "=" as "5";
                       Galbraith renders the minus sign as the digit 2, so -0.378 reads 20.378
    two-column text    a sentence's second half can be hundreds of lines from its first, with
                       the other column's text interleaved between them
    soft hyphenation   a word split across a line break as "regulari-" / "zation"

So the quotation check below never relies on one method. It reports a strict normalised substring
match AND a 5-gram coverage fraction, and anything it cannot match it reports as UNMATCHED rather
than as absent - because this script's own history is the reason to distrust a negative.

NO NETWORK. Read-only. Pure standard library apart from pdftotext, which is used when present.

usage: cs_citecheck.py <results-repo> <refs-dir> <out-dir>
"""
import sys, os, re, subprocess, shutil, tempfile, hashlib, unicodedata
from pathlib import Path

LIG = {"ﬀ": "ff", "ﬁ": "fi", "ﬂ": "fl", "ﬃ": "ffi", "ﬄ": "ffl",
       "‘": "'", "’": "'", "“": '"', "”": '"',
       "–": "-", "—": "-", "−": "-", " ": " "}
SHINGLE = 4
STRONG, WEAK = 0.80, 0.40
# An elided quotation - "A ... B" - can never match contiguously, and its shingles are broken at
# the elision. Each fragment is therefore verified separately and the quotation is classified by
# its WORST fragment, so an ellipsis cannot be used to smuggle an unsupported half past the check.
RE_ELIDE = re.compile(r"\s*(?:\.\.\.|\u2026)\s*")


def stop(m):
    print("STOP: %s" % m, file=sys.stderr)
    raise SystemExit(1)


def norm(s):
    for k, v in LIG.items():
        s = s.replace(k, v)
    s = unicodedata.normalize("NFKD", s)
    s = "".join(c for c in s if not unicodedata.combining(c))
    s = re.sub(r"(\w)-\s*\n\s*(\w)", r"\1\2", s)        # soft hyphenation across a line break
    return re.sub(r"\s+", " ", s)


def squash(s):
    return re.sub(r"[^a-z0-9]", "", norm(s).lower())


def shingles(s, n=SHINGLE):
    w = re.findall(r"[a-z0-9']+", norm(s).lower())
    return {" ".join(w[i:i + n]) for i in range(max(0, len(w) - n + 1))}


def extract(refs, cache):
    """Text of every PDF in the reference library. pdftotext when available; otherwise names only."""
    pdfs = sorted(p for p in Path(refs).rglob("*.pdf") if ".git" not in p.parts)
    if not pdfs:
        stop("no PDFs found under %s - is that the chorismate-refs checkout?" % refs)
    have = shutil.which("pdftotext") is not None
    out = {}
    for p in pdfs:
        if not have:
            out[p.name] = ""
            continue
        t = Path(cache) / (re.sub(r"[^A-Za-z0-9._-]", "_", p.stem) + ".txt")
        if not t.exists():
            subprocess.run(["pdftotext", "-layout", str(p), str(t)],
                           stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        try:
            out[p.name] = t.read_text(errors="replace") if t.exists() else ""
        except Exception:
            out[p.name] = ""
    return out, have, len(pdfs)


RE_DOI = re.compile(r"10\.[0-9]{4,9}/[A-Za-z0-9._;()/:+<>-]+")
RE_QUOTE = re.compile("[\"\u201c]([^\"\u201c\u201d]{30,220}?)[\"\u201d]")
RE_CITE = re.compile(r"(doi:|10\.\d{4}/|et al\.|\b(19|20)\d\d\b\s*,|\bp\.\s*\d|§\s*\d"
                     r"|thesis|dissertation|writes|states|reads|quoted|verbatim|report[s]?\b)", re.I)


def docs(res):
    out = {}
    for p in sorted(list(Path(res).rglob("*.md")) + list(Path(res).rglob("*.txt"))):
        if ".git" in p.parts:
            continue
        try:
            out[str(p.relative_to(res))] = p.read_text(errors="replace")
        except Exception:
            pass
    return out


def dois(D):
    cited = {}
    for name, t in D.items():
        for m in RE_DOI.finditer(t):
            d = m.group(0).rstrip(".,;)").lower()
            if len(d.split("/", 1)[1]) >= 4:
                cited.setdefault(d, set()).add(name)
    return cited


def present(doi, lib, libsq):
    suf = squash(doi.split("/", 1)[1])
    if len(suf) < 5:
        return None
    for n in lib:
        if suf in squash(n):
            return n
    for n, t in libsq.items():
        if t and doi in t:
            return n + " (in text)"
    return None


def quotes(D):
    """Quotes are matched on a FLATTENED copy of each document: markdown blockquote markers
    stripped and newlines collapsed. A quotation set as a multi-line block quote - which is how
    the Dittner passage of 26 October was set - is invisible to a line-bounded regex, and that
    quotation is the reason this script exists."""
    out = []
    for name, t in D.items():
        flat = norm(re.sub(r"(?m)^[ \t]*>[ \t]?", "", t))
        for m in RE_QUOTE.finditer(flat):
            q = m.group(1).strip()
            if len(q.split()) < 6:
                continue
            ctx = flat[max(0, m.start() - 220): m.end() + 220]
            if RE_CITE.search(ctx):
                out.append((name, q))
    return out


def cover(frag, index):
    """Best 4-gram coverage of one fragment over any file in the index."""
    h = shingles(frag)
    if not h:
        return 0.0, None
    best, bk = 0.0, None
    for n, S in index.items():
        f = len(h & S) / len(h)
        if f > best:
            best, bk = f, n
    return best, bk


def verify(q, libsq, libsh, ownsq, ownsh, src=None):
    """Strict first, tolerant second, own-documents last. Never reports 'absent' on one method."""
    qs = squash(q)
    for n, t in libsq.items():
        if t and qs in t:
            return "VERBATIM", n, 1.0
    frags = [f for f in RE_ELIDE.split(q) if len(f.split()) >= 4] or [q]
    scores = [cover(f, libsh) for f in frags]
    best = min(x[0] for x in scores)
    bk = max(scores, key=lambda x: x[0])[1]
    if best >= STRONG:
        return "VERIFIED", bk, best
    own = {n: t for n, t in ownsq.items() if n != src}   # never match the document it came from
    for n, t in own.items():
        if qs in t:
            return "INTERNAL", n, 1.0
    oidx = {n: S for n, S in ownsh.items() if n != src}
    oscores = [cover(f, oidx) for f in frags]
    ob = min(x[0] for x in oscores)
    if ob >= STRONG:
        return "INTERNAL", max(oscores, key=lambda x: x[0])[1], ob
    if best >= WEAK:
        return "PARTIAL", bk, best
    return "UNMATCHED", bk, best


def main(res, refs, out):
    res, refs, out = Path(res), Path(refs), Path(out)
    if str(out.resolve()).startswith(str(res.resolve())):
        stop("refusing to write inside the results repo: %s" % out)
    if not (res / ".git").is_dir():
        stop("%s is not a git repo" % res)
    cache = tempfile.mkdtemp(prefix="citecheck.")
    print("   extracting the reference library ...", flush=True)
    lib, have_pdftotext, npdf = extract(refs, cache)
    libsq = {n: squash(t) for n, t in lib.items()}
    libsh = {n: shingles(t) for n, t in lib.items()}

    print("   reading the committed documents ...", flush=True)
    D = docs(res)
    ownsq = {n: squash(t) for n, t in D.items()}
    ownsh = {n: shingles(t) for n, t in D.items()}

    C = dois(D)
    rows = [(d, present(d, lib, libsq), sorted(C[d])) for d in sorted(C)]
    have = [r for r in rows if r[1]]
    miss = [r for r in rows if not r[1]]

    print("   checking quotations ...", flush=True)
    Q = quotes(D)
    res_q = [(n, q) + verify(q, libsq, libsh, ownsq, ownsh, src=n) for n, q in Q]
    by = lambda k: [r for r in res_q if r[2] == k]

    O = ["# Citation integrity, %s" % res.name,
         "",
         "Produced by `cs_citecheck.py`. NO NETWORK, read-only, nothing inferred. Re-run it whenever a",
         "paper is added to the reference library and the counts below should move.",
         "",
         "Reference library: **%d PDFs**%s" % (npdf, "" if have_pdftotext
                                               else "  — **pdftotext was NOT available, so quotation checking is disabled and the DOI check is filename-only. Install poppler-utils and re-run.**"),
         "Committed documents scanned: **%d**" % len(D),
         "",
         "## 1. Cited DOIs, and whether the paper is in the library",
         "",
         "The project's rule is that a claim must come from the full paper in `chorismate-refs`. A DOI",
         "cited with no paper behind it marks a claim that cannot currently meet that rule. It does not",
         "mean the claim is wrong — it means it is unverifiable here.",
         "",
         "**%d distinct DOIs cited; %d present, %d absent.**" % (len(rows), len(have), len(miss)),
         "",
         "### Absent from the library",
         "",
         "| DOI | cited in |",
         "|---|---|"]
    for d, _, src in miss:
        O.append("| `%s` | %s |" % (d, ", ".join(s.split("/")[-1] for s in src[:3])
                                    + (" …" if len(src) > 3 else "")))
    O += ["", "### Present in the library", "", "| DOI | file |", "|---|---|"]
    for d, f, _ in have:
        O.append("| `%s` | %s |" % (d, f))

    O += ["", "## 2. Quotations attributed to a paper", "",
          "Four categories. **VERBATIM** matched as an exact normalised substring. **VERIFIED** matched",
          "at or above %.0f%% of its %d-word shingles, every elided fragment checked separately — the" % (100 * STRONG, SHINGLE),
          "normal result for a real quotation broken by two-column extraction or OCR. **INTERNAL** matched",
          "another committed document, i.e. the",
          "project quoting itself, not a paper. **PARTIAL** and **UNMATCHED** are reported for inspection",
          "and are NOT claims of fabrication: see this script's docstring for why that distinction is",
          "the whole reason it exists.",
          "",
          "**%d attributed quotations: %d verbatim, %d verified, %d internal, %d partial, %d unmatched.**"
          % (len(res_q), len(by("VERBATIM")), len(by("VERIFIED")), len(by("INTERNAL")),
             len(by("PARTIAL")), len(by("UNMATCHED"))),
          ""]
    for k, head in (("UNMATCHED", "Unmatched — inspect each by hand against the PDF"),
                    ("PARTIAL", "Partial — most will be OCR or line-break damage, not error")):
        g = by(k)
        if not g:
            O += ["### %s" % head, "", "None.", ""]
            continue
        O += ["### %s" % head, "", "| document | quotation | closest file | shingle match |", "|---|---|---|---|"]
        for n, q, _, bk, f in sorted(g, key=lambda r: -r[4]):
            O.append("| %s | %s | %s | %.0f%% |"
                     % (n.split("/")[-1], q[:90].replace("|", "/"), bk or "—", 100 * f))
        O.append("")
    O += ["## 3. What this cannot tell you", "",
          "It cannot tell you a quotation is fabricated. It can only say this library does not contain it,",
          "which is a different statement — and on 26 October 2026 that difference was got wrong in a",
          "committed document. An unmatched quotation may come from any of the absent papers in section 1.",
          "It also cannot check paraphrase, a number quoted without quotation marks, or whether a paper",
          "that IS present actually supports the claim made from it.", ""]
    out.mkdir(parents=True, exist_ok=True)
    (out / "CITATION_INTEGRITY_20261026.md").write_text("\n".join(O) + "\n")
    print("\n".join(O[:22]))
    print("   ...")
    print("   wrote %s" % (out / "CITATION_INTEGRITY_20261026.md"))
    shutil.rmtree(cache, ignore_errors=True)


if __name__ == "__main__":
    if len(sys.argv) != 4:
        stop("usage: cs_citecheck.py <results-repo> <refs-dir> <out-dir>")
    main(sys.argv[1], sys.argv[2], sys.argv[3])
