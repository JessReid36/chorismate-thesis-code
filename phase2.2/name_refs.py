#!/usr/bin/env python3
"""name_refs.py - identify newly downloaded PDFs, warn if the library already has
them, and give each the filename cs_citecheck.py's present() will match.

Reads each PDF's OWN text. Nothing is trusted from outside the files.

usage: python3 name_refs.py <refs-dir> <folder-or-pdf> [more ...]
   eg: python3 name_refs.py ~/Desktop/chorismate-refs ~/Downloads/newrefs

Renames nothing. Prints a report and the mv commands for you to check.
"""
import re, subprocess, sys, unicodedata
from pathlib import Path

LIG = {"ﬀ":"ff","ﬁ":"fi","ﬂ":"fl","ﬃ":"ffi","ﬄ":"ffl",
       "‘":"'","’":"'","“":'"',"”":'"',
       "–":"-","—":"-","−":"-"," ":" "}
def norm(s):
    for k,v in LIG.items(): s=s.replace(k,v)
    s=unicodedata.normalize("NFKD",s)
    return "".join(c for c in s if not unicodedata.combining(c))
def squash(s): return re.sub(r"[^a-z0-9]","",norm(s).lower())

RE_DOI = re.compile(r"10\.[0-9]{4,9}/[A-Za-z0-9._;()/:+<>-]+")

def text(p, page1=False):
    cmd=["pdftotext"]+(["-f","1","-l","1"] if page1 else [])+[str(p),"-"]
    try: return subprocess.run(cmd,capture_output=True,text=True,timeout=180).stdout
    except FileNotFoundError:
        sys.exit("STOP: pdftotext not found.  sudo apt install poppler-utils")
    except subprocess.TimeoutExpired: return ""

def clean(d):
    d=d.rstrip(".,;:)>]'\"")
    for tail in ("Thisarticle","Downloaded","http","www","Published","Received",
                 "Accepted","Cite","Author","Available"):
        i=d.find(tail)
        if i>0: d=d[:i]
    return d.rstrip(".,;:)>]'\"")

def dois_in(t):
    out=[]
    for m in RE_DOI.finditer(t):
        d=clean(m.group(0))
        if len(d.split("/",1)[-1])>=5 and d not in out: out.append(d)
    return out

def titleish(page1):
    """Longest plausible title line from page 1, used only for duplicate detection."""
    lines=[l.strip() for l in page1.splitlines() if l.strip()]
    cand=[l for l in lines[:25] if 25<=len(l)<=170
          and not re.search(r"@|^\s*(doi|https?|www|received|accepted|published|abstract|"
                            r"keywords|author manuscript|hhs public access|subscriber)",l,re.I)]
    return max(cand,key=len) if cand else ""

def main(refs, targets):
    refs=Path(refs).expanduser()
    if not refs.is_dir(): sys.exit("STOP: refs dir not found: %s" % refs)
    print("indexing the reference library at %s ..." % refs, flush=True)
    libtxt={p.name: squash(text(p)) for p in sorted(refs.glob("*.pdf"))}
    libnames=list(libtxt)
    print("  %d PDFs indexed\n" % len(libnames))

    pdfs=[]
    for f in targets:
        d=Path(f).expanduser()
        if d.is_file() and d.suffix.lower()==".pdf": pdfs.append(d)
        elif d.is_dir(): pdfs+=sorted(d.glob("*.pdf"))
        else: print("skipping (not a pdf or folder): %s" % d)
    if not pdfs: sys.exit("no PDFs found")

    moves=[]; dups=[]
    for p in pdfs:
        print("="*78); print("FILE: %s" % p)
        p1=text(p,True); full=text(p)
        for l in [l.strip() for l in p1.splitlines() if l.strip()][:8]:
            print("     %s" % l[:94])
        if not full.strip():
            print("  !! no extractable text (scanned?). Handle by hand."); continue

        # --- already in the library? ---
        t=titleish(p1); hit=None
        if len(squash(t))>=25:
            for n in libnames:
                if libtxt[n] and squash(t) in libtxt[n]: hit=n; break
        if hit:
            print("  ** ALREADY IN THE LIBRARY as %s" % hit)
            print("     matched on title: %s" % t[:88])
            dups.append((p,hit))

        # --- the paper's own DOI ---
        d1=dois_in(p1); dall=dois_in(full)
        own = d1[0] if d1 else (dall[0] if dall else None)
        where = "page 1" if d1 else ("later in the text" if dall else None)
        if not own:
            print("  !! no DOI in the file. Name it by hand."); continue
        print("  -- DOIs on page 1      : %s" % (", ".join(d1[:4]) or "none"))
        print("  -- own DOI taken as    : %s   (found on %s)" % (own,where))
        if not d1:
            print("     WARNING: not on page 1, so this may be a reference, not this paper's DOI.")
        suf=own.split("/",1)[1]; target="%s.pdf" % suf
        ok=squash(suf) in squash(target)
        print("  -- filename to use     : %s" % target)
        print("  -- present() matches   : %s" % ("YES" if ok else "NO - fix by hand"))
        clash = target in libnames
        if clash: print("     NOTE: %s already exists in the library." % target)
        if hit:
            print("     -> library already has this paper; a correctly-named copy still helps"
                  " cs_citecheck.py match it by filename.")
        if ok and not clash: moves.append((p,target))
        elif ok and clash: print("     -> skipping the mv: target name already taken.")

    print("\n"+"="*78)
    if dups:
        print("DUPLICATES (content already in the library):")
        for p,n in dups: print("   %s  ==  %s" % (p.name,n))
        print()
    if moves:
        print("RENAME + MOVE (check each target above first):\n")
        print('REFS="%s"' % refs)
        for p,t in moves: print('cp -i -- "%s" "$REFS/%s"' % (p,t))
    else:
        print("nothing to copy.")

if __name__=="__main__":
    if len(sys.argv)<3: sys.exit(__doc__)
    main(sys.argv[1], sys.argv[2:])
