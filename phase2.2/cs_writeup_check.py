#!/usr/bin/env python3
"""Validation of WRITEUP_CONSOLIDATED_<date>.md against the two repositories it is written from.

    python3 cs_writeup_check.py --res <results repo> --code <code repo> --doc <the .md> --figures <its figures dir>
                                [--full] [--selftest] [--provenance <file>]

It is run twice before anything is committed: once with the document and its figures in a staging directory
outside the repositories, and once after they have been placed at their final paths.  Both runs must pass.

THE TEN CHECKS (the numbering is the brief's)

  1 numbers      Every numeral in the document is traceable.  The document is cut into blocks at blank lines.  A block
                 carries tags: [src: KEY, ...] names committed files (keys resolved by the table of Appendix C) and
                 [dv: DVnn, ...] names rows of figures/derived.tsv.  Each numeral in the block must either occur, as a
                 whole number and not as part of a longer one, in one of the block's [src:] files (class SOURCED), or
                 equal the value of one of its [dv:] rows (class DERIVED), or stand in a block marked (estimated)
                 (class ESTIMATED).  Anything else fails.  NOT treated as numerals, and counted as such: dates;
                 headings; list numbering; cross-references ("section 4.6", "Table 1", "Stage 4e", "Amendment 5");
                 identifiers that join letters and whole numbers (D33, S6-20, def2-SVP); text in back-ticks, which is
                 checked as a path or a name instead; and fenced blocks.  An identifier that contains a decimal number
                 (5.6433e_n20) is not exempt: it must occur verbatim in a [src:] file.
  2 quotes       Every passage in curly quotation marks occurs verbatim, whitespace-normalised, in a [src:] file of
                 its own block.
  3 citations    Every file in Appendix C is tracked by git (asked of git ls-files) or, for the write-up's own figures
                 directory, exists there; every key used in a tag is defined; every back-ticked path is tracked.
  4 withdrawn    Each withdrawn figure in WITHDRAWN below appears only inside the section headed "withdrawn".  Each
                 registry entry carries the sentence of the committed erratum that withdraws it, and that sentence is
                 itself checked against the erratum, so the registry cannot drift from the record.
  5 selflabels   Where sections 1 to 11 use a result that labels itself PARTIAL, INCONCLUSIVE, EXPLORATORY, post hoc or
                 NOT ADOPTED, the label stands in the same section.  The entries are SELF_LABELS below, each with the
                 committed sentence in which the result labels itself.  The appendices list what was read and are not
                 uses, so they are outside this check.
  6 controls     The production reference barrier, every ddE of frozen/frozen.tsv and D33 stage A's corrected reference
                 are recomputed from raw committed outputs (figures/make_figures.py, through the committed readers)
                 and must match the record; the worst disagreement is reported.
  7 counts       Every run count in figures/derived.tsv equals an independent count of that diagnostic's outputs, made
                 here from git ls-files and not from the file system.
  8 not_tested   Every name between the document's not-tested markers stands on no ORCA keyword line of any input or
                 output echo in either repository.
  9 regeneration Every figure has .svg, .png and .tsv; and make_figures.py and sweep_pass1.py, re-run here, reproduce
                 every .tsv byte for byte.  sweep_pass2.py takes about four minutes and is re-run only with --full.
 10 heads        The two HEAD hashes in the document are the repositories' HEADs (or, once the write-up is committed,
                 ancestors of them with nothing but the write-up's own files changed since), and neither working tree
                 holds anything other than the write-up's own files.

SIX CHECKS ADDED ON 11 NOVEMBER 2026, when the write-up was revised for DIAGNOSTIC 34

 11 exp_heads    The write-up was rebuilt from results 0ede870 and code ee28301 (the revision of 18 November 2026 for
                 DIAGNOSTIC 43, AMENDMENTS 6 and 7 and the pre-registration of DIAGNOSTIC 44; before that 6c09d9c and
                 9f66ee3, on 16 November 2026; before that 21d631b and 512eead, on the evening of 15 November 2026;
                 before that 4af2afa and 1e47614, and before that f461bd1 and 23c2500).  The document states both, and each
                 repository's HEAD is that commit (or, once the write-up is committed, a descendant of it with nothing
                 but the write-up's own files changed since).
 12 d34_numbers  A block that cites DIAGNOSTIC 34 (a [src:] tag naming its result, its erratum, d34.tsv, or the tables
                 T3 and T4 recomputed from its outputs) may hold only numerals that occur in
                 STAGE6_DIAGNOSTIC34_RESULT_20261110.txt or d34.tsv, or that are [dv:] rows.  A [dv:] row is arithmetic
                 done by figures/make_figures.py, and each one used in such a block must come out the same when that
                 script is re-run here.  A number that stands only in a table of the write-up's own does not pass.
 13 d34_controls The control test, recomputed here from the raw outputs of d34/ and of frozen/sp_outputs.tar.xz through
                 the committed readers, by code in this file and not by make_figures.py: both arms' control barriers
                 are 26.8248, and every ddE of frozen/frozen.tsv is reproduced.  figures/controls.tsv must say the same.
 14 arm1_legs    Arm 1's legs of the criterion (Loewdin rise and dipole ratio in each state, whether a level shift or
                 TRAH was engaged, and ddE) are recomputed here from d34/a1_* and must equal what the document reports
                 of them: the a1 columns of figures/T3_d34_arms.tsv, the B1 rows of d34.tsv, and the [dv:] rows
                 that count or rank them.
 15 claims       The revision brief carried six claims from another session, C1 to C6, each unverified.  CLAIMS below
                 lists, for each, the [dv:] rows that hold the numbers behind it.  Every one of those rows is recomputed
                 here from committed raw files and must equal figures/derived.tsv; C4's four numbers are recomputed
                 from the raw outputs of standoff/ and must stand in DIAGNOSTIC 11's result.  The forms of C6 that the
                 recomputation rejected must not appear, and C3, which is post hoc, must not appear in the Conclusions.
 16 posthoc      Every block that discusses C3 (the sign or the size of the positive charge against relocation) carries
                 the words "post hoc".

FOUR CHECKS ADDED ON 15 NOVEMBER 2026, when the write-up was revised for DIAGNOSTICS 35 to 38

 17 benchmark    DIAGNOSTIC 36's DLPNO-CCSD(T) barrier and its three parts are recomputed here from the raw outputs of
                 d36a2/, by code in this file and not by make_figures.py, and must equal the [dv:] rows that hold them,
                 the barrier of d36.tsv, and the number the document states.  DIAGNOSTIC 38's two energy shifts for its
                 interpretable run are recomputed from d38/ in the same way.
 18 experiment   ERRATUM 7, E7.11: a computed barrier may stand beside an experimental one only if the sentence says
                 which quantity is which.  Every block that holds the experimental 24.5 says that it is a free energy,
                 and every block that holds it together with the correlated barrier also says that the two are not
                 the same quantity or that the closeness is not agreement.
 19 d38_a1       Run A1 of DIAGNOSTIC 38 is NOT INTERPRETABLE by that diagnostic's pre-registered rule.  Every block
                 that quotes one of A1's own numbers carries those words.
 20 posthoc_lvl  DIAGNOSTIC 36 labels the comparison of DIAGNOSTIC 37's five functionals with its benchmark POST HOC:
                 every block that quotes one of those five distances, or the spread and position built on them,
                 carries the words "post hoc".  The orbital energies of charged designs (figures/T9) are tabulated
                 as data only, because the record holds that an orbital energy is not a test of binding: every
                 block that quotes one carries "post hoc" and either "data only" or "not the test".

Check 8 was also widened on 15 November 2026: a name is refused on the not-tested list if any ORCA keyword that
requests it under another spelling (M062X for M06-2X, RI-PWPB95 for PWPB95, ...) stands on a keyword line.

FOUR CHECKS ADDED ON THE EVENING OF 15 NOVEMBER 2026, when the write-up was revised for DIAGNOSTICS 40 and 41 and
for DIAGNOSTIC 37's Amendment 1

 21 benchmark_sv DIAGNOSTIC 41's DLPNO-CCSD(T)/def2-SVPD barrier and its three parts are recomputed here from the raw
                 outputs of d38/N (the reactant leg, as committed) and d41/T, by code in this file, and must equal the
                 [dv:] rows that hold them, the barrier of d41.tsv and the number the document states.  The shift from
                 def2-SVPD to def2-TZVPD is recomputed from those two and DIAGNOSTIC 36's two, and must equal its row.
 22 detachment   Every detachment energy of DIAGNOSTIC 40 is recomputed here from the raw outputs: the monoanion of
                 d40/<run> minus the committed dianion run that d40.tsv names as its source, each read by code in this
                 file (the (QM/MM) total for a QM/MM run, the last final energy otherwise).  Each must equal d40.tsv and
                 the row of figures/T12 or T13 that holds it, and in T13 the word must be the one the sign gives:
                 BOUND above zero, METASTABLE at or below it.
 23 metastable   DIAGNOSTIC 40 names five lowerings that Probe 2's sentence reaches: a ddE measured where the dianion
                 is METASTABLE "is not a barrier effect and must not be reported as one".  D40 withdraws none of them,
                 so check 4 does not apply; until an erratum decides, every block of the document that quotes one of
                 the five carries the word METASTABLE.
 24 d39_held     DIAGNOSTIC 39's result is committed and is NOT ENTERED in this revision, at the author's instruction.
                 The document says so; no raw output or table of D39 is a source; no table of the figures directory
                 other than the file census and the two run counts reads a D39 output; and a block that cites D39's
                 result or amendment may hold no number that is not in another of its sources or a [dv:] row, so that
                 no number of D39 can be quoted under cover of the citation.

Two older checks changed with them.  Check 5: the self-label of revDSD-PBEP86-D4 was NOT MEASURED, which is DIAGNOSTIC
37's label for it; Amendment 1 has measured it, so the label that must stand wherever it is named is now Amendment 1.
Check 20: the def2-SVPD distance of revDSD-PBEP86-D4 from DIAGNOSTIC 41's benchmark, and the distances of def2-SVPD
barriers from the def2-TZVPD benchmark, are this write-up's own comparisons, and a block that quotes one of those
used here carries the words "post hoc".

TWO CHECKS ADDED ON 16 NOVEMBER 2026, when the write-up was revised for DIAGNOSTIC 42

 25 bound_d42    Every detachment energy of DIAGNOSTIC 42 is recomputed here from the raw outputs: the monoanion of
                 d42/<run> minus the committed dianion run of DIAGNOSTIC 31 that d42.tsv names as its source, each
                 read by code in this file.  Each must equal d42.tsv and the row of figures/T14 that holds it, and the
                 word must be the one DIAGNOSTIC 42's rule gives: METASTABLE at or below zero; BOUND above zero if the
                 monoanion's last stability verdict is stable.  The ddE of each condition is recomputed from DIAGNOSTIC
                 31's committed dianion outputs and its two charge-free controls and must equal figures/T14 and
                 d42.tsv; the reading of Probe 2's sentence in figures/T14 must be the one the two words give; the
                 number of conditions BOUND at both states must be its [dv:] row; and the verdict on H9 that the
                 document quotes must be the one those words give.
 26 not_quotable DIAGNOSTIC 42 computed a ddE for each of its thirteen conditions and says "no ddE listed here is made
                 quotable by D42".  Every block of the document that quotes one of the thirteen, as recomputed by
                 check 25, carries the words NOT QUOTABLE, and a block that quotes one whose condition has a METASTABLE
                 state carries the word METASTABLE as well.  And the comparisons of DIAGNOSTIC 42's numbers that are
                 this write-up's own (how thin the margins are beside the spread between functionals; the change per
                 unit of charge) are in [dv:] rows marked POST HOC, and a block that quotes one of those values carries
                 the words "post hoc".

TWO CHECKS ADDED ON 18 NOVEMBER 2026, when the write-up was revised for DIAGNOSTIC 43 and AMENDMENTS 6 and 7

 27 survivors_d43 DIAGNOSTIC 43 is recomputed here from the raw outputs of d43/ and the two charge-free controls of
                 DIAGNOSTIC 31.  For each of its 36 condition-states: the detachment energy, the monoanion run minus the
                 dianion run that d43.tsv names, with the word AMENDMENT 6 gives (METASTABLE at or below zero; BOUND
                 above zero if the monoanion's last stability verdict is stable); and the four legs of the dianion
                 through the committed readers (cs6_audit.read, cs6_d33b.stability and SHIFT) at the committed
                 thresholds.  For each condition: the four-leg verdict, whether the fifth leg is passed, and ddE.  Each
                 must equal d43.tsv and the row of figures/T15 that holds it.  The survivors, the designs of part B that
                 are RELIABLE and BOUND at both states, must be the ones figures/T15 and d43.tsv mark, and their number
                 must be its [dv:] row.  The verdicts on H3, H2 and H9 that the document quotes must be the ones those
                 readings give.  And the document's statement that DIAGNOSTIC 44 has no output must be true of the
                 repository: the number of job.out files under d44/ must be the [dv:] row that says so, and if any
                 exists the document may not say that none does.
 28 working_level AMENDMENT 7: "No ddE of a survivor is quoted as more than a working-level number until it has been
                 through A7.3."  Every block of the document that quotes the ddE of a condition of DIAGNOSTIC 43 that is
                 RELIABLE and BOUND at both states, as recomputed by check 27, carries the words "working-level number".
                 A block that quotes the ddE of a condition with a METASTABLE state carries the words "Probe 2" or the
                 word METASTABLE, because Probe 2's sentence applies to it.

THE SCOPED REPORT (added 11 November 2026).  With --scoped <report> --scoped-dir <its directory> the script validates
the second report, WRITEUP_SCOPED_<date>.md, and leaves the checks above alone.  Its eight checks:

  S1 numbers       as check 1, against the scoped report's own source table and scoped/derived.tsv
  S2 quotes        as check 2
  S3 citations     every key of the source table resolves to a tracked file, or to a table in the scoped directory
  S4 regeneration  scoped/make_scoped.py, re-run here, reproduces every .tsv of the scoped directory byte for byte
  S5 heads         the report names the two HEAD hashes it refers to, and nothing but its own files has changed since
  S6 levels        every block of the sections that hold claims (what is established, what must be remeasured, the
                   conclusions and hypotheses) carries a [level: ...] tag.  A tag names the level or levels of theory the
                   claim rests on; each must also stand in the sentence itself, and must be a level the sweep of committed
                   inputs finds in a directory the block's sources were measured in.  [level: comparison] is for claims
                   whose content is a comparison of methods and requires two levels in the sentence; [level: independent]
                   is for facts about the record.  So no claim drawn at one level can be stated without it.
  S7 d34           the report's statement about DIAGNOSTIC 34's two arms and the old level is recomputed here from the
                   raw outputs (d34_recompute, the same code as check 13 to 15) and must agree
  S8 conditional   the title and the first page say the report is conditional on D36 and that the level is not a decision

--selftest then proves each check fires: it perturbs a number, a quotation, a count, a citation, a withdrawn figure,
a self-label, a figure's .tsv, a HEAD hash and the not-tested list, and then, for the six added checks, the stated
commit, a DIAGNOSTIC 34 number, a raw control output, a recorded frozen ddE, an arm 1 leg, a claim's number, a rejected
claim and the post hoc label, and, for the four checks of 15 November 2026, the benchmark's derived value, the
denial beside the experimental number, the label of run A1, the post hoc label of a distance, the data-only label
of an orbital energy and a method listed as not tested under another spelling, and, for the four checks of that
evening, the derived value of the def2-SVPD benchmark, a detachment energy of a table, the word METASTABLE beside a
reached lowering, a number of DIAGNOSTIC 39 quoted under its citation and the post hoc label of a distance at
def2-SVPD, and, for the two checks of 16 November 2026, a detachment energy of DIAGNOSTIC 42 in its table, the count of
conditions bound at both states, the words NOT QUOTABLE beside one of its ddE, the word METASTABLE beside one that
Probe 2's sentence reaches and the post hoc label of a comparison, and, for the two checks of 18 November 2026, a
detachment energy of DIAGNOSTIC 43 in its table, the count of survivors, the count of outputs of DIAGNOSTIC 44, the
words working-level number beside a survivor's ddE and the mention of Probe 2 beside a ddE its sentence reaches, one at
a time, on a copy, and requires the right check to fail for the right reason.
"""
import argparse
import hashlib
import importlib.util
import os
import re
import shutil
import subprocess
import sys
import tarfile
import tempfile
from pathlib import Path

S6 = "phase2.2/cs_stage6"
MONTHS = "January|February|March|April|May|June|July|August|September|October|November|December"
XREF = (r"\b(?:sections?|Tables?|items?|Appendix|Figures?|Stages?|Phases?|Amendments?|Errata|Erratum|Probes?|Experiment|"
        r"Pass|pass|readings?|Readings?)\s+[A-Z]?\.?\d+[a-e]?(?:\.\d+)*(?:\s*(?:,|and|to)\s*\d+[a-e]?(?:\.\d+)*)*")

# ---- check 4: the withdrawal registry.  (pattern in the document, key of the withdrawing file, its sentence)
WITHDRAWN = [
    (r"0\.58", "D29P", "A figure of Sigma q^2 <= 0.58 e^2 was"),
    (r"3\.4 times", "D29P", "is about 3.4 times too high"),
    (r"0\.093(?!\d)", "D32E2", "+0.093 e IS IN NO OUTPUT AND NO RESULT."),      # the three-decimal figure, not a longer number that begins with it
    (r"factor of nine", "D32E2", "Both are WITHDRAWN."),
    (r"q_quad|0\.80 e\b", "E3", "q_quad = 0.80 e IS WITHDRAWN"),
    (r"5\.6165", "D27", "superseded for any future use by"),
    (r"110 conditions|16 controls", "D29E1", "The correct counts are 182 and 26."),
    (r"\b11 MB\b", "D38E1", "11 MB must not be used"),      # added 15 November 2026: DIAGNOSTIC 38's result, Erratum 1
]

# ---- check 5: the self-label registry.  (trigger in a section, label required there, key, the result's own sentence)
SELF_LABELS = [
    (r"\bD13\b", "PARTIAL", "D13", "PARTIAL: 3 of 6 sites failing at |q| = 1 pass at 0.25."),
    (r"(?<![-\w])D3\b(?!\.)", "INCONCLUSIVE", "D3", "VERDICT PER D3.3: INCONCLUSIVE"),
    (r"\bD24\b|[Ss]mear", "EXPLORATORY", "D24", "It is therefore EXPLORATORY and its status is a scoping"),
    (r"rescue|argon|ECP wall|Pauli wall|walled design", "NOT ADOPTED", "D33P", "THE WALL IS NOT ADOPTED AND IS NOT CALIBRATED."),
    (r"re-read", "post hoc", "D32E1", "THIS RECOMPUTATION IS POST HOC AND IS NOT A REPLACEMENT VERDICT."),
    (r"D34.{0,80}stage A|stage A needs no compute", "post hoc", "D34P", "STAGE A IS POST HOC AND IS LABELLED SO WHEREVER IT IS REPORTED."),
    # added 11 November 2026: the old-level numbers of S6-20 and S6-10, which two committed documents restrict
    (r"24\.826057|24\.8261|21\.816364|21\.8164", "STRAINED", "E3", "their barriers must not be quoted without that said in the same sentence"),
    (r"24\.826057|24\.8261", "PROVISIONAL", "D20", "must be quoted as provisional or not quoted"),
    # added 15 November 2026
    # changed on the evening of 15 November 2026: D37's label for it was NOT MEASURED; Amendment 1 has measured it
    (r"revDSD", "Amendment 1", "D37A1R", "revDSD-PBEP86-D4, which D37 could not run, run two ways."),
    (r"reaches the correlat|REACHES THE CORRELATION", "mechanism", "D38R", "A MECHANISM TEST. One geometry, the reactant only, no transition state, no ddE,"),
]

# ---- the revision of 11 November 2026 (DIAGNOSTIC 34): checks 11 to 16
EXPECTED_HEADS = {"chorismate-thesis-results": "0ede870", "chorismate-thesis-code": "ee28301"}   # 6c09d9c, 9f66ee3 on 16 November; 21d631b, 512eead on the evening of 15 November; 4af2afa, 1e47614 that morning; f461bd1, 23c2500 on 11 November
D34_RECORD = ("D34R", "D34T")                           # the committed result and its table: where a D34 number must stand
D34_TAGS = ("D34R", "D34T", "D34E1", "T3D", "T4W")      # a block tagged with any of these is a block about D34
D34_REF = "26.8248"                                     # the corrected reference barrier, kcal/mol
C3_TRIGGER = (r"largest positive charge|no positive charge|without a positive charge|onto the positive charge|"
              r"positive charges|sign[- ]swap|magnitude threshold")
# id, what the claim is, the [dv:] rows that hold the numbers behind it, the forms of it that the recomputation rejected
CLAIMS = [
    ("C1", "arm 1 already shows the relocation", ["DV113", "DV116", "DV117"], []),
    ("C2", "the old level shows none of it", ["DV115", "DV118"], []),
    ("C3", "the largest positive charge against the verdict (post hoc)",
     ["DV128", "DV129", "DV130", "DV140", "DV141", "DV142", "DV143", "DV144", "DV145", "DV146", "DV147", "DV149",
      "DV151", "DV152", "DV153", "DV154", "DV155", "DV156", "DV157"], []),
    ("C4", "DIAGNOSTIC 11's sign control, and DIAGNOSTIC 12's sign arm", ["DV158", "DV159", "DV160", "DV163"], []),
    ("C5", "S6-P and ADD-s1690_qp10000 are one electronic system", ["DV131", "DV132"], []),
    ("C6", "the size of the stability correction", ["DV119", "DV120", "DV121", "DV122", "DV123", "DV124"],
     [r"\b(?:8|eight) of (?:the )?(?:10|ten)\b[^.\n]{0,160}unstable", r"unstable[^.\n]{0,160}\b(?:8|eight) of (?:the )?(?:10|ten)\b",
      r"first (?:direct )?measurement[^.\n]{0,120}stability correction",
      r"stability correction[^.\n]{0,160}(?:never|not) (?:before )?been measured"]),
]
C4_NUMBERS = ("-0.0926", "+0.0002", "+3.9642", "+0.0195")   # q = -1 at dref and d450, then q = +1 at the same rungs

# ---- the revision of 15 November 2026 (DIAGNOSTICS 35 to 38): checks 17 to 20
H2K = 627.5094740631
# check 8: the ORCA keywords, in upper case, that request a name under another spelling
KW_ALIAS = {"M06-2X": ["M062X"], "PWPB95": ["RI-PWPB95"], "REVDSD-PBEP86": ["RI-REVDSD-PBEP86-D4/2021", "REVDSD-PBEP86-D4/2021"],
            "REVDSD-PBEP86-D4": ["RI-REVDSD-PBEP86-D4/2021", "REVDSD-PBEP86-D4/2021"],
            "D4": ["RI-REVDSD-PBEP86-D4/2021", "REVDSD-PBEP86-D4/2021"], "VV10": ["WB97M-V", "WB97X-V"]}
# check 17: the [dv:] rows that hold DIAGNOSTIC 36's barrier and its parts, and DIAGNOSTIC 38's two shifts for run B1
BENCH_DV = {"DV220": "E(CCSD(T))", "DV221": "E(0)", "DV222": "E(CCSD)", "DV223": "Triples Correction (T)"}
D38_DV = {"DV250": "E(0)", "DV251": "E(CORR)(corrected)"}
# check 19: numbers that belong to run A1 of DIAGNOSTIC 38 alone
A1_NUMBERS = (r"831\.648351321", r"2\.598042111", r"834\.351462380", r"0\.183558659", r"0\.001632748", r"0\.181542852",
              r"0\.016858947", r"\+0\.3539")
# check 20: DIAGNOSTIC 37's five distances from the benchmark and what is built on them; the orbital energies of T9
POSTHOC_LEVEL = (r"-2\.5196", r"-2\.9008", r"\+3\.5321", r"\+4\.5400", r"\+5\.8962", r"\b12\.45\b", r"\b53 per cent",
                 # added on the evening of 15 November 2026: revDSD-PBEP86-D4 at def2-SVPD against DIAGNOSTIC 41's benchmark,
                 # and def2-SVPD barriers against the def2-TZVPD benchmark; neither comparison is in a committed document
                 r"-1\.0554", r"-1\.6572", r"-1\.7032", r"-2\.0189", r"\+2\.9645")
T9_NUMBERS = (r"\+0\.7782", r"\+2\.4036", r"\+1\.3747", r"\+1\.6513", r"\+1\.6925", r"\+1\.3398", r"\+1\.1380")
EXPERIMENT = r"(?<![\d.])24\.5(?![\d])"
BENCHMARK = r"(?<![\d.])23\.86(?:03)?(?![\d])"

# ---- the revision of the evening of 15 November 2026 (DIAGNOSTICS 40 and 41, DIAGNOSTIC 37's Amendment 1): checks 21 to 24
EV = 27.211386245988                                    # eV per Eh, the constant DIAGNOSTIC 40 uses
# check 21: the [dv:] rows that hold DIAGNOSTIC 41's barrier and its parts
D41_DV = {"DV300": "E(CCSD(T))", "DV301": "E(0)", "DV302": "E(CCSD)", "DV303": "Triples Correction (T)"}
# check 23: the five lowerings DIAGNOSTIC 40's R40.7 says Probe 2's sentence reaches
REACHED = r"(?<![\d.])-(?:3\.3985|4\.1567|5\.3701|14\.3612|13\.0007)(?!\d)"
# check 24: what of DIAGNOSTIC 39 is held.  Its result and its amendment may be cited for the statement that they exist;
# its outputs and its table may not be sources at all
D39_DOCS = r"cs_stage6/STAGE6_DIAGNOSTIC39_(?:RESULT|AMENDMENT)"
D39_DATA = r"cs_stage6/(?:d39/|d39\.tsv$)"

# ---- the revision of 16 November 2026 (DIAGNOSTIC 42): checks 25 and 26
# check 26: the [dv:] rows that hold this write-up's own comparisons of DIAGNOSTIC 42's numbers; each is marked POST HOC there
D42_POSTHOC_DV = ("DV406", "DV407", "DV412", "DV413", "DV414", "DV415", "DV416", "DV417", "DV418", "DV419")
D42_NOTSTABLE = "ABOVE ZERO, MONOANION NOT STABLE"

NUM = re.compile(r"^[+-]?\d+(?:\.\d+)?(?:e[+-]?\d+)?$", re.I)
RANGE = re.compile(r"^(\d+(?:\.\d+)?)-(\d+(?:\.\d+)?)$")


def sha(p):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def git(repo, *args):
    return subprocess.run(["git", "-C", str(repo)] + list(args), capture_output=True, text=True)


def norm(s):
    return " ".join(s.replace("−", "-").split())


class Ctx:
    """Everything that does not depend on the document's text: the repositories, their tracked files, and the
    tables regenerated once from them."""

    def __init__(self, res, code, full):
        self.res, self.code, self.full = Path(res).resolve(), Path(code).resolve(), full
        self.tracked = {"R": set(git(self.res, "ls-files").stdout.splitlines()),
                        "C": set(git(self.code, "ls-files").stdout.splitlines())}
        self.cache = {}
        self.regen = None
        self.d34root = self.res / S6 / "d34"           # the negative tests point this at a perturbed copy
        self.d34cache = {}

    def text(self, repo, path, figdir):
        if repo == "W":
            p = Path(figdir) / path
        else:
            p = (self.res if repo == "R" else self.code) / path
        key = str(p)
        if repo == "W" or key not in self.cache:
            try:
                t = p.read_text(errors="replace").replace("−", "-")
            except OSError:
                t = None
            if repo == "W":
                return t
            self.cache[key] = t
        return self.cache[key]

    def regenerate(self, figdir):
        """Run the committed generators once, from the copies in the figures directory, into a scratch directory."""
        if self.regen is not None:
            return self.regen
        d = Path(tempfile.mkdtemp(prefix="cswu_regen_"))
        log = []
        for script, extra in (("make_figures.py", ["--tsv-only"]), ("sweep_pass1.py", [])) + ((("sweep_pass2.py", []),) if self.full else ()):
            r = subprocess.run([sys.executable, str(Path(figdir) / script), "--res", str(self.res), "--code", str(self.code),
                                "--out", str(d)] + extra, capture_output=True, text=True)
            log.append((script, r.returncode, (r.stdout + r.stderr).strip().splitlines()[-3:]))
        self.regen = (d, log)
        return self.regen


def parse_sources(doc):
    src = {}
    for m in re.finditer(r"^\|\s*([A-Z0-9]+)\s*\|\s*([RCW])\s*\|\s*(\S+)\s*\|\s*$", doc, re.M):
        src[m.group(1)] = (m.group(2), m.group(3))
    return src


def blocks_of(doc):
    """Yield (section heading, block text) with fenced blocks removed.  A table followed by a 'Sources for' block is
    joined to it, so the table's rows are checked against that tag."""
    out, sec, fence, cur = [], "(preamble)", False, []
    for line in doc.split("\n"):
        if line.strip().startswith("```"):
            fence = not fence
            continue
        if fence:
            continue
        if line.startswith("## "):
            if cur:
                out.append((sec, "\n".join(cur)))
                cur = []
            sec = line[3:].strip()
            out.append((sec, line))
            continue
        if not line.strip():
            if cur:
                out.append((sec, "\n".join(cur)))
                cur = []
        else:
            cur.append(line)
    if cur:
        out.append((sec, "\n".join(cur)))
    merged = []
    for sec, b in out:
        if b.lstrip().startswith("Sources for") and merged and merged[-1][1].lstrip().startswith("|"):
            merged[-1] = (merged[-1][0], merged[-1][1] + "\n" + b)
        else:
            merged.append((sec, b))
    return merged


def tags(block):
    keys, dvs = [], []
    block = re.sub(r"`[^`\n]+`", " ", block)          # a tag shown as an example in back-ticks is not a tag
    for m in re.finditer(r"\[src:\s*([^\]]+)\]", block):
        keys += [k.strip() for k in m.group(1).split(",") if k.strip()]
    for m in re.finditer(r"\[dv:\s*([^\]]+)\]", block):
        dvs += [k.strip() for k in m.group(1).split(",") if k.strip()]
    return keys, dvs


def occurs(num, text):
    """num as a whole number in text: not the middle or the truncation of a longer one."""
    return re.search(r"(?<![\d.])" + re.escape(num) + r"(?!\d)(?!\.\d)", text) is not None


def numerals_of(b):
    """The numerals of a block, as check 1 defines them.  Returns (items, exempt).  An item is (word, kind, candidates):
    'num' has one candidate, 'range' two, and 'ident' is an identifier that holds a decimal and must stand verbatim."""
    t = re.sub(r"<!--.*?-->", " ", b, flags=re.S)
    t = re.sub(r"\[(?:src|dv):[^\]]*\]", " ", t)
    t = re.sub(r"!\[[^\]]*\]\([^)]*\)", " ", t)
    ncode = len(re.findall(r"`[^`\n]+`", t))
    t = re.sub(r"`[^`\n]+`", " ", t)
    t, nd = re.subn(r"\b\d{1,2} (?:%s) \d{4}\b" % MONTHS, " ", t)
    t, nx = re.subn(XREF, " ", t)
    exempt = ncode + nd + nx
    lines = []
    for line in t.split("\n"):
        line, k = re.subn(r"^\s*\d+\.\s+", "", line)
        exempt += k
        lines.append(line)
    t = "\n".join(lines).replace("|", " ").replace("“", " ").replace("”", " ")
    items = []
    for w in t.split():
        w = w.strip("()[],;:\"'*").rstrip(".").strip("()[],;:\"'*")
        if w.startswith("+-"):
            w = w[2:]
        if not re.search(r"\d", w):
            continue
        if NUM.match(w):
            items.append((w, "num", [w.lstrip("+-")]))
        elif RANGE.match(w):
            items.append((w, "range", list(RANGE.match(w).groups())))
        elif re.search(r"[A-Za-z_]", w) and not re.search(r"\d\.\d", w):
            exempt += 1                                # an identifier: letters joined to whole numbers
        else:
            items.append((w, "ident", [w]))            # an identifier holding a decimal, or anything unclassified
    return items, exempt


def load(path, name):
    spec = importlib.util.spec_from_file_location(name, str(path))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def d34_recompute(ctx):
    """DIAGNOSTIC 34 again, by this file and not by make_figures.py: the 44 raw outputs under d34/, the 22 of
    frozen/sp_outputs.tar.xz, the 12 of standoff/ that carry DIAGNOSTIC 11's sign control and the 54 in-vacuo ones of
    DIAGNOSTIC 2's scan of both signs at site 1690 (qscan/outputs/qscan_outputs.tar.xz), through the committed
    readers (cs6_audit.read, cs6_d33b.stability and SHIFT, cs6_analyse.read_run).  Returns the quantities checks 13,
    14 and 15 compare; raises if a run cannot be read."""
    key = str(ctx.d34root)
    if key in ctx.d34cache:
        return ctx.d34cache[key]
    P = ctx.code / "phase2.2"
    audit, d33b, analyse = (load(P / (n + ".py"), "cswu_" + n) for n in ("cs6_audit", "cs6_d33b", "cs6_analyse"))
    H, root = analyse.H, ctx.res / S6
    frozen = {}
    for l in (root / "frozen" / "frozen.tsv").read_text().splitlines()[1:]:
        f = l.split("\t")
        if l.strip():
            frozen[f[0]] = (int(f[1]), float(f[2]))
    rec = {}
    for l in (root / "d34.tsv").read_text().splitlines()[1:]:
        f = l.split("\t")
        if len(f) >= 6:
            rec[(f[0], f[1])] = f
    names = list(frozen)
    if sorted(names) != sorted(n for (stg, n) in rec if stg == "A"):
        raise RuntimeError("d34.tsv and frozen/frozen.tsv do not name the same arrangements")

    def one(d, nsph):
        t = (Path(d) / "job.out").read_text(errors="replace")
        c = audit.read(t)
        if c is None or not c["norm"] or c["notconv"]:
            raise RuntimeError("%s is not a normally terminated, converged run" % d)
        e_tot, e_mm, _g = analyse.read_run(d, nsph)
        return dict(E=e_tot, Emm=e_mm, lw=c["maxlw"], dip=c["dip"], shifted=any(re.search(x, t) for x in d33b.SHIFT),
                    stab=d33b.stability(t))

    R = {}
    for arm in ("a1", "a2"):
        for name in ["S6-C"] + names:
            for st in ("R", "TS"):
                R[(arm, name, st)] = one(ctx.d34root / ("%s_%s_%s" % (arm, name, st)), 0 if name == "S6-C" else frozen[name][0])
    tmp = Path(tempfile.mkdtemp(prefix="cswu_frozen_"))
    try:
        with tarfile.open(root / "frozen" / "sp_outputs.tar.xz", "r:xz") as tf:
            tf.extractall(tmp, filter="data")
        for name in ["S6-C"] + names:
            for st in ("R", "TS"):
                d = tmp / name / ("%s_sp" % st)
                shutil.copy(root / "relax" / name / ("%s_sp" % st) / "combined.xyz", d / "combined.xyz")
                R[("A", name, st)] = one(d, 0 if name == "S6-C" else frozen[name][0])
    finally:
        shutil.rmtree(tmp, ignore_errors=True)

    def rise(arm, name, st):
        return R[(arm, name, st)]["lw"] - R[(arm, "S6-C", st)]["lw"]

    def dde(arm, name):
        return ((R[(arm, name, "TS")]["E"] - R[(arm, "S6-C", "TS")]["E"]) - (R[(arm, name, "R")]["E"] - R[(arm, "S6-C", "R")]["E"])) * H

    def over(arm, name):                               # at the four decimals the tables carry, as make_figures.py judges it
        return max(float("%.4f" % rise(arm, name, "R")), float("%.4f" % rise(arm, name, "TS"))) >= audit.D_LW

    def charges(name):
        ls = (root / "relax" / name / "R_sp" / "site.ORCAFF.prms").read_text().splitlines()
        k = next(i for i, l in enumerate(ls) if l.strip().startswith("$atoms"))
        cnt = int(ls[k + 1].split()[0])
        return [float(r.split()[2]) for r in ls[k + 2:k + 2 + cnt][-frozen[name][0]:]]

    D = dict(names=names, frozen=frozen, rec=rec, H=H)
    D["ctl"] = {arm: (R[(arm, "S6-C", "TS")]["E"] - R[(arm, "S6-C", "R")]["E"]) * H for arm in ("A", "a1", "a2")}
    D["ddeA"] = {n: dde("A", n) for n in names}
    D["a1"] = {n: dict(rise_R="%.4f" % rise("a1", n, "R"), rise_TS="%.4f" % rise("a1", n, "TS"),
                       dip_R="%.2f" % (R[("a1", n, "R")]["dip"] / R[("a1", "S6-C", "R")]["dip"]),
                       dip_TS="%.2f" % (R[("a1", n, "TS")]["dip"] / R[("a1", "S6-C", "TS")]["dip"]),
                       shifted="yes" if (R[("a1", n, "R")]["shifted"] or R[("a1", n, "TS")]["shifted"]) else "no",
                       dde="%.6f" % dde("a1", n)) for n in names}
    v = {}
    v["DV113"] = str(sum(1 for n in names if over("a1", n)))
    v["DV116"] = str(sum(1 for n in names if over("a1", n) != over("a2", n)))
    v["DV117"] = "%.4f" % max(rise("a1", n, st) for n in names for st in ("R", "TS"))
    v["DV134"] = str(sum(1 for n in names if D["a1"][n]["shifted"] == "yes"))
    v["DV115"] = str(sum(1 for n in names if over("A", n)))
    v["DV118"] = "%.4f" % max(rise("A", n, st) for n in names for st in ("R", "TS"))
    nopos = [n for n in names if max(charges(n)) <= 0.0]
    v["DV128"] = str(len(nopos))
    v["DV129"] = str(sum(1 for n in nopos if rec[("B2", n)][4] == "RELIABLE"))
    v["DV130"] = "%.4f" % max(rise("a2", n, st) for n in nopos for st in ("R", "TS"))
    qm = {n: [R[(arm, n, st)]["E"] - R[(arm, n, st)]["Emm"] for arm in ("A", "a1", "a2") for st in ("R", "TS")] for n in names}
    groups = []
    for n in names:
        for g in groups:
            if max(abs(a - b) for a, b in zip(qm[n], qm[g[0]])) < 1e-9:
                g.append(n)
                break
        else:
            groups.append([n])
    v["DV131"] = str(len(groups))
    v["DV132"] = "%.1e" % max(abs(a - b) for a, b in zip(qm["S6-P"], qm["ADD-s1690_qp10000"]))
    unst = [n for n in names if "UNSTABLE" in (R[("a2", n, "R")]["stab"]["first"], R[("a2", n, "TS")]["stab"]["first"])]
    v["DV119"] = str(len(unst))
    v["DV120"] = str(sum(1 for n in names for st in ("R", "TS") if R[("a2", n, st)]["stab"]["first"] == "UNSTABLE"))
    v["DV121"] = "/".join(sorted(set(str(R[("a2", n, st)]["stab"]["nneg"]) for n in unst for st in ("R", "TS"))))
    corr = max(abs(R[("a2", n, st)]["E"] - R[("a1", n, st)]["E"]) for n in names for st in ("R", "TS"))
    v["DV122"], v["DV123"] = "%.4f" % (corr * H), "%.9f" % corr
    v["DV124"] = str(sum(1 for n in names if float("%.6f" % dde("a1", n)) * float("%.6f" % dde("a2", n)) < 0))
    D["dv"] = v

    def e11(n):
        return audit.read((root / "standoff" / n / "job.out").read_text(errors="replace"))["E"]

    def d_shield(rung, q):                             # ddE bare minus ddE walled, both against the charge-free control
        x = {}
        for kind in ("bare", "wall"):
            a = "%s_%s_%s_ma" % (rung, kind, q)
            x[kind] = ((e11(a + "_TS") - e11("ctl_ma_TS")) - (e11(a + "_R") - e11("ctl_ma_R"))) * H
        return "%+.4f" % (x["bare"] - x["wall"])
    D["c4"] = (d_shield("dref", "m1"), d_shield("d450", "m1"), d_shield("dref", "p1"), d_shield("d450", "p1"))

    # DIAGNOSTIC 12's sign arm (C-A12.2), read from its committed table: d_shield at +1 and at -1 on ten sites
    t12 = (root / "sitescan" / "DIAGNOSTIC12_RESULT.txt").read_text()
    arm = t12[t12.index("C-A12.2"):t12.index("BASIS ARM")]
    sg = [(abs(float(a)), abs(float(b))) for a, b in re.findall(r"^\s+\d+\s+\d\.\d{3}\s+([+-]\d+\.\d{4})\s+([+-]\d+\.\d{4})\s*$", arm, re.M)]
    v.update(DV158="%.4f" % max(b for _a, b in sg), DV159="%.4f" % max(a for a, _b in sg), DV160=str(len(sg)),
             DV163=str(sum(1 for a, b in sg if a > b)))

    # DIAGNOSTIC 2: both signs at site 1690.  The Loewdin and dipole legs only; these runs have no stability analysis.
    tmp = Path(tempfile.mkdtemp(prefix="cswu_d2_"))
    try:
        with tarfile.open(root / "qscan" / "outputs" / "qscan_outputs.tar.xz", "r:xz") as tf:
            tf.extractall(tmp, filter="data")

        def leg(cond, run):
            o = []
            for st in ("R", "TS"):
                c = audit.read((tmp / cond / ("%s_%s" % (run, st)) / "job.out").read_text(errors="replace"))
                b = audit.read((tmp / cond / ("bare_%s" % st) / "job.out").read_text(errors="replace"))
                o.append((c["maxlw"] - b["maxlw"], c["dip"] / b["dip"]))
            return o
        fired = {"p": 0, "m": 0}
        nneg = 0
        for cond in ("C1", "C3", "C5"):
            for mag in ("0050", "0250", "0500", "1000"):
                for sign in ("p", "m"):
                    lg = leg(cond, "q%s%s" % (sign, mag))
                    fired[sign] += sum(1 for r_, d_ in lg if float("%.4f" % r_) >= audit.D_LW or float("%.2f" % d_) >= audit.R_DIP)
                    nneg += 2 * (sign == "m")
        pos, neg = leg("C5", "qp1000"), leg("C5", "qm1000")
        v.update(DV151="%.4f" % max(r_ for r_, _d in pos), DV152="%.4f" % max(r_ for r_, _d in neg),
                 DV153="%.2f" % max(d_ for _r, d_ in pos), DV154="%.2f" % max(d_ for _r, d_ in neg),
                 DV155=str(nneg), DV156=str(fired["m"]), DV157=str(fired["p"]))
    finally:
        shutil.rmtree(tmp, ignore_errors=True)
    ctx.d34cache[key] = D
    return D


# ===================================================================== the scoped report (11 November 2026)
LEVEL_RE = re.compile(r"(?:DLPNO-CCSD\(T\)|CAM-B3LYP-D3BJ|B3LYP-D3BJ|wB97X-D3|CAM-B3LYP|B3LYP|MP2|HF)/(?:ma-def2-SVP|def2-TZVPD|def2-SVPD|def2-SVP|6-31G\(d\))")
SCOPED_LEVEL_SECTIONS = ("3.", "4.", "7.")
SCOPED_ORDER = ["numbers", "quotes", "citations", "regeneration", "heads", "levels", "d34", "conditional"]


def scoped_sources(doc):
    """The scoped report's source table: | key | repo | path | measured in |.  Returns key -> (repo, path, [directories])."""
    src = {}
    for m in re.finditer(r"^\|\s*([A-Z0-9]+)\s*\|\s*([RCW])\s*\|\s*(\S+)\s*\|\s*([^|\n]*?)\s*\|\s*$", doc, re.M):
        src[m.group(1)] = (m.group(2), m.group(3), [d.strip() for d in m.group(4).split(",") if d.strip() not in ("", "-")])
    return src


def run_scoped(doc, sdir, ctx):
    F = {k: [] for k in SCOPED_ORDER}
    st = dict(sourced=0, derived=0, quotes=0, keys=0, tsv=0, blocks=0, tagged=0, d34=0)
    sdir = Path(sdir)
    src = scoped_sources(doc)
    st["keys"] = len(src)
    dv = {}
    if (sdir / "derived.tsv").is_file():
        for l in (sdir / "derived.tsv").read_text().splitlines():
            f = l.split("\t")
            if len(f) > 1 and not l.startswith("#") and f[0] != "id":
                dv[f[0]] = f
    # S4: regenerate, and read the sweep for S6
    tmp = Path(tempfile.mkdtemp(prefix="cswu_scoped_"))
    levels_in = {}
    try:
        r = subprocess.run([sys.executable, str(sdir / "make_scoped.py"), "--res", str(ctx.res), "--out", str(tmp)], capture_output=True, text=True)
        if r.returncode != 0:
            F["regeneration"].append("make_scoped.py exited %d: %s" % (r.returncode, (r.stderr or r.stdout).strip().replace("\n", " | ")[-200:]))
        for q in sorted(sdir.glob("*.tsv")):
            st["tsv"] += 1
            if not (tmp / q.name).is_file():
                F["regeneration"].append("%s is not produced when make_scoped.py is re-run here" % q.name)
            elif (tmp / q.name).read_bytes() != q.read_bytes():
                F["regeneration"].append("%s differs from the file regenerated from the committed sources" % q.name)
        if (tmp / "S1_inputs_by_directory.tsv").is_file():
            for l in (tmp / "S1_inputs_by_directory.tsv").read_text().splitlines():
                f = l.split("\t")
                if len(f) == 5 and not l.startswith("#") and f[0] != "directory":
                    levels_in[f[0]] = set(re.sub(r" \d+$", "", x) for x in f[4].split("; "))
    finally:
        shutil.rmtree(tmp, ignore_errors=True)
    # S3
    for k, (repo, path, dirs) in sorted(src.items()):
        if repo == "W":
            if not (sdir / path).is_file():
                F["citations"].append("%s: %s is not in the scoped directory" % (k, path))
        elif path not in ctx.tracked[repo]:
            F["citations"].append("%s: %s is not tracked by git in the %s repository" % (k, path, "results" if repo == "R" else "code"))
        for d in dirs:
            if levels_in and d not in levels_in:
                F["citations"].append("%s: 'measured in' names %s, which holds no committed input" % (k, d))
    # S1, S2, S6, block by block
    for sec, b in blocks_of(doc):
        if b.startswith("#") or "| key | repo | path |" in b or "HEAD this document refers to" in b:
            continue
        keys, dvs = tags(b)
        lt = re.search(r"\[level:([^\]]*)\]", b)
        prose = re.sub(r"\[(?:src|dv|level):[^\]]*\]", " ", b)
        st["blocks"] += 1
        for k in keys:
            if k not in src:
                F["citations"].append("section '%s': [src: %s] is not in the source table" % (sec, k))
        for d in dvs:
            if d not in dv:
                F["numbers"].append("section '%s': [dv: %s] is not a row of scoped/derived.tsv" % (sec, d))
        texts = [t for t in (ctx.text(src[k][0], src[k][1], sdir) for k in keys if k in src) if t]
        dvals = set(dv[d][1].lstrip("+-") for d in dvs if d in dv)
        for w, kind, cand in numerals_of(re.sub(r"\[level:[^\]]*\]", " ", b))[0]:
            if kind in ("range", "ident") and any(w in x for x in texts):
                st["sourced"] += 1
                continue
            if kind == "ident":
                F["numbers"].append("section '%s': '%s' is not verbatim in any [src:] file of its block (%s)" % (sec, w, ", ".join(keys) or "no src tag"))
                continue
            for c in cand:
                if c in dvals:
                    st["derived"] += 1
                elif any(occurs(c, x) for x in texts):
                    st["sourced"] += 1
                else:
                    F["numbers"].append("section '%s': %s is in no [src:] file of its block (%s) and is no [dv:] value (%s): ...%s..."
                                        % (sec, c, ", ".join(keys) or "no src tag", ", ".join(dvs) or "no dv tag", norm(b)[max(0, norm(b).find(c) - 50):norm(b).find(c) + 30]))
        for q in re.findall(r"“([^”]+)”", b):
            st["quotes"] += 1
            if not any(norm(q) in norm(x) for x in texts):
                F["quotes"].append("section '%s': not verbatim in any [src:] file of its block (%s): “%s”" % (sec, ", ".join(keys) or "no src tag", q[:70]))
        if sec.startswith(SCOPED_LEVEL_SECTIONS):
            if not lt:
                F["levels"].append("section '%s': a block carries no [level:] tag: %s..." % (sec, norm(b)[:80]))
                continue
            st["tagged"] += 1
            named = [x.strip() for x in lt.group(1).split(";") if x.strip()]
            inprose = set(LEVEL_RE.findall(prose))
            if named == ["independent"]:
                continue
            if named == ["comparison"]:
                if len(inprose) < 2:
                    F["levels"].append("section '%s': tagged a comparison but names fewer than two levels in the sentence: %s..." % (sec, norm(b)[:70]))
                continue
            cited = [d for k in keys if k in src for d in src[k][2]]
            union = set().union(*[levels_in.get(d, set()) for d in cited]) if cited else set()
            for lv in named:
                if not LEVEL_RE.fullmatch(lv):
                    F["levels"].append("section '%s': [level: %s] is not a functional/basis the script knows" % (sec, lv))
                elif lv not in inprose:
                    F["levels"].append("section '%s': %s is named in the tag but not in the sentence: %s..." % (sec, lv, norm(b)[:70]))
                elif cited and lv not in union:
                    F["levels"].append("section '%s': %s is not a level of any directory its sources were measured in (%s): %s..."
                                       % (sec, lv, ", ".join(sorted(set(cited))), norm(b)[:60]))
            for k in keys:
                if k in src and src[k][2]:
                    have = set().union(*[levels_in.get(d, set()) for d in src[k][2]])
                    if have and not (have & set(named)):
                        F["levels"].append("section '%s': cites %s, measured in %s, and names none of its levels (%s)"
                                           % (sec, k, ", ".join(src[k][2]), "; ".join(sorted(have))[:120]))
    # S5
    own = {"chorismate-thesis-results": re.compile(r"^phase2\.2/cs_stage6/(WRITEUP_SCOPED_\d{8}\.md|scoped/.*)$"),
           "chorismate-thesis-code": re.compile(r"^phase2\.2/cs_writeup_check\.py$")}
    for name, repo in (("chorismate-thesis-results", ctx.res), ("chorismate-thesis-code", ctx.code)):
        m = re.search(r"\|\s*%s\s*\|[^|]*\|\s*`([0-9a-f]{40})`" % re.escape(name), doc)
        head = git(repo, "rev-parse", "HEAD").stdout.strip()
        if not m:
            F["heads"].append("the report records no HEAD for %s" % name)
            continue
        if m.group(1) != head:
            anc = git(repo, "merge-base", "--is-ancestor", m.group(1), head).returncode == 0
            changed = git(repo, "diff", "--name-only", m.group(1), head).stdout.split() if anc else None
            if not anc or any(not own[name].match(c) for c in changed):
                F["heads"].append("%s: the report records %s and HEAD is %s" % (name, m.group(1), head))
        for l in git(repo, "status", "--porcelain", "--untracked-files=all").stdout.splitlines():
            if not own[name].match(l[3:].strip()):
                F["heads"].append("%s: the working tree holds a change that is not the report's own: %s" % (name, l.strip()))
    # S7
    try:
        D = d34_recompute(ctx)
        for sv, d in (("SV60", "DV113"), ("SV61", "DV115"), ("SV62", "DV117"), ("SV63", "DV118"), ("SV67", "DV116")):
            st["d34"] += 1
            if sv not in dv or dv[sv][1] != D["dv"][d]:
                F["d34"].append("%s is %s in scoped/derived.tsv; recomputed here from the raw outputs: %s" % (sv, dv[sv][1] if sv in dv else "absent", D["dv"][d]))
    except (SystemExit, Exception) as e:
        F["d34"].append("DIAGNOSTIC 34 could not be recomputed: %s" % str(e)[:200])
    # S8
    title = doc.split("\n", 1)[0]
    first = doc.split("\n## ", 1)[0]
    if "conditional on D36" not in title:
        F["conditional"].append("the title does not say the report is conditional on D36")
    for need in ("conditional on D36", "not a decision", "no amendment"):
        if need not in first:
            F["conditional"].append("the first page does not say '%s'" % need)
    return F, st


def report_scoped(F, st, out=sys.stdout):
    note = {"numbers": "%d sourced in a named committed file, %d derived (scoped/derived.tsv)" % (st["sourced"], st["derived"]),
            "quotes": "%d quotations" % st["quotes"], "citations": "%d keys in the source table" % st["keys"],
            "regeneration": "%d tables regenerated" % st["tsv"], "heads": "two repositories",
            "levels": "%d blocks carry a level tag" % st["tagged"], "d34": "%d values recomputed from raw outputs" % st["d34"],
            "conditional": "title and first page"}
    bad = 0
    for i, k in enumerate(SCOPED_ORDER, 1):
        print("S%d %-12s %s   (%s)" % (i, k, "PASS" if not F[k] else "FAIL %2d" % len(F[k]), note[k]), file=out)
        for m in F[k][:40]:
            print("      - " + m, file=out)
        bad += len(F[k])
    print("RESULT: %s" % ("ALL EIGHT SCOPED CHECKS PASS" if not bad else "%d FAILURE(S)" % bad), file=out)
    return bad


def selftest_scoped(doc, sdir, ctx):
    def sub(a, b):
        def f(d, sd):
            if a not in d:
                raise RuntimeError("selftest target not found: %s" % a[:60])
            return d.replace(a, b, 1), sd
        return f

    def tsv(name, a, b):
        def f(d, sd):
            q = Path(sd) / name
            t = q.read_text()
            if a not in t:
                raise RuntimeError("selftest target not found in %s: %s" % (name, a[:60]))
            q.write_text(t.replace(a, b, 1))
            return d, sd
        return f
    def swap_rsafe(d, sd):
        a1 = "It was measured by D11 along one axis from one site at B3LYP-D3BJ/ma-def2-SVP"
        a2 = "[dv: SV10] [level: B3LYP-D3BJ/ma-def2-SVP]"
        if a1 not in d or a2 not in d:
            raise RuntimeError("selftest target not found: r_safe")
        return d.replace(a1, a1.replace("B3LYP-D3BJ/ma-def2-SVP", "CAM-B3LYP-D3BJ/def2-SVPD"), 1).replace(a2, "[dv: SV10] [level: CAM-B3LYP-D3BJ/def2-SVPD]", 1), sd
    tests = [
        ("a number", "numbers", "26.8348", sub("**The reference barrier is 26.8248 kcal/mol", "**The reference barrier is 26.8348 kcal/mol")),
        ("a quotation", "quotes", "not verbatim", sub("“the production functional is NOT changed by this document”", "“the production functional is not changed by this document”")),
        ("a regenerated table", "regeneration", "S1_inputs_by_directory.tsv", tsv("S1_inputs_by_directory.tsv", "d34\t44\t44\t", "d34\t44\t43\t")),
        ("a HEAD hash", "heads", "records", sub("`133e212c", "`233e212c")),
        ("a conclusion stripped of its level tag", "levels", "carries no [level:] tag",
         sub("it is the number every lowering at this level is quoted against. [src: D33A] [level: CAM-B3LYP-D3BJ/def2-SVPD]",
             "it is the number every lowering at this level is quoted against. [src: D33A]")),
        ("a level left in the tag but taken out of the sentence", "levels", "named in the tag but not in the sentence",
         sub("**The reference barrier is 26.8248 kcal/mol at CAM-B3LYP-D3BJ/def2-SVPD.**", "**The reference barrier is 26.8248 kcal/mol.**")),
        ("r_safe restated at a level it was not measured at", "levels", "is not a level of any directory its sources were measured in", swap_rsafe),
        ("a comparison that names one level", "levels", "fewer than two levels",
         sub("nearer B3LYP-D3BJ/def2-TZVPD than CAM-B3LYP-D3BJ/def2-TZVPD, D28's direction", "nearer B3LYP than CAM-B3LYP, D28's direction")),
        ("the arm 1 count", "d34", "SV60", tsv("derived.tsv", "SV60\t7\t", "SV60\t6\t")),
        ("the condition in the title", "conditional", "title", sub("a scoped report, conditional on D36", "a scoped report")),
    ]
    lines, fired = ["NEGATIVE TESTS (scoped report)"], 0
    for name, check, reason, fn in tests:
        tmp = Path(tempfile.mkdtemp(prefix="cswu_sneg_"))
        try:
            sd = tmp / "scoped"
            shutil.copytree(sdir, sd, ignore=shutil.ignore_patterns("__pycache__"))
            d2, _ = fn(doc, sd)
            if d2 == doc and not any((sd / q.name).read_bytes() != q.read_bytes() for q in Path(sdir).glob("*.tsv")):
                raise RuntimeError("selftest '%s' changed nothing" % name)
            F, _s = run_scoped(d2, sd, ctx)
            hit = [m for m in F[check] if reason in m]
            fired += bool(hit)
            lines.append("   %-52s -> check '%s' %s%s" % (name, check, "FIRED" if hit else "DID NOT FIRE", (": " + hit[0][:100]) if hit else ""))
        finally:
            shutil.rmtree(tmp, ignore_errors=True)
    lines.append("   %d of %d negative tests fired" % (fired, len(tests)))
    print("\n".join(lines))
    return fired == len(tests), lines


def run_checks(doc, figdir, ctx, docpath=None, mode="pre"):
    F = {k: [] for k in ("numbers", "quotes", "citations", "withdrawn", "selflabels", "controls", "counts", "not_tested",
                         "regeneration", "heads", "exp_heads", "d34_numbers", "d34_controls", "arm1_legs", "claims", "posthoc",
                         "benchmark", "experiment", "d38_a1", "posthoc_lvl",
                         "benchmark_sv", "detachment", "metastable", "d39_held", "bound_d42", "not_quotable", "survivors_d43", "working_level")}
    stats = dict(sourced=0, derived=0, estimated=0, exempt=0, quotes=0, files=0, controls=0, worst="", counts=0, names=0, tsv=0,
                 d34_blocks=0, d34_rec=0, d34_dv=0, d34_ctl="", legs=0, claim_dvs=0, c3_blocks=0,
                 bench="", exp_blocks=0, a1_blocks=0, ph_blocks=0, t9_blocks=0,
                 bench_sv="", vde=0, vde_neg=0, reached_blocks=0, d39_blocks=0,
                 vde42=0, vde42_neg=0, both42=0, cond42=0, nq_blocks=0, ph42_blocks=0,
                 vde43=0, vde43_neg=0, cond43=0, rel43=0, both43=0, surv43=0, out44=0, wl_blocks=0, p2_blocks=0)
    doc = doc.replace("−", "-")
    src = parse_sources(doc)

    # ---------------------------------------------------------------- 3 citations
    for k, (repo, path) in sorted(src.items()):
        stats["files"] += 1
        if repo == "W":
            if not (Path(figdir) / path).is_file():
                F["citations"].append("%s: figures/%s does not exist" % (k, path))
        elif path not in ctx.tracked[repo]:
            F["citations"].append("%s: %s is not tracked by git in the %s repository" % (k, path, "results" if repo == "R" else "code"))
    used = set()
    for _, b in blocks_of(doc):
        ks, _d = tags(b)
        used |= set(ks)
    for k in sorted(used - set(src)):
        F["citations"].append("tag uses the key %s, which Appendix C does not define" % k)
    for m in re.finditer(r"`([^`\n]+)`", doc):
        s = m.group(1)
        if "/" in s and re.search(r"\.(py|tsv|txt|md|png|svg)$", s):
            if s == "phase2.2/cs_writeup_check.py" or s == "figures/PROVENANCE.txt":
                continue                               # this script, and the record it writes
            if s.startswith("figures/"):
                if not (Path(figdir) / s[len("figures/"):]).is_file():
                    F["citations"].append("back-ticked path %s does not exist in the figures directory" % s)
            elif s not in ctx.tracked["R"] and s not in ctx.tracked["C"]:
                F["citations"].append("back-ticked path %s is tracked in neither repository" % s)
        elif s == "figures/":
            continue
        elif s.endswith("/") or (s.startswith("phase2.2/") and "." not in s.split("/")[-1]):
            if not any(t.startswith(s.rstrip("/") + "/") for t in ctx.tracked["R"] | ctx.tracked["C"]):
                F["citations"].append("back-ticked directory %s holds no tracked file" % s)

    # ---------------------------------------------------------------- derived values
    dv = {}
    dvt = Path(figdir) / "derived.tsv"
    if dvt.is_file():
        for l in dvt.read_text().splitlines():
            if l.startswith("#") or l.startswith("id\t") or not l.strip():
                continue
            f = l.split("\t")
            dv[f[0]] = f
    else:
        F["regeneration"].append("figures/derived.tsv is missing")

    # ---------------------------------------------------------------- 1 numbers and 2 quotes
    for sec, b in blocks_of(doc):
        if b.startswith("#") or "| key | repo | path |" in b or "HEAD this document refers to" in b:
            continue
        keys, dvs = tags(b)
        texts = [ctx.text(src[k][0], src[k][1], figdir) or "" for k in keys if k in src]
        ntexts = None
        for d in dvs:
            if d not in dv:
                F["numbers"].append("section '%s': tag names %s, which figures/derived.tsv does not hold" % (sec, d))
        dvals = set(dv[d][1].lstrip("+-") for d in dvs if d in dv)
        est = "(estimated" in b
        # quotations first, on the untouched block
        for q in re.findall(r"“([^”]+)”", b):
            stats["quotes"] += 1
            if ntexts is None:
                ntexts = [norm(t) for t in texts]
            if not any(norm(q) in t for t in ntexts):
                F["quotes"].append("section '%s': not verbatim in any [src:] file of its block (%s): “%s”"
                                   % (sec, ", ".join(keys) or "no src tag", q[:90]))
        items, ex = numerals_of(b)
        stats["exempt"] += ex
        for w, kind, cand in items:
            if kind == "range" and any(w in x for x in texts):
                stats["sourced"] += 1
                continue
            if kind == "ident":
                if any(w in x for x in texts):
                    stats["sourced"] += 1
                else:
                    F["numbers"].append("section '%s': '%s' is not verbatim in any [src:] file of its block (%s)"
                                        % (sec, w, ", ".join(keys) or "no src tag"))
                continue
            for c in cand:
                if c in dvals:
                    stats["derived"] += 1
                elif any(occurs(c, x) for x in texts):
                    stats["sourced"] += 1
                elif est:
                    stats["estimated"] += 1
                else:
                    F["numbers"].append("section '%s': %s is in no [src:] file of its block (%s) and is no [dv:] value (%s): ...%s..."
                                        % (sec, c, ", ".join(keys) or "no src tag", ", ".join(dvs) or "no dv tag",
                                           norm(b)[max(0, norm(b).find(c) - 50):norm(b).find(c) + 30]))

    # ---------------------------------------------------------------- 4 withdrawn
    secs = {}
    for sec, b in blocks_of(doc):
        secs.setdefault(sec, []).append(b)
    for pat, key, sentence in WITHDRAWN:
        if key not in src:
            F["withdrawn"].append("registry key %s is not in Appendix C" % key)
            continue
        if norm(sentence) not in norm(ctx.text(src[key][0], src[key][1], figdir) or ""):
            F["withdrawn"].append("registry sentence is not in %s: %s" % (src[key][1], sentence))
        for sec, bs in secs.items():
            if re.search("withdrawn", sec, re.I):
                continue
            body = "\n".join(x for x in bs if "| key | repo | path |" not in x)
            m = re.search(pat, body)
            if m:
                F["withdrawn"].append("the withdrawn figure '%s' appears in section '%s', outside the withdrawal section" % (m.group(0), sec))

    # ---------------------------------------------------------------- 5 self-labels
    for pat, label, key, sentence in SELF_LABELS:
        if key not in src:
            F["selflabels"].append("registry key %s is not in Appendix C" % key)
            continue
        if norm(sentence) not in norm(ctx.text(src[key][0], src[key][1], figdir) or ""):
            F["selflabels"].append("registry sentence is not in %s: %s" % (src[key][1], sentence))
        for sec, bs in secs.items():
            if not re.match(r"\d+\.", sec):
                continue
            body = "\n".join(bs)
            if re.search(pat, body) and not re.search(re.escape(label), body, re.I):
                F["selflabels"].append("section '%s' uses %s (matched /%s/) and does not carry its self-label %s" % (sec, key, pat, label))

    # ---------------------------------------------------------------- 9 regeneration, 6 controls, 7 counts, 8 not tested
    rd, log = ctx.regenerate(figdir)
    for script, rc, tail in log:
        if rc != 0:
            F["regeneration"].append("%s exited %d: %s" % (script, rc, " | ".join(tail)))
    for p in sorted(rd.glob("*.tsv")):
        stats["tsv"] += 1
        q = Path(figdir) / p.name
        if not q.is_file():
            F["regeneration"].append("%s is regenerated but not present in the figures directory" % p.name)
        elif q.read_bytes() != p.read_bytes():
            F["regeneration"].append("%s differs from the file regenerated from committed sources" % p.name)
    for stem in sorted(set(x.name.split(".")[0] for x in Path(figdir).glob("F[0-9]*_*.png")) |
                       set(x.name[:-4] for x in rd.glob("F[0-9]*_*.tsv") if not x.name.startswith("F1b_"))):
        for ext in (".svg", ".png", ".tsv"):
            if not (Path(figdir) / (stem + ext)).is_file():
                F["regeneration"].append("figure %s has no %s" % (stem, ext))
    if not ctx.full and not (Path(figdir) / "T2_names.tsv").is_file():
        F["regeneration"].append("T2_names.tsv is missing")

    ct = rd / "controls.tsv"
    need = {"C04": "D33 stage A corrected reference", "C06": "production reference barrier"}
    seen, nfro, worst = set(), 0, (0.0, "")
    if ct.is_file():
        for l in ct.read_text().splitlines():
            if l.startswith("#") or l.startswith("id\t"):
                continue
            f = l.split("\t")
            stats["controls"] += 1
            seen.add(f[0])
            nfro += f[0].startswith("C07.")
            rel = float(f[5]) / float(f[6]) if float(f[6]) else (0.0 if float(f[5]) == 0 else 9e9)
            if rel >= worst[0]:
                worst = (rel, "%s %s: |diff| %s against a tolerance of %s" % (f[0], f[1], f[5], f[6]))
            if f[7] != "PASS":
                F["controls"].append("%s %s: recomputed %s, record %s (%s)" % (f[0], f[1], f[2], f[3], f[4]))
        stats["worst"] = worst[1]
        for k, what in need.items():
            if k not in seen:
                F["controls"].append("the %s control (%s) was not computed" % (what, k))
        nrec = sum(1 for l in ((ctx.res / S6 / "frozen" / "frozen.tsv").read_text().splitlines()[1:]) if l.strip())
        if nfro != nrec:
            F["controls"].append("frozen.tsv has %d ddE values and %d were recomputed" % (nrec, nfro))
    else:
        F["controls"].append("controls.tsv was not regenerated")

    for did, f in sorted(dv.items()):
        if f[2] != "count":
            continue
        m = re.match(r"number of committed (job\.out|job\.inp) files under (\S+)$", f[3])
        m2 = re.match(r"number of committed job\.out\.xz under phase2\.2/cs_stage4e/", f[3])
        if m:
            n = sum(1 for t in ctx.tracked["R"] if t.startswith(m.group(2) + "/") and t.endswith("/" + m.group(1)))
        elif m2:
            n = sum(1 for t in ctx.tracked["R"] if re.match(r"phase2\.2/cs_stage4e/L[1-5]_[^/]+/job\.out\.xz$", t))
        else:
            continue
        stats["counts"] += 1
        if str(n) != f[1]:
            F["counts"].append("%s: derived.tsv says %s, git ls-files gives %d (%s)" % (did, f[1], n, f[3]))

    m = re.search(r"<!-- not-tested-begin -->(.*?)<!-- not-tested-end -->", doc, re.S)
    tok = {}
    tb = rd / "T1b_keyword_tokens.tsv"
    if tb.is_file():
        for l in tb.read_text().splitlines():
            if l.startswith("#") or l.startswith("token\t"):
                continue
            a, i, o = l.split("\t")
            tok[a] = int(i) + int(o)
    if not m:
        F["not_tested"].append("the document has no not-tested block")
    elif not tok:
        F["not_tested"].append("the keyword-line census was not regenerated")
    else:
        for name in re.findall(r"`([^`]+)`", m.group(1)):
            stats["names"] += 1
            n_kw = tok.get(name.upper(), 0)
            via = [a for a in KW_ALIAS.get(name.upper(), []) if tok.get(a, 0)]
            if n_kw or via:
                F["not_tested"].append("%s is listed as not tested and stands on %d keyword line(s)%s"
                                       % (name, n_kw + sum(tok[a] for a in via), (" as " + ", ".join(via)) if via else ""))

    # ---------------------------------------------------------------- 10 heads
    own_res = re.compile(r"^phase2\.2/cs_stage6/(WRITEUP_CONSOLIDATED_\d{8}\.md|figures/.*)$")
    own_code = re.compile(r"^phase2\.2/cs_writeup_check\.py$")
    for name, repo, own in (("chorismate-thesis-results", ctx.res, own_res), ("chorismate-thesis-code", ctx.code, own_code)):
        m = re.search(r"\|\s*%s\s*\|[^|]*\|\s*`([0-9a-f]{40})`" % re.escape(name), doc)
        head = git(repo, "rev-parse", "HEAD").stdout.strip()
        if not m:
            F["heads"].append("the document records no HEAD for %s" % name)
        elif m.group(1) != head:
            anc = git(repo, "merge-base", "--is-ancestor", m.group(1), head).returncode == 0
            changed = git(repo, "diff", "--name-only", m.group(1), head).stdout.split() if anc else None
            if not anc or any(not own.match(c) for c in changed):
                F["heads"].append("%s: the document records %s and HEAD is %s" % (name, m.group(1), head))
        for l in git(repo, "status", "--porcelain", "--untracked-files=all").stdout.splitlines():
            p = l[3:].strip().strip('"')
            if not own.match(p):
                F["heads"].append("%s: the working tree holds a change that is not the write-up's own: %s" % (name, l.strip()))

    # ================================================================ the six checks of 11 November 2026
    # ---------------------------------------------------------------- 11 the commits the revision was built from
    for name, repo, own in (("chorismate-thesis-results", ctx.res, own_res), ("chorismate-thesis-code", ctx.code, own_code)):
        pre = EXPECTED_HEADS[name]
        m = re.search(r"\|\s*%s\s*\|[^|]*\|\s*`([0-9a-f]{40})`" % re.escape(name), doc)
        if not m or not m.group(1).startswith(pre):
            F["exp_heads"].append("the document does not state %s as the HEAD of %s (it states %s)" % (pre, name, m.group(1)[:12] if m else "none"))
        full = git(repo, "rev-parse", "--verify", "--quiet", pre + "^{commit}").stdout.strip()
        head = git(repo, "rev-parse", "HEAD").stdout.strip()
        if not full:
            F["exp_heads"].append("%s holds no commit %s" % (name, pre))
        elif full != head:
            anc = git(repo, "merge-base", "--is-ancestor", full, head).returncode == 0
            changed = git(repo, "diff", "--name-only", full, head).stdout.split() if anc else None
            if not anc or any(not own.match(c) for c in changed):
                F["exp_heads"].append("%s: HEAD is %s, which is not %s nor that commit with only the write-up's own files added" % (name, head[:12], pre))

    # ---------------------------------------------------------------- 12 every DIAGNOSTIC 34 number
    recs = [ctx.text(src[k][0], src[k][1], figdir) or "" for k in D34_RECORD if k in src]
    for k in D34_RECORD:
        if k not in src:
            F["d34_numbers"].append("Appendix C does not define %s" % k)
    regen_dv = {}
    if (rd / "derived.tsv").is_file():
        for l in (rd / "derived.tsv").read_text().splitlines():
            f = l.split("\t")
            if len(f) > 1 and not l.startswith("#"):
                regen_dv[f[0]] = f[1]
    for sec, b in blocks_of(doc):
        if b.startswith("#") or "| key | repo | path |" in b or "HEAD this document refers to" in b:
            continue
        keys, dvs = tags(b)
        if not set(keys) & set(D34_TAGS):
            continue
        stats["d34_blocks"] += 1
        dvals = set(dv[d][1].lstrip("+-") for d in dvs if d in dv)
        for d in dvs:
            if d in dv and regen_dv.get(d) != dv[d][1]:
                F["d34_numbers"].append("section '%s': %s is %s in figures/derived.tsv and %s when make_figures.py is re-run here"
                                        % (sec, d, dv[d][1], regen_dv.get(d, "absent")))
        for w, kind, cand in numerals_of(b)[0]:
            if kind in ("range", "ident"):
                if any(w in x for x in recs):
                    stats["d34_rec"] += 1
                    continue
                if kind == "ident":
                    F["d34_numbers"].append("section '%s': '%s' is in neither the D34 result nor d34.tsv" % (sec, w))
                    continue
            for c in cand:
                if c in dvals:
                    stats["d34_dv"] += 1
                elif any(occurs(c, x) for x in recs):
                    stats["d34_rec"] += 1
                else:
                    F["d34_numbers"].append("section '%s': %s is in neither the D34 result nor d34.tsv and is no [dv:] value of its block (%s): ...%s..."
                                            % (sec, c, ", ".join(dvs) or "no dv tag", norm(b)[max(0, norm(b).find(c) - 50):norm(b).find(c) + 30]))

    # ---------------------------------------------------------------- 13, 14, 15: DIAGNOSTIC 34 recomputed here
    D = None
    try:
        D = d34_recompute(ctx)
    except (SystemExit, Exception) as e:               # the committed readers stop by SystemExit
        msg = "DIAGNOSTIC 34 could not be recomputed from %s: %s" % (ctx.d34root, str(e)[:200])
        for k in ("d34_controls", "arm1_legs", "claims"):
            F[k].append(msg)
    fig_ctl = {}
    fc = Path(figdir) / "controls.tsv"
    if fc.is_file():
        for l in fc.read_text().splitlines():
            f = l.split("\t")
            if len(f) >= 8 and not l.startswith("#"):
                fig_ctl[f[0]] = f
    if D:
        # 13
        for cid, arm in (("C20", "a1"), ("C21", "a2")):
            mine = "%.4f" % D["ctl"][arm]
            if mine != D34_REF:
                F["d34_controls"].append("arm %s control barrier recomputed from %s is %s, not %s" % (arm[1], ctx.d34root.name + "/", mine, D34_REF))
            if cid not in fig_ctl:
                F["d34_controls"].append("figures/controls.tsv has no row %s for the arm %s control barrier" % (cid, arm[1]))
            elif fig_ctl[cid][2] != mine or fig_ctl[cid][7] != "PASS":
                F["d34_controls"].append("figures/controls.tsv records %s (%s) for the arm %s control barrier; recomputed here: %s"
                                         % (fig_ctl[cid][2], fig_ctl[cid][7], arm[1], mine))
        if not re.search(r"control barrier[^.\n]{0,200}%s|%s[^.\n]{0,200}control" % (re.escape(D34_REF), re.escape(D34_REF)), doc):
            F["d34_controls"].append("the document does not state the control barrier %s" % D34_REF)
        nf, worstf = 0, 0.0
        for i, n in enumerate(D["names"]):
            mine, recv = D["ddeA"][n], D["frozen"][n][1]
            nf += 1
            worstf = max(worstf, abs(mine - recv))
            if abs(mine - recv) > 5e-7:
                F["d34_controls"].append("frozen ddE of %s recomputed here is %.6f; frozen/frozen.tsv records %.6f" % (n, mine, recv))
            cid = "C07.%02d" % (i + 1)
            if cid not in fig_ctl:
                F["d34_controls"].append("figures/controls.tsv has no row %s for the frozen ddE of %s" % (cid, n))
            elif fig_ctl[cid][2] != "%.6f" % mine or fig_ctl[cid][7] != "PASS":
                F["d34_controls"].append("figures/controls.tsv records %s (%s) for the frozen ddE of %s (%s); recomputed here: %.6f"
                                         % (fig_ctl[cid][2], fig_ctl[cid][7], n, cid, mine))
        stats["d34_ctl"] = "arm 1 %.4f, arm 2 %.4f; %d frozen ddE, worst |diff| %.1e" % (D["ctl"]["a1"], D["ctl"]["a2"], nf, worstf)

        # 14
        t3 = {}
        tp = Path(figdir) / "T3_d34_arms.tsv"
        if tp.is_file():
            hdr = None
            for l in tp.read_text().splitlines():
                if l.startswith("#") or not l.strip():
                    continue
                f = l.split("\t")
                if hdr is None:
                    hdr = f
                else:
                    t3[f[0]] = dict(zip(hdr, f))
        else:
            F["arm1_legs"].append("figures/T3_d34_arms.tsv is missing")
        cols = (("rise_R", "a1_rise_R"), ("rise_TS", "a1_rise_TS"), ("dip_R", "a1_dip_R"), ("dip_TS", "a1_dip_TS"),
                ("shifted", "a1_shifted"), ("dde", "a1_ddE_kcal"))
        for n in D["names"]:
            if t3 and n not in t3:
                F["arm1_legs"].append("figures/T3_d34_arms.tsv has no row for %s" % n)
                continue
            for mk, tk in cols:
                stats["legs"] += 1
                if t3 and t3[n].get(tk) != D["a1"][n][mk]:
                    F["arm1_legs"].append("%s %s: figures/T3_d34_arms.tsv says %s; recomputed here from d34/a1_%s_*: %s"
                                          % (n, tk, t3[n].get(tk), n, D["a1"][n][mk]))
            b1 = D["rec"].get(("B1", n))
            if not b1 or abs(float(b1[5]) - float(D["a1"][n]["dde"])) > 5e-7:
                F["arm1_legs"].append("%s: d34.tsv's arm 1 ddE is %s; recomputed here: %s" % (n, b1[5] if b1 else "absent", D["a1"][n]["dde"]))
        for d in ("DV113", "DV117", "DV134"):
            stats["legs"] += 1
            if d not in dv or dv[d][1] != D["dv"][d]:
                F["arm1_legs"].append("%s is %s in figures/derived.tsv; recomputed here from d34/a1_*: %s" % (d, dv[d][1] if d in dv else "absent", D["dv"][d]))

        # 15, the numbers behind the claims
        mine = dict(D["dv"])
        t5 = rd / "T5_positive_charge_vs_verdict.tsv"
        if t5.is_file():
            rows = [l.split("\t") for l in t5.read_text().splitlines() if l.strip() and not l.startswith("#")][1:]
            two = [r for r in rows if r[0] != "D31" and r[1] != "argon ECP centre"]
            lo, hi = [r for r in two if float(r[4]) <= 0.50], [r for r in two if float(r[4]) >= 0.55]
            d31 = [r for r in rows if r[0] == "D31"]
            mine.update(DV140=str(len(two)), DV141=str(len(lo)), DV142=str(sum(1 for r in lo if r[6] == "UNRELIABLE")),
                        DV143=str(len(hi)), DV144=str(sum(1 for r in hi if r[6] != "UNRELIABLE")), DV145=str(len(d31)),
                        DV146=str(sum(1 for r in d31 if float(r[4]) >= 0.55 and r[6] == "RELIABLE")),
                        DV147=str(sum(1 for r in d31 if float(r[4]) <= 0.50 and r[6] == "UNRELIABLE")),
                        DV149=str(sum(1 for r in d31 if float(r[4]) == 0.0)))
        used = set()
        for _sec, b in blocks_of(doc):
            used |= set(tags(b)[1])
        for cid, what, ids, _rej in CLAIMS:
            for d in ids:
                stats["claim_dvs"] += 1
                if d not in mine:
                    if d in used:
                        F["claims"].append("%s (%s): the document uses %s and this script has no recomputation of it" % (cid, what, d))
                elif d not in dv or dv[d][1] != mine[d]:
                    F["claims"].append("%s (%s): %s is %s in figures/derived.tsv and %s recomputed here"
                                       % (cid, what, d, dv[d][1] if d in dv else "absent", mine[d]))
        d11 = ctx.text("R", S6 + "/standoff/DIAGNOSTIC11_RESULT.txt", figdir) or ""
        for want, got in zip(C4_NUMBERS, D["c4"]):
            stats["claim_dvs"] += 1
            if want != got or not occurs(want.lstrip("+-"), d11):
                F["claims"].append("C4 (DIAGNOSTIC 11's sign control): %s recomputed from standoff/ where the claim and D11's result have %s" % (got, want))
    # 15, the forms the recomputation rejected, and C3's place
    flat = norm(re.sub(r"`[^`\n]+`", " ", doc))
    for cid, what, _ids, rejected in CLAIMS:
        for pat in rejected:
            m = re.search(pat, flat, re.I)
            if m:
                F["claims"].append("%s (%s): a form of the claim that the recomputation rejected appears: '%s'" % (cid, what, m.group(0)[:100]))
    # ---------------------------------------------------------------- 16 the post hoc label, wherever C3 is discussed
    for sec, b in blocks_of(doc):
        if b.startswith("#") or not re.search(C3_TRIGGER, b):
            continue
        stats["c3_blocks"] += 1
        if re.match(r"10\.", sec):
            F["claims"].append("C3 is post hoc and appears in the Conclusions: ...%s..." % norm(b)[:90])
        if not re.search(r"post hoc", b, re.I):
            F["posthoc"].append("section '%s': a block discusses C3 (matched /%s/) without the words post hoc: %s..."
                                % (sec, re.search(C3_TRIGGER, b).group(0), norm(b)[:90]))

    # ================================================================ the four checks of 15 November 2026
    # ---------------------------------------------------------------- 17 the correlated barrier, recomputed here
    def cc_read(path):
        """One DLPNO-CCSD(T) output, read by this file's own code.  A label printed with two values is refused."""
        t = path.read_text(errors="replace")
        if "****ORCA TERMINATED NORMALLY****" not in t:
            return None
        o = {}
        for lab in ("E(0)", "E(CORR)(corrected)", "E(CCSD)", "Triples Correction (T)", "E(CCSD(T))"):
            v = set(re.findall(r"^" + re.escape(lab) + r"\s+\.\.\.\s+(-?\d+\.\d+)", t, re.M))
            if len(v) != 1:
                return None
            o[lab] = float(v.pop())
        o["final"] = float(re.findall(r"FINAL SINGLE POINT ENERGY\s+(-?\d+\.\d+)", t)[-1])
        return o
    d36 = {st: cc_read(ctx.res / S6 / "d36a2" / ("B_%s" % st) / "job.out") if (ctx.res / S6 / "d36a2" / ("B_%s" % st) / "job.out").is_file()
           else None for st in ("R", "TS")}
    if not d36["R"] or not d36["TS"]:
        F["benchmark"].append("d36a2/B_R or d36a2/B_TS is not a readable, normally terminated DLPNO-CCSD(T) output")
    else:
        for st in ("R", "TS"):
            if abs(d36[st]["final"] - d36[st]["E(CCSD(T))"]) > 1e-9:
                F["benchmark"].append("d36a2/B_%s: the last FINAL SINGLE POINT ENERGY is not its E(CCSD(T))" % st)
        for did, lab in sorted(BENCH_DV.items()):
            mine = "%.4f" % ((d36["TS"][lab] - d36["R"][lab]) * H2K)
            if did not in dv or dv[did][1] != mine:
                F["benchmark"].append("%s (%s part of the barrier) is %s in figures/derived.tsv; recomputed here from d36a2/: %s"
                                      % (did, lab, dv[did][1] if did in dv else "absent", mine))
        bench = "%.4f" % ((d36["TS"]["E(CCSD(T))"] - d36["R"]["E(CCSD(T))"]) * H2K)
        t36 = ctx.res / S6 / "d36.tsv"
        rec = [l.split("\t") for l in t36.read_text().splitlines() if l.startswith("B\t")] if t36.is_file() else []
        if not rec or abs(float(rec[0][4]) - float(bench)) > 6e-5:
            F["benchmark"].append("d36.tsv records %s for arm B; recomputed here: %s" % (rec[0][4] if rec else "nothing", bench))
        if not re.search(r"(?<![\d.])" + re.escape(bench) + r"(?![\d])", doc):
            F["benchmark"].append("the document does not state the correlated barrier %s" % bench)
        stats["bench"] = bench
    d38 = {rn: cc_read(ctx.res / S6 / "d38" / rn / "job.out") if (ctx.res / S6 / "d38" / rn / "job.out").is_file() else None for rn in ("N", "B1")}
    if not d38["N"] or not d38["B1"]:
        F["benchmark"].append("d38/N or d38/B1 is not a readable, normally terminated DLPNO-CCSD(T) output")
    else:
        for did, lab in sorted(D38_DV.items()):
            mine = "%+.9f" % (d38["B1"][lab] - d38["N"][lab])
            if did not in dv or dv[did][1] != mine:
                F["benchmark"].append("%s (%s, run B1 minus run N) is %s in figures/derived.tsv; recomputed here from d38/: %s"
                                      % (did, lab, dv[did][1] if did in dv else "absent", mine))

    # ---------------------------------------------------------------- 18 to 20: labels that travel with a number
    for sec, b in blocks_of(doc):
        if b.startswith("#"):
            continue
        flat = norm(b)
        if re.search(EXPERIMENT, b):
            stats["exp_blocks"] += 1
            if not re.search(r"free energy", b):
                F["experiment"].append("section '%s': a block holds the experimental 24.5 and does not say it is a free energy: %s..." % (sec, flat[:90]))
            if re.search(BENCHMARK, b) and not re.search(r"NOT AGREEMENT|not the same quantity|a different quantity", b):
                F["experiment"].append("section '%s': a block sets the correlated barrier beside the experimental 24.5 without saying they are "
                                       "different quantities: %s..." % (sec, flat[:90]))
        hit = next((x for x in A1_NUMBERS if re.search(x, b)), None)
        if hit:
            stats["a1_blocks"] += 1
            if "NOT INTERPRETABLE" not in b:
                F["d38_a1"].append("section '%s': a block quotes a number of run A1 (matched /%s/) without the words NOT INTERPRETABLE: %s..."
                                   % (sec, hit, flat[:90]))
        hit = next((x for x in POSTHOC_LEVEL if re.search(x, b)), None)
        if hit and not b.lstrip().startswith("|"):
            stats["ph_blocks"] += 1
            if not re.search(r"post hoc", b, re.I):
                F["posthoc_lvl"].append("section '%s': a block quotes a post hoc distance from the benchmark (matched /%s/) without the words "
                                        "post hoc: %s..." % (sec, hit, flat[:90]))
        hit = next((x for x in T9_NUMBERS if re.search(x, b)), None)
        if hit:
            stats["t9_blocks"] += 1
            if not (re.search(r"post hoc", b, re.I) and re.search(r"data only|not the test", b, re.I)):
                F["posthoc_lvl"].append("section '%s': a block quotes an orbital energy of a charged design (matched /%s/) without 'post hoc' and "
                                        "'data only': %s..." % (sec, hit, flat[:90]))

    # ================================================================ the four checks of the evening of 15 November 2026
    # ---------------------------------------------------------------- 21 the correlated barrier at def2-SVPD, recomputed here
    p41 = {"R": ctx.res / S6 / "d38" / "N" / "job.out", "TS": ctx.res / S6 / "d41" / "T" / "job.out"}
    d41 = {st: cc_read(q) if q.is_file() else None for st, q in p41.items()}
    if not d41["R"] or not d41["TS"]:
        F["benchmark_sv"].append("d38/N or d41/T is not a readable, normally terminated DLPNO-CCSD(T) output")
    else:
        for st, name in (("R", "d38/N"), ("TS", "d41/T")):
            if abs(d41[st]["final"] - d41[st]["E(CCSD(T))"]) > 1e-9:
                F["benchmark_sv"].append("%s: the last FINAL SINGLE POINT ENERGY is not its E(CCSD(T))" % name)
        for did, lab in sorted(D41_DV.items()):
            mine = "%.4f" % ((d41["TS"][lab] - d41["R"][lab]) * H2K)
            if did not in dv or dv[did][1] != mine:
                F["benchmark_sv"].append("%s (%s part of the def2-SVPD barrier) is %s in figures/derived.tsv; recomputed here from d38/N and d41/T: %s"
                                         % (did, lab, dv[did][1] if did in dv else "absent", mine))
        b41 = (d41["TS"]["E(CCSD(T))"] - d41["R"]["E(CCSD(T))"]) * H2K
        if d36["R"] and d36["TS"]:
            mine = "%+.4f" % ((d36["TS"]["E(CCSD(T))"] - d36["R"]["E(CCSD(T))"]) * H2K - b41)
            if "DV305" not in dv or dv["DV305"][1] != mine:
                F["benchmark_sv"].append("DV305 (the correlated barrier, def2-TZVPD minus def2-SVPD) is %s in figures/derived.tsv; recomputed here: %s"
                                         % (dv["DV305"][1] if "DV305" in dv else "absent", mine))
        t41 = ctx.res / S6 / "d41.tsv"
        rec = [l.split("\t") for l in t41.read_text().splitlines() if l.startswith("E(CCSD(T))\t")] if t41.is_file() else []
        if not rec or abs(float(rec[0][3]) - b41) > 6e-5:
            F["benchmark_sv"].append("d41.tsv records %s for the def2-SVPD barrier; recomputed here: %.4f" % (rec[0][3] if rec else "nothing", b41))
        if not re.search(r"(?<![\d.])" + re.escape("%.4f" % b41) + r"(?![\d])", doc):
            F["benchmark_sv"].append("the document does not state the def2-SVPD correlated barrier %.4f" % b41)
        stats["bench_sv"] = "%.4f" % b41

    # ---------------------------------------------------------------- 22 every detachment energy of DIAGNOSTIC 40, recomputed here
    def e_last(q, qmmm):
        """The energy of one run by DIAGNOSTIC 40's rule: the (QM/MM) total of a QM/MM run, which must be printed once;
        otherwise the last FINAL SINGLE POINT ENERGY.  None if the run did not terminate normally."""
        if not q.is_file():
            return None
        x = q.read_text(errors="replace")
        if "****ORCA TERMINATED NORMALLY****" not in x:
            return None
        if qmmm:
            v = re.findall(r"FINAL SINGLE POINT ENERGY \(QM/MM\)\s+(-?\d+\.\d+)", x)
            return float(v[0]) if len(v) == 1 else None
        v = re.findall(r"FINAL SINGLE POINT ENERGY\s+(-?\d+\.\d+)", x)
        return float(v[-1]) if v else None

    def table(name, keycol):
        q = Path(figdir) / name
        rows, cols = {}, None
        if q.is_file():
            for l in q.read_text().splitlines():
                if l.startswith("#") or not l.strip():
                    continue
                f = l.split("\t")
                if cols is None:
                    cols = f
                    continue
                r = dict(zip(cols, f))
                rows[r[keycol]] = r
        return rows
    t12, t13 = table("T12_d40_bare_functionals.tsv", "run"), table("T13_d40_arrangements.tsv", "run")
    t40 = ctx.res / S6 / "d40.tsv"
    rows40 = []
    if t40.is_file():
        ls = t40.read_text().splitlines()
        rows40 = [dict(zip(ls[0].split("\t"), l.split("\t"))) for l in ls[1:] if l.strip()]
    if not rows40:
        F["detachment"].append("d40.tsv is missing or empty")
    for r in rows40:
        rn, d_m, d_d = r["run"], ctx.res / S6 / "d40" / r["run"], ctx.res / S6 / r["source"]
        qmmm = (d_m / "combined.xyz").is_file()
        if qmmm != (d_d / "combined.xyz").is_file():
            F["detachment"].append("%s and its source %s are not the same kind of run" % (rn, r["source"]))
            continue
        em, ed = e_last(d_m / "job.out", qmmm), e_last(d_d / "job.out", qmmm)
        if em is None or ed is None:
            F["detachment"].append("%s or its source %s is not a readable, normally terminated run" % (rn, r["source"]))
            continue
        vde = (em - ed) * EV
        stats["vde"] += 1
        stats["vde_neg"] += vde <= 0
        if abs(vde - float(r["VDE_eV"])) > 6e-5:
            F["detachment"].append("%s: d40.tsv records %s eV; recomputed here from the raw outputs: %+.4f" % (rn, r["VDE_eV"], vde))
        tr = (t12 if r["part"] == "A" else t13).get(rn)
        if tr is None:
            F["detachment"].append("%s is in neither figures/T12 nor figures/T13" % rn)
            continue
        if tr["VDE_eV"] != "%+.4f" % vde:
            F["detachment"].append("%s: the figures table says %s eV; recomputed here from the raw outputs: %+.4f" % (rn, tr["VDE_eV"], vde))
        if r["part"] == "B":
            word = "BOUND" if vde > 0 else "METASTABLE"
            if tr.get("word") != word or r.get("word") != word:
                F["detachment"].append("%s: the sign gives %s; figures/T13 says %s and d40.tsv says %s" % (rn, word, tr.get("word"), r.get("word")))
    if rows40 and len(t12) + len(t13) != len(rows40):
        F["detachment"].append("figures/T12 and T13 hold %d runs and d40.tsv %d" % (len(t12) + len(t13), len(rows40)))

    # ---------------------------------------------------------------- 23 the five lowerings Probe 2's sentence reaches
    for sec, b in blocks_of(doc):
        if b.startswith("#") or not re.search(REACHED, b):
            continue
        stats["reached_blocks"] += 1
        if "METASTABLE" not in b:
            F["metastable"].append("section '%s': a block quotes a lowering that DIAGNOSTIC 40 reaches (matched %s) without the word METASTABLE: %s..."
                                   % (sec, re.search(REACHED, b).group(0), norm(b)[:90]))

    # ---------------------------------------------------------------- 24 DIAGNOSTIC 39 is committed and is not entered
    held = set()
    for k, (repo, path) in sorted(src.items()):
        if repo == "R" and re.search(D39_DATA, path):
            F["d39_held"].append("Appendix C names %s (%s): no output or table of DIAGNOSTIC 39 may be a source while it is not entered" % (k, path))
        if repo == "R" and (re.search(D39_DOCS, path) or re.search(D39_DATA, path)):
            held.add(k)
    if not re.search(r"NOT ENTERED", doc):
        F["d39_held"].append("the document does not say that DIAGNOSTIC 39's result is NOT ENTERED")
    for q in sorted(Path(figdir).glob("*.tsv")):
        if q.name == "derived.tsv" or re.match(r"T1[a-e]?_|T2_", q.name):
            continue                                   # the two run counts, and the census of every committed file
        if re.search(r"cs_stage6/d39/|(?<![\w.])d39/", q.read_text(errors="replace")):
            F["d39_held"].append("figures/%s reads or names an output of DIAGNOSTIC 39" % q.name)
    for sec, b in blocks_of(doc):
        if b.startswith("#") or "| key | repo | path |" in b:
            continue
        keys, dvs = tags(b)
        if not any(k in held for k in keys):
            continue
        stats["d39_blocks"] += 1
        others = [ctx.text(src[k][0], src[k][1], figdir) or "" for k in keys if k in src and k not in held]
        dvals = set(dv[x][1].lstrip("+-") for x in dvs if x in dv)
        items, _ex = numerals_of(b)
        for w, kind, cand in items:
            if kind == "ident":
                if not any(w in x for x in others):
                    F["d39_held"].append("section '%s': '%s' stands in a block that cites DIAGNOSTIC 39's result or amendment and is in none of "
                                         "the block's other sources" % (sec, w))
                continue
            for c in cand:
                if c not in dvals and not any(occurs(c, x) for x in others):
                    F["d39_held"].append("section '%s': %s stands in a block that cites DIAGNOSTIC 39's result or amendment and is in none of the "
                                         "block's other sources; a number of D39 may not be quoted while it is NOT ENTERED" % (sec, c))

    # ================================================================ the two checks of 16 November 2026
    # ---------------------------------------------------------------- 25 every detachment energy and ddE of DIAGNOSTIC 42, recomputed here
    t14 = table("T14_d42_conditions.tsv", "run")
    t42 = ctx.res / S6 / "d42.tsv"
    rows42 = []
    if t42.is_file():
        ls = t42.read_text().splitlines()
        rows42 = [dict(zip(ls[0].split("\t"), l.split("\t"))) for l in ls[1:] if l.strip()]
    if not rows42:
        F["bound_d42"].append("d42.tsv is missing or empty")
    ctl42 = {st: e_last(ctx.res / S6 / "d31" / ("ctl_%s" % st) / "job.out", False) for st in ("R", "TS")}
    if None in ctl42.values():
        F["bound_d42"].append("a charge-free control of DIAGNOSTIC 31 is not a readable, normally terminated run")
    words42, edi42, class42, dde42 = {}, {}, {}, {}
    for r in rows42:
        rn = r["run"]
        q_m = ctx.res / S6 / "d42" / rn / "job.out"
        if r["source"] != "d31/B/" + rn:
            F["bound_d42"].append("%s: d42.tsv names %s as its source, not the run of DIAGNOSTIC 31 of the same name" % (rn, r["source"]))
            continue
        em, ed = e_last(q_m, False), e_last(ctx.res / S6 / r["source"] / "job.out", False)
        if em is None or ed is None:
            F["bound_d42"].append("%s or its source %s is not a readable, normally terminated run" % (rn, r["source"]))
            continue
        vde = (em - ed) * EV
        stats["vde42"] += 1
        stats["vde42_neg"] += vde <= 0
        if abs(vde - float(r["VDE_eV"])) > 6e-5:
            F["bound_d42"].append("%s: d42.tsv records %s eV; recomputed here from the raw outputs: %+.4f" % (rn, r["VDE_eV"], vde))
        verd = re.findall(r"stability analysis[^\n]*?wave ?function is (stable|unstable)", q_m.read_text(errors="replace"), re.I)
        word = "METASTABLE" if vde <= 0 else ("BOUND" if verd and verd[-1].lower() == "stable" else D42_NOTSTABLE)
        words42[(r["condition"], r["state"])] = word
        edi42[(r["condition"], r["state"])] = ed
        class42[r["condition"]] = r["class"]
        tr = t14.get(rn)
        if tr is None:
            F["bound_d42"].append("%s is not in figures/T14" % rn)
            continue
        if tr["VDE_eV"] != "%+.4f" % vde:
            F["bound_d42"].append("%s: figures/T14 says %s eV; recomputed here from the raw outputs: %+.4f" % (rn, tr["VDE_eV"], vde))
        if tr.get("word") != word or r.get("word") != word:
            F["bound_d42"].append("%s: DIAGNOSTIC 42's rule gives %s; figures/T14 says %s and d42.tsv says %s" % (rn, word, tr.get("word"), r.get("word")))
    if rows42 and len(t14) != len(rows42):
        F["bound_d42"].append("figures/T14 holds %d runs and d42.tsv %d" % (len(t14), len(rows42)))
    both42 = set()
    for c in sorted(class42):
        if (c, "R") not in edi42 or (c, "TS") not in edi42 or None in ctl42.values():
            F["bound_d42"].append("%s: the two states and the two controls are not all readable, so its ddE is not recomputed" % c)
            continue
        dde = ((edi42[(c, "TS")] - ctl42["TS"]) - (edi42[(c, "R")] - ctl42["R"])) * H2K
        dde42[c] = dde
        w = (words42[(c, "R")], words42[(c, "TS")])
        if w == ("BOUND", "BOUND"):
            both42.add(c)
        want = "APPLIES" if "METASTABLE" in w else ("does not reach" if w == ("BOUND", "BOUND") else "neither")
        for st in ("R", "TS"):
            rn = "%s_%s" % (c, st)
            rr = next(x for x in rows42 if x["run"] == rn)
            if abs(float(rr["ddE_kcal"]) - dde) > 6e-5:
                F["bound_d42"].append("%s: d42.tsv records a ddE of %s kcal/mol; recomputed here from DIAGNOSTIC 31's outputs: %+.4f" % (rn, rr["ddE_kcal"], dde))
            tr = t14.get(rn)
            if tr is not None and tr["ddE_kcal"] != "%+.4f" % dde:
                F["bound_d42"].append("%s: figures/T14 says a ddE of %s kcal/mol; recomputed here from DIAGNOSTIC 31's outputs: %+.4f" % (rn, tr["ddE_kcal"], dde))
            if tr is not None and tr["probe2"] != want:
                F["bound_d42"].append("%s: figures/T14 reads Probe 2's sentence as '%s'; the two words give '%s'" % (rn, tr["probe2"], want))
    stats["cond42"], stats["both42"] = len(class42), len(both42)
    if rows42 and ("DV386" not in dv or dv["DV386"][1] != str(len(both42))):
        F["bound_d42"].append("DV386 (conditions BOUND at both states) is %s in figures/derived.tsv; recomputed here from the raw outputs: %d"
                              % (dv["DV386"][1] if "DV386" in dv else "absent", len(both42)))
    h42 = [c for c in class42 if class42[c] == "H"]
    if h42:
        fals = any("METASTABLE" in (words42.get((c, "R")), words42.get((c, "TS"))) for c in h42)
        if fals and "H9 IS NOT FALSIFIED" in doc:
            F["bound_d42"].append("the document quotes H9 IS NOT FALSIFIED, and a condition whose charges sum to zero or more is METASTABLE here")
        if not fals and "H9 IS FALSIFIED" in doc:
            F["bound_d42"].append("the document quotes H9 IS FALSIFIED, and no condition whose charges sum to zero or more is METASTABLE here")
        if not fals and all(c in both42 for c in h42) and ("H9 IS NOT FALSIFIED BY THESE %d CONDITIONS." % len(h42)) not in doc:
            F["bound_d42"].append("the document does not quote DIAGNOSTIC 42's verdict, H9 IS NOT FALSIFIED BY THESE %d CONDITIONS." % len(h42))

    # ---------------------------------------------------------------- 26 the thirteen ddE of DIAGNOSTIC 42 are NOT QUOTABLE
    nq42 = dict((c, "%.4f" % v) for c, v in dde42.items())
    reached42 = set(c for c in dde42 if "METASTABLE" in (words42.get((c, "R")), words42.get((c, "TS"))))
    ph42 = [dv[k][1] for k in D42_POSTHOC_DV if k in dv]
    if rows42 and len(ph42) != len(D42_POSTHOC_DV):
        F["not_quotable"].append("figures/derived.tsv does not hold all of %s" % ", ".join(D42_POSTHOC_DV))
    for k in D42_POSTHOC_DV:
        if k in dv and "POST HOC" not in dv[k][2]:
            F["not_quotable"].append("%s is a comparison of this write-up's own and its row in figures/derived.tsv is not marked POST HOC" % k)
    for sec, b in blocks_of(doc):
        if b.startswith("#") or "| key | repo | path |" in b:
            continue
        hit = sorted(c for c, s in nq42.items() if re.search(r"(?<![\d.])" + re.escape(s) + r"(?!\d)", b))
        if hit:
            stats["nq_blocks"] += 1
            if "NOT QUOTABLE" not in b:
                F["not_quotable"].append("section '%s': a block quotes a ddE of DIAGNOSTIC 42 (%s, %s) without the words NOT QUOTABLE: %s..."
                                         % (sec, hit[0], nq42[hit[0]], norm(b)[:90]))
            r_ = [c for c in hit if c in reached42]
            if r_ and "METASTABLE" not in b:
                F["not_quotable"].append("section '%s': a block quotes a ddE that Probe 2's sentence reaches (%s, %s) without the word METASTABLE: %s..."
                                         % (sec, r_[0], nq42[r_[0]], norm(b)[:90]))
        hitp = [s for s in ph42 if re.search(r"(?<![\d.+-])" + re.escape(s) + r"(?!\d)", b)]
        if hitp and not b.lstrip().startswith("|"):
            stats["ph42_blocks"] += 1
            if not re.search(r"post hoc", b, re.I):
                F["not_quotable"].append("section '%s': a block quotes a comparison of DIAGNOSTIC 42's numbers that is this write-up's own (%s) without "
                                         "the words post hoc: %s..." % (sec, hitp[0], norm(b)[:90]))

    # ================================================================ the two checks of 18 November 2026
    # ---------------------------------------------------------------- 27 DIAGNOSTIC 43, recomputed here from the raw outputs
    P43 = ctx.code / "phase2.2"
    audit43, d33b43 = load(P43 / "cs6_audit.py", "cswu43_audit"), load(P43 / "cs6_d33b.py", "cswu43_d33b")
    t43 = ctx.res / S6 / "d43.tsv"
    rows43 = []
    if t43.is_file():
        ls = t43.read_text().splitlines()
        rows43 = [dict(zip(ls[0].split("\t"), l.split("\t"))) for l in ls[1:] if l.strip()]
    if not rows43:
        F["survivors_d43"].append("d43.tsv is missing or empty")
    t15 = {}
    q15 = Path(figdir) / "T15_d43_conditions.tsv"
    if q15.is_file():
        cols15 = None
        for l in q15.read_text().splitlines():
            if l.startswith("#") or not l.strip():
                continue
            f = l.split("\t")
            if cols15 is None:
                cols15 = f
                continue
            r_ = dict(zip(cols15, f))
            t15[(r_["condition"], r_["state"])] = r_
    if rows43 and len(t15) != len(rows43):
        F["survivors_d43"].append("figures/T15 holds %d condition-states and d43.tsv %d" % (len(t15), len(rows43)))

    def read43(q):
        """One dianion run of DIAGNOSTIC 43 or one control of DIAGNOSTIC 31, through the committed readers."""
        if not q.is_file():
            return None
        x = q.read_text(errors="replace")
        c = audit43.read(x)
        if c is None or not c["norm"] or c["notconv"]:
            return None
        c["shifted"] = any(re.search(s_, x) for s_ in d33b43.SHIFT)
        c["stab"] = d33b43.stability(x)
        return c
    c43 = {st: read43(ctx.res / S6 / "d31" / ("ctl_%s" % st) / "job.out") for st in ("R", "TS")}
    if None in c43.values():
        F["survivors_d43"].append("a charge-free control of DIAGNOSTIC 31 is not a readable, normally terminated run")
    st43, part43 = {}, {}
    for r in rows43:
        k = (r["condition"], r["state"])
        part43[r["condition"]] = (r["part"], float(r["sum_q_e"]))
        di = read43(ctx.res / S6 / "d43" / r["di_run"] / "job.out")
        q_m = ctx.res / S6 / "d43" / r["mo_run"] / "job.out"
        em = e_last(q_m, False)
        if di is None or em is None or c43[r["state"]] is None:
            F["survivors_d43"].append("%s %s: the dianion run, the monoanion run or the control is not a readable, normally terminated run" % k)
            continue
        ctl = c43[r["state"]]
        vde = (em - di["E"]) * EV
        stats["vde43"] += 1
        stats["vde43_neg"] += vde <= 0
        verd = re.findall(r"stability analysis[^\n]*?wave ?function is (stable|unstable)", q_m.read_text(errors="replace"), re.I)
        word = "METASTABLE" if vde <= 0 else ("BOUND" if verd and verd[-1].lower() == "stable" else D42_NOTSTABLE)
        rise, dip = di["maxlw"] - ctl["maxlw"], di["dip"] / ctl["dip"]
        flag = rise >= audit43.D_LW or dip >= audit43.R_DIP or (di["stab"]["ran"] and di["stab"]["last"] == "UNSTABLE")
        held = di["shifted"] or di["gap"] is None or ctl["gap"] is None or not di["stab"]["ran"]
        if not held and di["gap"] < audit43.GAP_ABS and di["gap"] < audit43.R_GAP * ctl["gap"]:
            flag = True
        st43[k] = dict(E=di["E"], vde=vde, word=word, flag=flag, held=held, rise=rise)
        tr = t15.get(k)
        for what, rec in (("d43.tsv", r), ("figures/T15", tr)):
            if rec is None:
                F["survivors_d43"].append("%s %s is not in %s" % (k[0], k[1], what))
                continue
            if abs(float(rec["VDE_eV"]) - vde) > 6e-5:
                F["survivors_d43"].append("%s %s: %s says %s eV; recomputed here from the raw outputs: %+.4f" % (k[0], k[1], what, rec["VDE_eV"], vde))
            if rec.get("word") != word:
                F["survivors_d43"].append("%s %s: AMENDMENT 6's rule gives %s; %s says %s" % (k[0], k[1], word, what, rec.get("word")))
            if abs(float(rec["lowdin_rise_e"]) - rise) > 6e-5:
                F["survivors_d43"].append("%s %s: %s says a Loewdin rise of %s e; recomputed here: %+.4f" % (k[0], k[1], what, rec["lowdin_rise_e"], rise))
    dde43, both43, rel43, surv43 = {}, set(), set(), set()
    for c in sorted(part43):
        if (c, "R") not in st43 or (c, "TS") not in st43 or None in c43.values():
            F["survivors_d43"].append("%s: the two states and the two controls are not all readable, so it is not recomputed" % c)
            continue
        x = {st: st43[(c, st)] for st in ("R", "TS")}
        dde = ((x["TS"]["E"] - c43["TS"]["E"]) - (x["R"]["E"] - c43["R"]["E"])) * H2K
        dde43[c] = dde
        four = "UNRELIABLE" if (x["R"]["flag"] or x["TS"]["flag"]) else ("INCONCLUSIVE" if (x["R"]["held"] or x["TS"]["held"]) else "RELIABLE")
        w = (x["R"]["word"], x["TS"]["word"])
        fifth = "PASSED" if w == ("BOUND", "BOUND") else ("NOT PASSED: METASTABLE" if "METASTABLE" in w else "NOT PASSED")
        if four == "RELIABLE":
            rel43.add(c)
        if w == ("BOUND", "BOUND"):
            both43.add(c)
        s_ = part43[c][0] == "B" and four == "RELIABLE" and fifth == "PASSED"
        if s_:
            surv43.add(c)
        for st in ("R", "TS"):
            rr = next(y for y in rows43 if (y["condition"], y["state"]) == (c, st))
            tr = t15.get((c, st))
            for what, rec, ys in (("d43.tsv", rr, "yes"), ("figures/T15", tr, "SURVIVOR")):
                if rec is None:
                    continue
                if abs(float(rec["ddE_kcal"]) - dde) > 6e-5:
                    F["survivors_d43"].append("%s %s: %s says a ddE of %s kcal/mol; recomputed here from the raw outputs: %+.4f" % (c, st, what, rec["ddE_kcal"], dde))
                if rec["four_legs"] != four or rec["fifth_leg"] != fifth:
                    F["survivors_d43"].append("%s %s: read here %s on four legs and the fifth leg %s; %s says %s and %s" % (c, st, four, fifth, what, rec["four_legs"], rec["fifth_leg"]))
                if (rec["survivor"] == ys) != s_:
                    F["survivors_d43"].append("%s %s: %s marks it %s; by the readings formed here it %s a survivor" % (c, st, what, rec["survivor"], "is" if s_ else "is not"))
    stats["cond43"], stats["rel43"], stats["both43"], stats["surv43"] = len(part43), len(rel43), len(both43), len(surv43)
    if rows43 and ("DV436" not in dv or dv["DV436"][1] != str(len(surv43))):
        F["survivors_d43"].append("DV436 (the survivors of the search) is %s in figures/derived.tsv; recomputed here from the raw outputs: %d"
                                  % (dv["DV436"][1] if "DV436" in dv else "absent", len(surv43)))
    if len(dde43) == len(part43) and part43:
        pa = [c for c in part43 if part43[c][0] == "A"]
        for name, cs in (("H3", pa), ("H2", [c for c in pa if part43[c][1] < 0]), ("H9", [c for c in part43 if part43[c][1] > -5e-7 and c in rel43])):
            ok = all(c in rel43 for c in cs) if name != "H9" else all(c in both43 for c in cs)
            said = "%s IS NOT FALSIFIED BY THESE %d CONDITIONS." % (name, len(cs))
            if ok and said not in doc:
                F["survivors_d43"].append("the document does not quote DIAGNOSTIC 43's verdict, %s" % said)
            if not ok and re.search(name + r" IS NOT FALSIFIED BY THESE \d+ CONDITIONS", doc) and name != "H9":
                F["survivors_d43"].append("the document quotes %s IS NOT FALSIFIED, and a condition of its test is not RELIABLE here" % name)
            if not ok and name == "H9" and ("H9 IS NOT FALSIFIED BY THESE %d CONDITIONS." % len(cs)) in doc:
                F["survivors_d43"].append("the document quotes H9 IS NOT FALSIFIED for DIAGNOSTIC 43, and a RELIABLE condition whose charges sum to zero or more is not BOUND at both states here")
    d44 = ctx.res / S6 / "d44"
    n44 = len(list(d44.rglob("job.out"))) if d44.is_dir() else 0
    stats["out44"] = n44
    if "DV470" not in dv or dv["DV470"][1] != str(n44):
        F["survivors_d43"].append("DV470 (the outputs of DIAGNOSTIC 44 in the record) is %s in figures/derived.tsv; counted here: %d"
                                  % (dv["DV470"][1] if "DV470" in dv else "absent", n44))
    if n44 and re.search(r"D44[^.]*?(?:has no output|holds no output)", doc):
        F["survivors_d43"].append("the document says DIAGNOSTIC 44 has no output in the record, and d44/ holds %d job.out files" % n44)

    # ---------------------------------------------------------------- 28 a ddE of DIAGNOSTIC 43 is a working-level number, or is reached by Probe 2's sentence
    wl43 = dict((c, "%.4f" % v) for c, v in dde43.items() if c in rel43 and c in both43)
    p243 = dict((c, "%.4f" % v) for c, v in dde43.items() if "METASTABLE" in (st43[(c, "R")]["word"], st43[(c, "TS")]["word"]))
    for sec, b in blocks_of(doc):
        if b.startswith("#") or "| key | repo | path |" in b:
            continue
        hit = sorted(c for c, s in wl43.items() if re.search(r"(?<![\d.])" + re.escape(s) + r"(?!\d)", b))
        if hit:
            stats["wl_blocks"] += 1
            if "working-level number" not in b:
                F["working_level"].append("section '%s': a block quotes a ddE of DIAGNOSTIC 43 (%s, %s) without the words working-level number: %s..."
                                          % (sec, hit[0], wl43[hit[0]], norm(b)[:90]))
        hit = sorted(c for c, s in p243.items() if re.search(r"(?<![\d.])" + re.escape(s) + r"(?!\d)", b))
        if hit:
            stats["p2_blocks"] += 1
            if "Probe 2" not in b and "METASTABLE" not in b:
                F["working_level"].append("section '%s': a block quotes a ddE that Probe 2's sentence reaches (%s, %s) without naming Probe 2 or the word "
                                          "METASTABLE: %s..." % (sec, hit[0], p243[hit[0]], norm(b)[:90]))
    return F, stats


def report(F, stats, out=sys.stdout):
    order = ["numbers", "quotes", "citations", "withdrawn", "selflabels", "controls", "counts", "not_tested", "regeneration", "heads",
             "exp_heads", "d34_numbers", "d34_controls", "arm1_legs", "claims", "posthoc",
             "benchmark", "experiment", "d38_a1", "posthoc_lvl",
             "benchmark_sv", "detachment", "metastable", "d39_held", "bound_d42", "not_quotable", "survivors_d43", "working_level"]
    notes = {
        "numbers": "%d sourced in a named committed file, %d derived (figures/derived.tsv), %d estimated; %d items not treated as numerals"
                   % (stats["sourced"], stats["derived"], stats["estimated"], stats["exempt"]),
        "quotes": "%d quotations" % stats["quotes"], "citations": "%d files in Appendix C" % stats["files"],
        "withdrawn": "%d registry entries" % len(WITHDRAWN), "selflabels": "%d registry entries" % len(SELF_LABELS),
        "controls": "%d recomputed; worst: %s" % (stats["controls"], stats["worst"]),
        "counts": "%d run counts" % stats["counts"], "not_tested": "%d names" % stats["names"],
        "regeneration": "%d tables regenerated" % stats["tsv"], "heads": "two repositories",
        "exp_heads": "results %s, code %s" % (EXPECTED_HEADS["chorismate-thesis-results"], EXPECTED_HEADS["chorismate-thesis-code"]),
        "d34_numbers": "%d blocks cite D34: %d numerals in its result or d34.tsv, %d derived" % (stats["d34_blocks"], stats["d34_rec"], stats["d34_dv"]),
        "d34_controls": stats["d34_ctl"] or "not recomputed",
        "arm1_legs": "%d values recomputed from d34/a1_*" % stats["legs"],
        "claims": "%d numbers behind C1 to C6 recomputed" % stats["claim_dvs"],
        "posthoc": "%d blocks discuss C3" % stats["c3_blocks"],
        "benchmark": ("DLPNO-CCSD(T) barrier %s recomputed from d36a2/; two shifts recomputed from d38/" % stats["bench"]) if stats["bench"] else "not recomputed",
        "experiment": "%d blocks hold the experimental 24.5" % stats["exp_blocks"],
        "d38_a1": "%d blocks quote a number of run A1" % stats["a1_blocks"],
        "posthoc_lvl": "%d blocks quote a post hoc distance; %d an orbital energy of a charged design" % (stats["ph_blocks"], stats["t9_blocks"]),
        "benchmark_sv": ("DLPNO-CCSD(T)/def2-SVPD barrier %s recomputed from d38/N and d41/T" % stats["bench_sv"]) if stats["bench_sv"] else "not recomputed",
        "detachment": "%d detachment energies recomputed from d40/ and their committed dianion runs; %d at or below zero" % (stats["vde"], stats["vde_neg"]),
        "metastable": "%d blocks quote a lowering that DIAGNOSTIC 40 reaches" % stats["reached_blocks"],
        "d39_held": "%d blocks cite DIAGNOSTIC 39's result or amendment; none quotes a number of it" % stats["d39_blocks"],
        "bound_d42": "%d detachment energies recomputed from d42/ and DIAGNOSTIC 31's dianion runs; %d at or below zero; %d of %d conditions "
                     "BOUND at both states" % (stats["vde42"], stats["vde42_neg"], stats["both42"], stats["cond42"]),
        "not_quotable": "%d blocks quote a ddE of DIAGNOSTIC 42; %d a post hoc comparison of its numbers" % (stats["nq_blocks"], stats["ph42_blocks"]),
        "survivors_d43": "%d detachment energies recomputed from d43/; %d at or below zero; of %d conditions %d RELIABLE on four legs and %d BOUND at "
                         "both states; %d survivors; %d outputs of DIAGNOSTIC 44 in the record"
                         % (stats["vde43"], stats["vde43_neg"], stats["cond43"], stats["rel43"], stats["both43"], stats["surv43"], stats["out44"]),
        "working_level": "%d blocks quote a ddE of DIAGNOSTIC 43 that is a working-level number; %d one that Probe 2's sentence reaches"
                         % (stats["wl_blocks"], stats["p2_blocks"])}
    bad = 0
    for i, k in enumerate(order, 1):
        if i == 11:
            print("   -- added 11 November 2026 for DIAGNOSTIC 34", file=out)
        if i == 17:
            print("   -- added 15 November 2026 for DIAGNOSTICS 35 to 38", file=out)
        if i == 21:
            print("   -- added on the evening of 15 November 2026 for DIAGNOSTICS 40 and 41 and DIAGNOSTIC 37's Amendment 1", file=out)
        if i == 25:
            print("   -- added 16 November 2026 for DIAGNOSTIC 42", file=out)
        if i == 27:
            print("   -- added 18 November 2026 for DIAGNOSTIC 43 and AMENDMENTS 6 and 7", file=out)
        print("%2d %-12s %s   (%s)" % (i, k, "PASS" if not F[k] else "FAIL  %d" % len(F[k]), notes[k]), file=out)
        for msg in F[k][:40]:
            print("      - " + msg, file=out)
        bad += len(F[k])
    print("RESULT: %s" % ("ALL TWENTY-EIGHT CHECKS PASS" if not bad else "%d FAILURE(S)" % bad), file=out)
    return bad


def selftest(doc, figdir, ctx):
    """Negative tests.  Each perturbs one thing on a copy and requires the named check to fail for the named reason."""
    tests = []

    def sub(old, new, count=1):
        def f(d, fd):
            if old not in d:
                raise RuntimeError("selftest target not found: %r" % old[:50])
            return d.replace(old, new, count), fd
        return f

    def in_section(heading, old, new):
        def f(d, fd):
            i = d.index("\n## " + heading)
            j = d.find("\n## ", i + 4)
            seg = d[i:j]
            if old not in seg:
                raise RuntimeError("selftest target not found in section %s: %r" % (heading, old))
            return d[:i] + seg.replace(old, new) + d[j:], fd
        return f

    def tsv(name, old, new):
        def f(d, fd):
            p = Path(fd) / name
            t = p.read_text()
            if old not in t:
                raise RuntimeError("selftest target not found in %s: %r" % (name, old))
            p.write_text(t.replace(old, new, 1))
            return d, fd
        return f

    tests.append(("a number", "numbers", "26.8348", in_section("6. What", "26.8248 at CAM-B3LYP", "26.8348 at CAM-B3LYP")))
    tests.append(("a quotation (full stop to comma)", "quotes", "interpretable,",
                  sub("“No design measured inside r_safe is interpretable.”", "“No design measured inside r_safe is interpretable,”")))
    tests.append(("a run count in the text", "numbers", "67", sub("D32 is 66.", "D32 is 67.")))
    tests.append(("a run count in derived.tsv", "counts", "DV72", tsv("derived.tsv", "DV72\t66\t", "DV72\t67\t")))
    tests.append(("a citation", "citations", "not tracked", sub("STAGE6_DIAGNOSTIC27_RESULT_20261105.txt |", "STAGE6_DIAGNOSTIC27_RESULT_20261106.txt |")))
    tests.append(("a withdrawn figure used as evidence", "withdrawn", "0.58",
                  in_section("6. What", "Each item states", "A budget of 0.58 applies. Each item states")))
    tests.append(("a self-label", "selflabels", "PARTIAL", in_section("4. The approaches", "PARTIAL", "mixed")))
    tests.append(("a figure's tsv", "regeneration", "F5_distance_ladder.tsv", tsv("F5_distance_ladder.tsv", "+0.3994", "+0.3894")))
    tests.append(("a HEAD hash", "heads", "records", sub("`ee283016", "`fe283016")))
    tests.append(("the not-tested list", "not_tested", "CAM-B3LYP", sub("`PBE0`,", "`PBE0`, `CAM-B3LYP`,")))

    # ---- the six checks added on 11 November 2026
    def raw_control(d, fd):
        """A copy of d34/ in which the arm 1 control's TS energy is 1e-4 Eh higher, in job.out and job.engrad alike."""
        real, alt = ctx.res / S6 / "d34", Path(fd).parent / "d34"
        alt.mkdir()
        for x in sorted(real.iterdir()):
            if x.name != "a1_S6-C_TS":
                os.symlink(x, alt / x.name)
        shutil.copytree(real / "a1_S6-C_TS", alt / "a1_S6-C_TS")
        old = re.findall(r"FINAL SINGLE POINT ENERGY\s+(-?\d+\.\d+)", (alt / "a1_S6-C_TS" / "job.out").read_text(errors="replace"))[-1]
        new = "%.12f" % (float(old) + 1e-4)
        for name in ("job.out", "job.engrad"):
            q = alt / "a1_S6-C_TS" / name
            t = q.read_text(errors="replace")
            if old not in t:
                raise RuntimeError("selftest target not found in %s: %s" % (name, old))
            q.write_text(t.replace(old, new))
        ctx.d34root = alt
        return d, fd

    tests.append(("the commit the revision was built from", "exp_heads", "does not state", sub("`0ede870f6", "`fede870f6")))
    tests.append(("a D34 number that another cited file holds", "d34_numbers", "14.3612",
                  sub("-4.1567 and -5.3701, all single negative spheres. Scope: inside r_safe; not a catalytic claim; no arrangement of more "
                      "than one sphere that the optimiser placed is quotable. [src: D34R, D34T]",
                      "-4.1567 and -14.3612, all single negative spheres. Scope: inside r_safe; not a catalytic claim; no arrangement of more "
                      "than one sphere that the optimiser placed is quotable. [src: D34R, D34T, D33B]"), "numbers"))
    tests.append(("a raw control output of d34/", "d34_controls", "arm 1 control barrier recomputed", raw_control))
    tests.append(("a recorded frozen ddE", "d34_controls", "frozen ddE of S6-P",
                  tsv("controls.tsv", "C07.01\tddE(frozen) S6-P\t-7.675715", "C07.01\tddE(frozen) S6-P\t-7.675815")))
    tests.append(("an arm 1 leg", "arm1_legs", "a1_rise_R", tsv("T3_d34_arms.tsv", "\t1.6692\t", "\t1.7692\t")))
    tests.append(("the number behind a claim", "claims", "C6", tsv("derived.tsv", "DV124\t2\t", "DV124\t3\t")))
    tests.append(("a claim the recomputation rejected", "claims", "rejected",
                  in_section("6. What", "This is not the first stability correction in the record.",
                             "In 8 of the 10 arrangements the first solution is unstable. This is not the first stability correction in the record.")))
    tests.append(("C3 stated as a conclusion", "claims", "appears in the Conclusions",
                  in_section("10. Conclusions", "Scope: one frame; r_safe is one axis from one site.",
                             "No arrangement with no positive charge relocates. Scope: one frame; r_safe is one axis from one site.")))
    tests.append(("the post hoc label", "posthoc", "without the words post hoc",
                  in_section("11. Hypotheses", "*Against, and limiting; this paragraph too is post hoc.*", "*Against, and limiting.*")))
    # ---- the four checks added on 15 November 2026, and the widening of check 8
    tests.append(("the benchmark's derived value", "benchmark", "DV220", tsv("derived.tsv", "DV220\t23.8603\t", "DV220\t23.8703\t")))
    tests.append(("the denial beside the experimental number", "experiment", "different quantities",
                  in_section("6. What", "a different quantity, and D36 says “ITS CLOSENESS TO 24.5 IS NOT AGREEMENT.”", "and the two are close.")))
    tests.append(("the experimental number without its kind", "experiment", "does not say it is a free energy",
                  in_section("6. What", "the record's 24.5 kcal/mol is a free energy in water", "the record's 24.5 kcal/mol is measured in water")))
    tests.append(("the label of run A1", "d38_a1", "without the words NOT INTERPRETABLE",
                  in_section("3. Levels", "Run A1 is “NOT INTERPRETABLE” under a rule", "Run A1 is set aside under a rule")))
    tests.append(("the post hoc label of a distance", "posthoc_lvl", "post hoc distance",
                  in_section("6. What", "Set against the benchmark they are post hoc, by D36's own label:", "Set against the benchmark:")))
    tests.append(("the data-only label of an orbital energy", "posthoc_lvl", "orbital energy of a charged design",
                  in_section("6. What", "from -0.1563 to -1.6998 eV. The sign is the reading; the size is not.",
                             "from -0.1563 to -1.6998 eV, the orbital lying at +0.7782 eV. The sign is the reading; the size is not.")))
    tests.append(("a method run under another spelling", "not_tested", "M06-2X", sub("`PBE0`,", "`PBE0`, `M06-2X`,")))
    # ---- the four checks added on the evening of 15 November 2026, and the widening of check 20
    tests.append(("the def2-SVPD benchmark's derived value", "benchmark_sv", "DV300", tsv("derived.tsv", "DV300\t22.8968\t", "DV300\t22.9068\t")))
    tests.append(("a detachment energy in a table", "detachment", "B_ADD-s540_qm06500_R",
                  tsv("T13_d40_arrangements.tsv", "\t-0.6414\tMETASTABLE\t", "\t+0.6414\tBOUND\t")))
    tests.append(("the word METASTABLE beside a reached lowering", "metastable", "without the word METASTABLE",
                  in_section("10. Conclusions", "D40 has since found the dianion METASTABLE with that sphere present.",
                             "D40 has since measured the dianion with that sphere present.")))
    tests.append(("a number of D39 quoted under its citation", "d39_held", "may not be quoted",
                  sub("6 inputs and 6 outputs are in the repository (derived counts).",
                      "6 inputs and 6 outputs are in the repository (derived counts), and its first route gives +1.7702 eV. [src: D39R]"), "numbers"))
    tests.append(("the post hoc label of a distance at def2-SVPD", "posthoc_lvl", "post hoc distance",
                  in_section("2. The reference barrier", "is post hoc, as are these counts: D41 named eight functionals and it is not one of them, and",
                             "was not named by D41, and")))
    # ---- the two checks added on 16 November 2026
    tests.append(("a detachment energy of D42 in its table", "bound_d42", "n10_q0.25_R",
                  tsv("T14_d42_conditions.tsv", "\t+0.1591\tBOUND\t", "\t-0.1591\tMETASTABLE\t")))
    tests.append(("the count of conditions bound at both states", "bound_d42", "DV386", tsv("derived.tsv", "DV386\t10\t", "DV386\t11\t")))
    tests.append(("the words NOT QUOTABLE beside a ddE of D42", "not_quotable", "without the words NOT QUOTABLE",
                  in_section("10. Conclusions", "-16.3574 kcal/mol, and those are NOT QUOTABLE: they are test conditions",
                             "-16.3574 kcal/mol: they are test conditions")))
    tests.append(("the word METASTABLE beside a reached ddE of D42", "not_quotable", "without the word METASTABLE",
                  in_section("5. The reliability criterion", "For the 3 conditions with a METASTABLE state Probe 2's sentence applies",
                             "For the 3 conditions with a state at or below zero Probe 2's sentence applies")))
    tests.append(("the post hoc label of a comparison of D42's numbers", "not_quotable", "without the words post hoc",
                  in_section("5. The reliability criterion", "which is how H9 was worded. Post hoc.", "which is how H9 was worded.")))
    # ---- the two checks added on 18 November 2026
    tests.append(("a detachment energy of D43 in its table", "survivors_d43", "B_n02_q1.00 R",
                  tsv("T15_d43_conditions.tsv", "\t+0.2272\tBOUND\t", "\t-0.2272\tMETASTABLE\t")))
    tests.append(("the count of survivors", "survivors_d43", "DV436", tsv("derived.tsv", "DV436\t6\t", "DV436\t5\t")))
    tests.append(("the count of outputs of D44", "survivors_d43", "DV470", tsv("derived.tsv", "DV470\t0\t", "DV470\t1\t")))
    tests.append(("the words working-level number beside a survivor's ddE", "working_level", "without the words working-level number",
                  in_section("6. What", "against a charge-free barrier of 26.8248, are working-level numbers: Amendment 7 allows",
                             "against a charge-free barrier of 26.8248, are lowerings: Amendment 7 allows")))
    tests.append(("the mention of Probe 2 beside a ddE of D43 that its sentence reaches", "working_level", "without naming Probe 2",
                  in_section("7. What has been withdrawn", "What D43 adds to the numbers Probe 2's sentence reaches is the ddE",
                             "What D43 adds to the record is the ddE")))
    fired = 0
    lines = []

    def say(x):
        print(x)
        lines.append(x)
    say("NEGATIVE TESTS. Each line below is a fault planted on purpose, on a copy; nothing real is changed.")
    say("   'caught' is the good outcome: the check named refused the planted fault.")
    d34root = ctx.d34root
    for name, check, reason, fn, *quiet in tests:
        tmp = Path(tempfile.mkdtemp(prefix="cswu_neg_"))
        try:
            fd = tmp / "figures"
            shutil.copytree(figdir, fd, ignore=shutil.ignore_patterns("__pycache__"))
            d2, _ = fn(doc, fd)
            F, _s = run_checks(d2, fd, ctx)
            hit = [m for m in F[check] if reason in m]
            ok = bool(hit) and not any(F[q] for q in quiet)        # 'quiet' names a check that must not have caught it
            fired += ok
            say("   %-44s -> check '%s' %s%s%s" % (name, check, "caught" if ok else "NOT CAUGHT",
                                                    (": " + hit[0][:110]) if ok else "",
                                                    ("  [check '%s' let it through]" % quiet[0]) if ok and quiet else ""))
        finally:
            ctx.d34root = d34root
            shutil.rmtree(tmp, ignore_errors=True)
    say("   %d of %d planted faults were caught%s" % (fired, len(tests), ", as they should be" if fired == len(tests) else ": THE CHECKER HAS A FAULT"))
    return fired == len(tests), lines


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--res", required=True)
    ap.add_argument("--code", required=True)
    ap.add_argument("--doc")
    ap.add_argument("--figures")
    ap.add_argument("--scoped", help="validate the scoped report at this path, and not the consolidated write-up")
    ap.add_argument("--scoped-dir", help="the scoped report's own directory of tables (scoped/)")
    ap.add_argument("--full", action="store_true", help="also re-run sweep_pass2.py (about four minutes)")
    ap.add_argument("--selftest", action="store_true", help="run the negative tests after the clean run")
    ap.add_argument("--provenance", help="write the provenance record to this file")
    a = ap.parse_args()
    ctx = Ctx(a.res, a.code, a.full)
    if a.scoped:
        if not a.scoped_dir:
            ap.error("--scoped needs --scoped-dir")
        sdoc = Path(a.scoped).read_text()
        F, st = run_scoped(sdoc, a.scoped_dir, ctx)
        print("cs_writeup_check.py --scoped on %s" % Path(a.scoped).name)
        ok = report_scoped(F, st) == 0
        neg = []
        if a.selftest:
            okn, neg = selftest_scoped(sdoc, a.scoped_dir, ctx)
            ok = okn and ok
        if a.provenance:
            import io
            buf = io.StringIO()
            report_scoped(F, st, buf)
            lines = ["PROVENANCE of %s" % Path(a.scoped).name, "",
                     "chorismate-thesis-results  main           HEAD %s" % git(ctx.res, "rev-parse", "HEAD").stdout.strip(),
                     "chorismate-thesis-code     tier1-realism  HEAD %s" % git(ctx.code, "rev-parse", "HEAD").stdout.strip(), "",
                     "sha256 of every script used",
                     "   %s  scoped/make_scoped.py" % sha(Path(a.scoped_dir) / "make_scoped.py"),
                     "   %s  code: phase2.2/cs_writeup_check.py" % sha(Path(__file__).resolve())]
            lines += ["   %s  code: phase2.2/%s  (committed reader, imported)" % (sha(ctx.code / "phase2.2" / n), n) for n in ("cs6_analyse.py", "cs6_audit.py", "cs6_d33b.py")]
            lines += ["", "validation output", ""] + ["   " + l for l in buf.getvalue().rstrip().splitlines()]
            if neg:
                lines += [""] + ["   " + l for l in neg]
            Path(a.provenance).write_text("\n".join(lines) + "\n")
        return 0 if ok else 1
    if not a.doc or not a.figures:
        ap.error("--doc and --figures are needed unless --scoped is given")
    doc = Path(a.doc).read_text()
    F, stats = run_checks(doc, a.figures, ctx, a.doc)
    print("cs_writeup_check.py on %s" % Path(a.doc).name)
    bad = report(F, stats)
    ok = bad == 0
    neg = []
    if a.selftest:
        okn, neg = selftest(doc, a.figures, ctx)
        ok = okn and ok
    if a.provenance:
        import io
        buf = io.StringIO()
        report(F, stats, buf)
        lines = ["PROVENANCE of %s" % Path(a.doc).name, "",
                 "chorismate-thesis-results  main           HEAD %s" % git(ctx.res, "rev-parse", "HEAD").stdout.strip(),
                 "chorismate-thesis-code     tier1-realism  HEAD %s" % git(ctx.code, "rev-parse", "HEAD").stdout.strip(), "",
                 "sha256 of every script used"]
        for label, p in ([("figures/" + n, Path(a.figures) / n) for n in ("make_figures.py", "sweep_pass1.py", "sweep_pass2.py")] +
                         [("code: phase2.2/cs_writeup_check.py", Path(__file__).resolve())] +
                         [("code: phase2.2/" + n + "  (committed reader, imported)", ctx.code / "phase2.2" / n)
                          for n in ("cs6_analyse.py", "cs6_audit.py", "cs6_d33b.py")]):
            lines.append("   %s  %s" % (sha(p), label))
        lines += ["", "validation output", ""] + ["   " + l for l in buf.getvalue().rstrip().splitlines()]
        if neg:
            lines += [""] + ["   " + l for l in neg]
        Path(a.provenance).write_text("\n".join(lines) + "\n")
    if ctx.regen:
        shutil.rmtree(ctx.regen[0], ignore_errors=True)
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
