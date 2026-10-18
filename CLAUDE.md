# CLAUDE.md

Context notes for the bistable-oscillator geodynamo series (Hoyng, Schmitt & Ossendrijver, 2001–2002).

## Read the markdown, not the PDFs

- The three papers in `papers/` are the source of truth, but do **not** read them for routine work. Read the verified resources in `resources/` instead. They contain every equation, parameter, table, figure caption and result with journal page references, and each item was checked against the PDFs (see `resources/VERIFICATION.md`).
- `papers/` is git-ignored and is usually absent from a clone. A missing `papers/` directory is not a problem; do not try to fetch or reconstruct the PDFs.
- Open a PDF only to settle a specific suspected transcription error, and only if the file is present locally. If that leads to a change in a resource file, add an entry to `resources/VERIFICATION.md` saying what was checked and what changed.

## Reading order

1. `resources/README.md`: index, what each paper contributes, and a "which file to read for what" table.
2. `resources/00-series-overview-and-shared-notation.md`: the common model, the notation concordance (the papers use different symbols for the same quantity), and the headline numbers side by side.
3. The paper file you need:
   - `resources/hoyng-2001-gafd-geodynamo-bistable-oscillator.md`, Paper 1 (HOS): the full model, simulations and Fokker–Planck derivation.
   - `resources/schmitt-2001-pepi-bistable-geodynamo-reversals.md`, Paper 2 (SOH): the six-page summary of the method.
   - `resources/hoyng-2002-pepi-variability-geomagnetic-dipole.md`, Paper 3 (HSO): the quantitative comparison with the Sint-800 record.

Every paper file has the same 15 sections (model §4, methods §5, results §6, figures §7, implementation guide §11, glossary §12). The files are long, so load only the sections you need.

## Conventions in the resource files

- Equation numbers are the papers' own. Page numbers are journal pages; the frontmatter of each file gives the PDF-page mapping.
- Sentences beginning **Note:** are the summariser's commentary. Everything else is the paper's own content. Keep that separation when editing.
- `[uncertain]` marks values read off figure axes. Typesetting slips in the originals are transcribed as printed and flagged with a Note, never silently corrected.
- Notation differs between papers (turbulent diffusivity is β in Papers 1 and 3 but η in Paper 2; the forcing parameter is written differently in each). State which paper's notation you are using; the overview file has the conversion table.

## Git

- `papers/` is ignored and must never be committed, not even with `git add -f`. It holds the publishers' PDFs. This must stay true when the repository is pushed to GitHub. The extracted text in `resources/raw-text/` is tracked.
- Do not rename or move the PDFs: `source_pdf` in each resource's frontmatter refers to them relative to the repository root.
