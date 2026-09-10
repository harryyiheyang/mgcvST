# Project instructions

## Source and Git synchronization

- This Git checkout is the package development and release source. The separate `magicST/mgcvST` directory contains historical research source and must not be copied over this repository as a whole.
- The author works on multiple computers. Before editing, fetch and reconcile the remote branch with local work. Preserve the author's changes; do not force-push to resolve divergence.
- After completing and validating an authorized package update, commit and push it during the same task. The author has provided standing authorization for this workflow. Report the pushed commit and local checkout path; explicitly report any push failure.
- Preserve running research jobs, saved input data, and historical simulation source unless the author explicitly requests their replacement.

## R code

- Preserve the author's successfully run analysis semantics. Keep research scripts flat and debuggable; use short variable names and minimal comments. Do not introduce error suppression or version-suffix backup files.
- Package helpers may implement reusable numerical or API mechanics. Prefer CppMatrix for supported matrix operations. Validate user inputs with clear errors.

## Academic writing

- Before drafting or revising academic prose, read `C:/Users/yxy1234/Documents/Y_Yang_Academic_Writing_Style_Guide.md` and `C:/Users/yxy1234/Downloads/SKILL.md` when available. The personal guide and the author's explicit instructions take precedence.
- Preserve functional redundancy, explicit statistical referents and separate discourse moves. Use `we` for intellectual actions and natural passive voice for procedures.
- Lead with what the evidence establishes and the conditions under which it applies. State material validity and reproducibility boundaries accurately without organizing prose around defensive disclaimers.
