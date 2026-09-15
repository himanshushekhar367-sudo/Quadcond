# Contributing

## Setting up

```bash
python -m pip install -e ".[dev]"
npm --prefix web install
```

Model and atlas binaries are not in the repository. `quadcond assets status`
reports what is missing and `quadcond assets fetch` retrieves and verifies it.
The test suite does not need them.

## Before opening a pull request

```bash
python -m pytest tests -q
ruff check quadcond tests
npm --prefix web run typecheck
npm --prefix web test
npm --prefix web run build
```

## House rules

These are the conventions the codebase is built on. A change that breaks one
of them will be asked to change, however small it is.

1. **A claim needs a basis.** Every head carries an `evidence_tier`,
   `target_semantics` and `biophysically_grounded` in the claim registry. A new
   head, or a new output, arrives with those fields filled in and with
   `docs/CLAIMS.md` updated in the same change.

2. **Refusal is a result, not an error.** Outside its applicability domain a
   head refuses with a reason. Do not substitute a default, a zero, or an
   extrapolation. On the frontend, refusal and score are separate members of a
   discriminated union so that reading one as the other does not compile.

3. **Do not invent a delta.** When one side of a comparison has no applicable
   motif, report that. Do not subtract two numbers that do not describe the
   same thing.

4. **Exact joins only.** Variants are matched on
   `(chromosome, position, ref, alt)`. Nearest-position matching is not an
   acceptable fallback.

5. **Provenance travels with the result.** Anything that leaves the process —
   a JSON response, a CSV, a report — carries the model file hash, the
   manifest comparison and the environment. CSV exports put the column names on
   line 1.

6. **Assets are verified before they are opened.** A hash mismatch is an error,
   not a warning, and an asset with neither a recorded hash nor a content
   fingerprint is refused rather than installed.

7. **Say what the code does.** Comments explain why a decision was made, not
   what the next line says. Documentation states limits before capabilities.

## Reporting a problem

Open an issue with the command you ran, the output of `quadcond info`, and the
provenance block from the result. Those three pieces identify the build, the
assets and the environment, which is usually enough to reproduce it.
