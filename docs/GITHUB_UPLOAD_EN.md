# Uploading the Repository to GitHub

[English run guide](../README_EN.md) | [Chinese version](GITHUB_UPLOAD.md)

## Files to Upload

Place the contents of `paper1-from-public` at the repository root. Visitors should see `README.md`, `README_EN.md`, `run.py`, `analysis/`, and `pipeline/` directly on the repository homepage.

Upload the source directory, not only a ZIP file. A ZIP can additionally be attached to a GitHub Release.

The source package excludes raw data, donor-level data objects, aggregated matrices, installed R libraries, caches, logs, temporary test files, and obsolete archives. `.gitignore` excludes default output locations; keep custom work directories outside the repository or explicitly ignore them before committing.

Include all eight English documents in `docs/`, alongside the Chinese originals. `README_EN.md` links to the English versions.

## Browser Upload

1. Create an empty GitHub repository, for example `paper1-model-adjustment`.
2. Use the upload-files option to add the unpacked files and subdirectories to the repository root.
3. Include hidden `.gitignore` and `.gitattributes` files, configuration, environment files, tests, and documentation.
4. Check the resulting structure before committing.

Browser drag-and-drop can omit hidden files. Terminal-based upload is often easier to verify.

## Linux Terminal Upload

Run these commands in a new local repository containing `run.py`. Create the empty remote repository first.

```bash
git init -b main
git add .
git diff --cached --stat
git status --short
```

Inspect the staged files. They should contain source, configuration, tests, and documentation, not `raw`, `prepared`, `results`, `R_library`, or private credentials.

After inspection, substitute the real repository URL:

```bash
git commit -m "Release Paper1 analysis v1.2.2"
git remote add origin https://github.com/YOUR_ACCOUNT/YOUR_REPOSITORY.git
git push -u origin main
```

For an existing repository, copy the intended updates into its checkout, inspect `git diff` and `git status`, then make a normal commit and push. Do not force-push or initialize over unrelated history merely to upload this package.

## English Homepage

GitHub normally displays the root `README.md`. This package keeps that Chinese file unchanged and provides `README_EN.md` separately. To expose both languages from the homepage, the authors can add a link such as:

```markdown
[English installation and execution guide](README_EN.md)
```

An English homepage can instead use the English content as the root README, with a link to a separately retained Chinese version. Update relative links and release checksums if changing file names or content.

## Releases, Citations, and Licensing

After checking that the tag is unused, create a `v1.2.2` tag and release for the intended package. Do not move an already published tag to different contents; publish subsequent documentation changes as an ordinary commit or a separately identified documentation release.

`CITATION.cff` supplies author and software-version information. Add the actual repository URL and archival DOI when available. If using Zenodo, obtain the DOI through a real archive release before adding it to the manuscript or citation metadata; a GitHub URL is not a DOI.

The authors must select a software license. Dataset usage terms and software licensing are separate. These instructions do not create or publish a remote repository automatically.

## Checksums and Later Changes

`release_manifest.json` lists packaged-file SHA-256 values and excludes itself. `.gitattributes` preserves the byte representation of the validated analysis source by disabling automatic line-ending conversion for those files.

After changing packaged files, update the manifest. Keep analysis-version changes separate from documentation-only updates. When altering analysis source or settings, use a new result directory rather than replacing code underneath a completed or resumable run whose records identify an older implementation.
