# Conda packaging

The conda recipe is `recipe/meta.yaml` (conda-build format). It follows bioconda's conventions so it can be submitted there as-is, apart from its `source`:
- `noarch: python`: one pure-Python build for every platform
- `run_exports` pinned with `pin_subpackage(..., max_pin="x.x")`
- tests: package imports, `--help` for both commands, and `pip check`
- MIT license, shipped as `LICENSE`

It installs two commands: `stola` and `stola-compare-cross-orthology`. The cluster mode of the second needs matplotlib, pandas and toytree (3.x API), so those are run dependencies of the conda package. With pip they're the optional `cluster` extra instead.

## Version
The version lives in two places, which must match:
- `src/stola/__init__.py` (`__version__`, read by pyproject.toml)
- `recipe/meta.yaml` (`{% set version = ... %}`)

## Local build and install
Build in a dedicated environment holding conda-build, so other environments are untouched:

```
conda create -n conda-build -c conda-forge conda-build
conda run -n conda-build conda build recipe -c conda-forge -c bioconda --output-folder dist/conda
conda create -n stola -c ./dist/conda -c conda-forge -c bioconda stola
```

`conda build` also runs the recipe's tests on the built package. `dist/` is gitignored.

## Releasing to bioconda
Bioconda builds from a public, versioned source tarball, so these must be in place first:
- the GitHub repo is public
- the release is tagged, e.g. `v0.3.0`, matching the recipe version

First release:
1. Tag and push the release (`git tag v0.3.0 && git push origin v0.3.0`), and create a GitHub release from the tag.
2. Get the tarball's checksum:
   `curl -sL https://github.com/danielle-khost/stola/archive/refs/tags/v0.3.0.tar.gz | shasum -a 256`
3. Fork `bioconda/bioconda-recipes` and add `recipes/stola/meta.yaml`, copied from `recipe/meta.yaml` but with `source:` replaced by the commented-out `url:` and `sha256:` lines.
4. Open a pull request. Bioconda's CI lints, builds and tests the recipe, then a bioconda reviewer merges it.

Once merged, anyone can install with:
`conda create -n stola -c conda-forge -c bioconda stola`

Later releases: tag a new version; bioconda's bot normally opens the version-bump PR automatically. If a release changes dependencies, update `requirements` in the bioconda recipe as well as here.

Contributor guide: https://bioconda.github.io/contributor/index.html
