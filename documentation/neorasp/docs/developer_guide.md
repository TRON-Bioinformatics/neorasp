# Developer Guide

## Software environment

The project uses [Pixi](https://pixi.prefix.dev/latest/) to manage all
development and runtime environments. After installing Pixi, enter the project
shell with:

```
pixi shell
```

This provides Snakemake and all dependencies needed to run the pipeline. For
deployments to systems without internet access or without Pixi, see
[Pixi pack](https://pixi.prefix.dev/latest/deployment/pixi_pack/).

The available Pixi environments and tasks are defined in
[`pixi.toml`](https://github.com/TRON-Bioinformatics/neorasp/blob/dev/pixi.toml).
The most relevant tasks are:

| Task                     | Description                                                          |
| ------------------------ | -------------------------------------------------------------------- |
| `pixi run test`          | Run the CI test suite (see [Tests](#tests)).                         |
| `pixi run test-sytax`    | Run the syntax test suite (see [Tests](#tests)).                     |
| `pixi run test-local`    | Run the test suite for HPC (see [Tests](#tests)).                    |
| `pixi run lint`          | Run all linters (Snakemake, Python, R, Markdown, shell, YAML, TOML). |
| `pixi run style`         | Auto-format all files (same scope as `lint`).                        |
| `pixi run lint-workflow` | Snakemake `--lint` of `workflow/Snakefile`.                          |
| `pixi run build-docs`    | Build the MkDocs site under `documentation/neorasp/`.                |

## Tests

CI currently includes a dry-run and a full integration test on a minimal samples. These verify workflow syntax, rule wiring and runtime of integrated tools.

Run the tests locally with:

```
pixi run test-sytax
pixi run test-local
```

## Code styling

Run `pixi run style` to auto-format code and docs. `pixi run lint` performs the
same checks without modifying files and exits non-zero on failure; it is used in
CI. Sub-tasks (`style-python`, `lint-snakemake`, ...) are available for
individual file types — see
[`pixi.toml`](https://github.com/TRON-Bioinformatics/neorasp/blob/dev/pixi.toml).

## Release

Before creating a new release:

- Bump the version in
  [`pixi.toml`](https://github.com/TRON-Bioinformatics/neorasp/blob/dev/pixi.toml).

## Contribute

[CONTRIBUTING.md](https://github.com/TRON-Bioinformatics/neorasp/blob/dev/CONTRIBUTING.md).
