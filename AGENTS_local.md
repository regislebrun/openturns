# AGENTS_local.md

## Machine purpose

This machine has two roles:

1. **OpenTURNS development**: each branch gets its own clone on the SSD and a
   build directory in `$HOME/tmp/<branch-name>/`. Install target goes into
   `<build-dir>/install/`. PYTHONPATH is set per-invocation when testing:
   ```bash
   PYTHONPATH="$HOME/tmp/<branch>/install/lib64/python3.13/site-packages" \
     python3 python/test/t_DataEngine_std.py
   ```
   Mount point: `/run/media/maison/Extreme SSD/openturns-nomade` (note the
   space — always quote PYTHONPATH values).

2. **Uncertainty quantification (UQ)**: uses upstream master, installed into
   `$HOME/openturns-adm/build/install/`. The PYTHONPATH entry in `~/.bashrc`
   makes this the default when no override is set:
   ```bash
   export PYTHONPATH="/home/maison/openturns-adm/build/install/lib64/python3.13/site-packages":...
   ```
   Development overrides take precedence because `PYTHONPATH` entries earlier
   in the string win, or the user sets it explicitly before invoking python.

When working on a branch, the per-invocation PYTHONPATH should be set
*before* any python call so it takes priority over the bashrc default.

## Configure
```bash
cmake -S /run/media/maison/Extreme\ SSD/openturns-nomade \
      -B $HOME/tmp/<branch-name> \
      -DCMAKE_INSTALL_PREFIX=$HOME/tmp/<branch-name>/install
```
Only run when the build directory does not exist yet.

## Build
```bash
cmake --build $HOME/tmp/<branch-name> --target install --parallel 20
```

## Library
- C++ library is located into lib/src
- Mention new classes into the ChangeLog file

## Python bindings
- SWIG files entries must be created into python/src/*.i
- New classes must be attached to experimental_module.i

## Tests
- C++ tests in `lib/test/`, Python tests in `python/test/`
- If a test has a corresponding `.expout` file; then the test output diffs against it
- Testcases must be added in the Python folder preferably, new C++ tests are not necessary but existing ones must keep passing
- Build a specific C++test: `cmake --build --target t_Axial_std $(( $(nproc) / 2 ))`
- Run a specific C++ test: `ctest -R cppcheck_Axial_std -V`
- Run a specific Python test: `ctest -R pyinstallcheck_Axial_std -V`
- Tests use `ott.assert_almost_equal` for floating-point numeric assertions
- Tests use `with ott.assert_raises(...)` for exception checks

## Documentation
- Docstring API documentation is located in python/src/*_doc.i files with SWIG directives
- Use R"RAW(...)RAW" delimiters in docstrings when backslashes are needed (ie latex formulas)
- Use latex macros from python/doc/math_notations.sty wherever possible to uniformize notations
- New classes documentation must contain the preamble warning: "This class is experimental..."
- Mathematical details of a new or refactored class go before the parameter
  description of the constructor, immediately below the experimental warning.
- Sphinx documentation is located in python/doc
- For new classes API entries need to be added in python/doc/user_manual
- The "Notes" numpydoc section will appear after the methods list,
  so only use it for specific elements like ResourceMap keys, corner cases etc
  to avoid important elements to be separated from the main documentation body.
- Relevant ResourceMap keys are listed as a bullet list in the Notes section.

## Conventions
- Use `<Classname|Area>: <Short task overview>` as commit title template if possible
- Commit message can end with `Closes #<issue number>` if it must reference a GH issue
- Avoid unicode characters in source, commit messages
- Add Python testcases for numeric issues
- Python files are linted with flake8, see `utils/lint.sh`
- No hard-coded numerical values in the code: use either relevant constants in
  SpecFunc (mainly `SpecFunc::Precision` for numerical tests of closeness) or
  keys in ResourceMap
- Any new entry in ResourceMap must have its counterpart in `openturns.conf.in`

## Lint tools
- `flake8`: installed at `~/.local/bin/flake8` (pip --user)
- `doc8`: install via `pip install --user doc8`
- `spellintian`: Debian-only (part of `lintian`), not available on Mageia 10.
  That check in `utils/lint.sh` can only run in CI.

## Local constraints
- Build directory: `$HOME/tmp/<branch-name>/` (suffixed by branch name)
- Installation prefix: `<build-dir>/install/`
- Before SWIG compilation: delete the `python` subdirectory in the build dir
- Use 20 threads to compile the C++ library
- Use 20 threads to compile the SWIG interface
- For sphinx: add `-DUSE_SPHINX=ON` at configure time and use exactly 1 thread,
  otherwise it freezes
- For HMAT: add `-DUSE_HMAT=ON` at configure time when needed
- Wait for the computer to be idle (load below 0.7 during at least 3mn) before starting any benchmark. For non-benchmark situations (builds, tests), proceed without waiting.
