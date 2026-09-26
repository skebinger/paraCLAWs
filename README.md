# paraCLAWs
parallel conservation law solver package

# Temporary readme
- linking requires the lz library: zlib and zlib-devel! (on fedora, search for equivalent packages on other OS)

## Building with local source overrides

The Makefile can build the solver from a base source tree without copying every
source file into the user's working directory. Set `BASE_DIR` inside the
copied `Makefile` to the directory containing the distributed `src/` and `mk/`
directories. The build discovers user `.f90` files anywhere under `USER_DIR`
automatically, so users do not need to reproduce the base tree or edit the
Makefile to register their sources.

The build constructs one complete source list from the base source list and
all `.f90` files under `USER_DIR`. It matches declared Fortran modules and
submodules: a local file that declares the same unit as a base source replaces
that base source in the final list, while other local files are added. Thus a
replacement may be kept flat in `USER_DIR` without reproducing the base
directory structure or editing `USER_SRCS`:

```text
USER_DIR/
  source_function.f90
```

If `source_function.f90` declares `SMOD_source_function`, it automatically
replaces the base source that declares that submodule. The user source must
implement the same Fortran interface expected by the base code. A source that
does not replace a base unit is still included as an additional source.

Generated objects, module files, and the executable are placed in
`USER_DIR/build`, keeping the base solver directory free of
case-specific build artifacts.
