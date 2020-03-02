# Central todo file

## documentation

- add docstrings and type annotations for all functions (continuous process)
- deal with sphinx in order to make a better index
- check for pdf export option (especially for Gernot and Falko)

## code

- add async features to improve import speed for sql and reading files
- use the pathlib correctly and consequently, especially when accessing directories
- identify performance bottlenecks and fix them
- refactor SFG spectrum, remove unnecessary code
- move all classes intended for daily use in separate modules **work.py**
- apply all PEP suggestions of PyCharm, fix all warnings and errors
- add support for calculation of surface coverage etc (already present in office)

## architectural decisions
- implement the way spectorm handles spectra (metaclass)
- design a way to build up the tables from config files
- "unhardcode" the import process: the importer must be controlled by config files
- implement handling for units with pint package (if it is useful)


