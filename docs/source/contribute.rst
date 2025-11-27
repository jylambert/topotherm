.. _contributing:

Contributing to Topotherm
=========================

Thank you for considering a contribution! We welcome code, documentation,
examples, bug reports, reviews, and ideas.

If you're new to contributing, a great first step is to open a small issue
(or draft PR) describing what you plan to do. Maintainers will help you
get it into shape.

What can I contribute?
----------------------

- **Bug reports & fixes**: minimal reproducer + expected vs. actual behaviour
- **Features**: small, incremental features are easiest to review
- **Documentation**: tutorials, examples, clarifications, typos
- **Tests**: add tests that cover a bug or new behaviour
- **Review**: test and comment on open pull requests

Where to get help
-----------------

- **Issues**: use GitHub Issues for bugs, requests, and questions.
- **Discussions**: design/architecture questions fit well in a “discussion”
  issue before implementation.
- **Security**: report privately (do not open a public issue).

Set up a development environment
--------------------------------

We recommend conda/mamba, but any virtual env works.
# TODO check this
.. code-block:: bash

   # clone your fork (or the main repo if you have push rights)
   git clone https://github.com/<your-user>/topotherm.git
   cd topotherm

   # create & activate an environment
   mamba env create -f environment.yml -n topotherm-dev  # or: conda env create ...
   mamba activate topotherm-dev

   # install in editable mode with dev & docs extras
   pip install -e ".[dev,docs]"

   # (optional) install pre-commit hooks to auto-format/check on every commit
   pip install pre-commit
   pre-commit install

Coding style & conventions
--------------------------

- Python version: 3.10+ (match project config/CI)
- Type hints: required in public APIs. Keep them descriptive.
- Docstrings: NumPy-style (Sphinx is configured with ``napoleon``).
- Imports: standard library → third-party → local; avoid unused imports.
- Formatting/linting (recommended):
  - `black` for code formatting
  - `isort` for imports
  - `ruff` (or `flake8`) for linting

If you enabled ``pre-commit``, these run automatically. Otherwise:
# TODO check this
.. code-block:: bash

   ruff check .           # or: flake8
   isort .
   black .

Run the test suite
------------------

We use ``pytest``. Always run tests locally before opening a PR.

.. code-block:: bash

   # install test deps (already included by .[dev])
   pip install -e ".[dev]"

   # run tests
   pytest

   # (optional) coverage report
   pytest --maxfail=1 --disable-warnings -q --cov=topotherm --cov-report=term-missing
   # TODO CHECK THIS

Add or update documentation
---------------------------

- User-facing changes should update the docs.
- API docs are built automatically from docstrings (AutoAPI).
- The docs homepage includes ``README.md``; section pages
  (Install, Usage, …) include slices of the README.

Build docs locally:

.. code-block:: bash

   pip install -e ".[docs]"
   sphinx-build -b html docs/source docs/_build/html
   # open docs/_build/html/index.html

Git workflow
------------

1. Fork the repository and create a topic branch:
   ``git checkout -b feature/short-description``.
2. Keep changes focused; small PRs are faster to review.
3. Sync with main before you open the PR:
   ``git fetch upstream && git rebase upstream/main`` (or merge).
4. Write good commit messages (explain why, not just what).
   Conventional Commits (e.g. ``fix:``, ``feat:``, ``docs:``) are welcome but not required.
5. Push and open a Draft PR early if you want feedback.

Pull request checklist
----------------------

- [ ] Tests added/updated and pass locally (``pytest``)
- [ ] Docs updated if user-facing behaviour changed
- [ ] Type hints and NumPy-style docstrings for public APIs
- [ ] ``pre-commit`` (or ``ruff/isort/black``) ran cleanly
- [ ] brief summary in PR description

Reporting bugs
--------------

Please include:

- Topotherm version and Python version
- OS, solver, and any relevant dependencies
- Minimal reproducible example (code + tiny data)
- Expected vs. actual behaviour and full traceback

License & contributor certificate
---------------------------------

By contributing, you agree your work will be released under the project license
(MIT). Make sure you have the right to contribute the code and that third-party
snippets are compatible with the license.

MIT license: https://en.wikipedia.org/wiki/MIT_License, see `LICENSE` file.

Credits & inspiration
---------------------

Parts of this guide were inspired by high-quality open-source projects and their
contribution guides, such as PyPSA's “Contributing” documentation. See their docs
for further examples and best practices.
