Contributing
============

We welcome contributions to **qp-planck**! This page explains how to get
involved and contribute to the project.

Ways to Contribute
------------------

* **Report bugs** — File issues on GitHub
* **Request features** — Suggest new functionality or improvements
* **Submit code** — Fix bugs, add features, improve documentation
* **Improve documentation** — Clarify explanations, add examples
* **Share use cases** — Tell us how you're using qp-planck

Getting Started
---------------

Fork and Clone
^^^^^^^^^^^^^^

1. Fork the repository on GitHub:
   `https://github.com/paganol/qp_planck <https://github.com/paganol/qp_planck>`_

2. Clone your fork locally:

   .. code-block:: bash
   
      git clone https://github.com/YOUR-USERNAME/qp_planck.git
      cd qp_planck

3. Add the upstream repository:

   .. code-block:: bash
   
      git remote add upstream https://github.com/paganol/qp_planck.git

Development Installation
^^^^^^^^^^^^^^^^^^^^^^^^^

Install in editable mode with development dependencies:

.. code-block:: bash

   pip install -e ".[dev]"

This installs:

* Core dependencies (numpy, scipy, astropy, healpy, matplotlib)
* Development tools (pytest, ruff, pre-commit)
* Documentation tools (sphinx, sphinx-rtd-theme)

Set Up Pre-Commit Hooks
^^^^^^^^^^^^^^^^^^^^^^^^

We use pre-commit hooks to maintain code quality:

.. code-block:: bash

   pre-commit install

This will automatically run linting and formatting checks before each commit.

Development Workflow
--------------------

Create a Branch
^^^^^^^^^^^^^^^

.. code-block:: bash

   git checkout -b feature/my-new-feature

Use descriptive branch names:

* ``feature/xyz`` — New features
* ``bugfix/xyz`` — Bug fixes
* ``docs/xyz`` — Documentation updates
* ``refactor/xyz`` — Code refactoring

Make Changes
^^^^^^^^^^^^

Edit the code following our style guidelines (see below).

Run Tests
^^^^^^^^^

Before committing, run the test suite:

.. code-block:: bash

   pytest

Check code coverage:

.. code-block:: bash

   pytest --cov=qp_planck --cov-report=html

View coverage report at ``htmlcov/index.html``.

Code Style
^^^^^^^^^^

We use **ruff** for linting and formatting:

.. code-block:: bash

   # Check for issues
   ruff check .
   
   # Auto-fix issues
   ruff check --fix .
   
   # Format code
   ruff format .

Our style follows:

* PEP 8 with 90-character line length
* Google-style docstrings
* Type hints where appropriate
* Import sorting with ``isort`` rules

Commit Changes
^^^^^^^^^^^^^^

Write clear, descriptive commit messages:

.. code-block:: bash

   git add .
   git commit -m "Add feature: brief description
   
   Longer explanation of what changed and why.
   
   Fixes #123"

Reference related issues with ``Fixes #123`` or ``Closes #123``.

Push and Create Pull Request
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: bash

   git push origin feature/my-new-feature

Then open a pull request on GitHub:

1. Navigate to your fork
2. Click "New Pull Request"
3. Select your branch
4. Fill in the PR template
5. Submit for review

Pull Request Guidelines
------------------------

**Before Submitting:**

* [ ] Tests pass (``pytest``)
* [ ] Code follows style guidelines (``ruff check``)
* [ ] Documentation is updated if needed
* [ ] CHANGELOG.md is updated (for notable changes)
* [ ] Commit messages are clear and descriptive

**PR Description Should Include:**

* Summary of changes
* Motivation / context
* Related issues (``Fixes #123``)
* Testing performed
* Any breaking changes

**Review Process:**

* Maintainers will review your PR
* Address any feedback
* Once approved, your PR will be merged

Code Guidelines
---------------

Docstrings
^^^^^^^^^^

Use Google-style docstrings for all public functions and classes:

.. code-block:: python

   def my_function(param1: int, param2: str) -> bool:
       """Brief one-line description.
       
       Longer description explaining the function's purpose,
       behavior, and any important details.
       
       Parameters
       ----------
       param1 : int
           Description of param1.
       param2 : str
           Description of param2.
       
       Returns
       -------
       bool
           Description of return value.
       
       Raises
       ------
       ValueError
           When invalid input is provided.
       
       Examples
       --------
       >>> my_function(42, "test")
       True
       """
       # Implementation
       pass

Type Hints
^^^^^^^^^^

Add type hints to function signatures:

.. code-block:: python

   from typing import List, Optional, Tuple
   
   def process_detectors(
       detectors: List[str],
       rimo_path: str,
       threshold: Optional[float] = None
   ) -> Tuple[np.ndarray, dict]:
       """..."""
       pass

Error Handling
^^^^^^^^^^^^^^

Raise informative exceptions:

.. code-block:: python

   def load_data(filepath: str) -> np.ndarray:
       """Load data from file."""
       if not os.path.exists(filepath):
           raise FileNotFoundError(
               f"Data file not found: {filepath}"
           )
       
       try:
           data = np.load(filepath)
       except Exception as e:
           raise ValueError(
               f"Failed to load {filepath}: {e}"
           ) from e
       
       return data

Testing
-------

Writing Tests
^^^^^^^^^^^^^

Tests go in the ``tests/`` directory. Use pytest:

.. code-block:: python

   # tests/test_utilities.py
   
   import pytest
   from qp_planck import list_planck
   
   def test_list_planck_143ghz():
       """Test listing 143 GHz detectors."""
       dets = list_planck("143GHz")
       
       assert len(dets) == 11  # 4 PSBs + 3 SWBs
       assert "143-1a" in dets
       assert "143-5" in dets
   
   def test_list_planck_invalid():
       """Test invalid detector set raises error."""
       result = list_planck("INVALID")
       assert result == -1

Run specific tests:

.. code-block:: bash

   pytest tests/test_utilities.py
   pytest tests/test_utilities.py::test_list_planck_143ghz

Test Coverage Goals
^^^^^^^^^^^^^^^^^^^

* Aim for >80% code coverage
* Focus on testing public APIs
* Include edge cases and error conditions
* Test MPI code paths when possible

Documentation
-------------

Building Documentation
^^^^^^^^^^^^^^^^^^^^^^

Build HTML documentation locally:

.. code-block:: bash

   cd docs
   make html

View at ``docs/_build/html/index.html``.

Clean build files:

.. code-block:: bash

   make clean

Adding Documentation
^^^^^^^^^^^^^^^^^^^^

* Update relevant ``.rst`` files in ``docs/``
* Add examples to docstrings
* Include plots/figures if helpful
* Cross-reference related docs with ``:doc:`path```

Documentation lives in:

* ``docs/`` — RST documentation source
* Docstrings in ``qp_planck/*.py`` — API reference
* ``README.md`` — Project overview
* ``CHANGELOG.md`` — Version history

Reporting Issues
----------------

Bug Reports
^^^^^^^^^^^

When reporting bugs, include:

1. **Description** — What went wrong?
2. **Steps to reproduce** — Minimal example to trigger the bug
3. **Expected behavior** — What should have happened?
4. **Actual behavior** — What actually happened?
5. **Environment** — Python version, OS, package versions
6. **Traceback** — Full error message/stack trace

Example issue:

.. code-block:: text

   **Bug**: `list_planck()` fails for LFI detectors
   
   **Steps to reproduce**:
   ```python
   from qp_planck import list_planck
   dets = list_planck("LFI")
   ```
   
   **Expected**: List of all LFI detectors
   
   **Actual**: `KeyError: 'LFI27M'`
   
   **Environment**: 
   - Python 3.11.4
   - qp-planck 0.1.0
   - Ubuntu 22.04
   
   **Traceback**:
   ```
   Traceback (most recent call last):
     ...
   KeyError: 'LFI27M'
   ```

Feature Requests
^^^^^^^^^^^^^^^^

For feature requests, describe:

1. **Motivation** — Why is this feature needed?
2. **Use case** — How would you use it?
3. **Proposed solution** — Ideas for implementation
4. **Alternatives** — Other approaches considered

Release Process
---------------

Maintainers follow this process for releases:

1. Update version in ``pyproject.toml``
2. Update ``CHANGELOG.md``
3. Create release commit: ``git commit -m "Release vX.Y.Z"``
4. Tag release: ``git tag vX.Y.Z``
5. Push: ``git push && git push --tags``
6. GitHub Actions builds and publishes to PyPI (when configured)

Communication
-------------

* **GitHub Issues** — Bug reports, feature requests
* **Pull Requests** — Code contributions, discussion
* **Discussions** — General questions, ideas (if enabled)

Code of Conduct
---------------

We follow the Python Community Code of Conduct:
`https://www.python.org/psf/conduct/ <https://www.python.org/psf/conduct/>`_

In summary:

* Be respectful and inclusive
* Assume good faith
* Give and accept constructive feedback
* Focus on what's best for the community

License
-------

By contributing, you agree that your contributions will be licensed under
the MIT License. See ``LICENSE`` file for details.

Questions?
----------

If you're unsure about anything:

* Open an issue to ask
* Look at existing PRs for examples
* Reach out to maintainers

We appreciate all contributions, big and small!

Thank You!
----------

Your contributions help make qp-planck better for everyone. We're grateful
for your time and effort! 🚀
