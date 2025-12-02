# qp-planck Documentation

This directory contains the Sphinx documentation source for **qp-planck**.

## Building Locally

1. Install dependencies (from project root):
   ```bash
   # Using uv (recommended)
   uv sync --all-extras
   
   # Or using pip
   pip install -e ".[dev]"
   ```

2. Build HTML documentation:
   ```bash
   cd docs
   uv run sphinx-build -b html . _build/html
   # Or: make html
   ```

3. View the documentation:
   ```bash
   open _build/html/index.html
   ```
   (Or use your browser to navigate to `_build/html/index.html`)

## Building on Read the Docs

The documentation is automatically built on Read the Docs when changes are
pushed to the repository. Configuration is in `../.readthedocs.yaml`.

## Documentation Structure

```
docs/
├── index.rst               # Main landing page
├── installation.rst        # Installation guide
├── quickstart.rst         # Quick start tutorial
├── physics_background.rst # Theory and formalism
├── configuration.rst      # Python API configuration
├── yaml_config.rst        # YAML configuration guide
├── examples.rst           # Example workflows
├── contributing.rst       # Contribution guidelines
├── changelog.rst          # Version history
├── references.rst         # Bibliography and links
├── api/                   # API reference
│   ├── utilities.rst
│   ├── hmap2mat.rst
│   ├── mat2fits.rst
│   └── pipeline.rst
├── conf.py               # Sphinx configuration
└── Makefile              # Build automation (Unix)
```

## Writing Documentation

### RST Syntax

We use reStructuredText (.rst) format. Key syntax:

**Headers:**
```rst
Title
=====

Section
-------

Subsection
^^^^^^^^^^
```

**Code blocks:**
```rst
.. code-block:: python

   from qp_planck import load_RIMO
   RIMO = load_RIMO("path/to/file.fits")
```

**Cross-references:**
```rst
See :doc:`installation` for setup instructions.
See :func:`qp_planck.load_RIMO` for API details.
```

**Math:**
```rst
.. math::

   C_\ell^{\text{obs}} = C_\ell^{\text{sky}} \cdot W_\ell
```

### Adding New Pages

1. Create `new_page.rst` in `docs/`
2. Add to table of contents in `index.rst`:
   ```rst
   .. toctree::
      :maxdepth: 2
      
      new_page
   ```
3. Build and verify

### Docstrings

API documentation is automatically generated from docstrings in the source code.
Use Google-style docstrings:

```python
def my_function(param: str) -> int:
    """Brief one-line description.
    
    Longer description with details.
    
    Parameters
    ----------
    param : str
        Description of parameter.
    
    Returns
    -------
    int
        Description of return value.
    
    Examples
    --------
    >>> my_function("test")
    42
    """
    pass
```

## Cleaning Build Artifacts

```bash
make clean
```

## Troubleshooting

**Import errors during build:**
- Ensure all dependencies are installed: `uv sync --all-extras`
- Check that the package is installed in development mode

**Missing modules in API docs:**
- Run `make clean` then `make html` to rebuild from scratch
- Check that modules are listed in `conf.py` extensions

**Math not rendering:**
- Ensure MathJax is configured in `conf.py`
- Use proper LaTeX syntax in math blocks

## Resources

- [Sphinx documentation](https://www.sphinx-doc.org/)
- [RST primer](https://www.sphinx-doc.org/en/master/usage/restructuredtext/basics.html)
- [Read the Docs](https://docs.readthedocs.io/)
- [Sphinx RTD theme](https://sphinx-rtd-theme.readthedocs.io/)
