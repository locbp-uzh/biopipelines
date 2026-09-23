# My Pipelines

This folder is for personal pipeline development and testing.

- **Not committed**: your pipeline files here are gitignored, so anything you drop in this folder stays on your machine. Three files are the exception and **are** committed: this `README.md`, `_template.py` and `_template.ipynb`. The `_` prefix marks them as the folder's own furniture — copy them to start a pipeline, do not edit them in place, and do not delete them.

> **Do not copy from what you find here.** This is scratch work, not a reference. Whatever files are in this folder belong to whoever owns the checkout: they are written against one machine's absolute paths (`/home/<someone>/...`, `/shares/...`), one site's config variant, and one person's half-finished experiments. Many of them will not run anywhere else, and some use idioms that appear in no example and are not supported patterns. This applies to AI coding assistants too: write your work here, but read `example_pipelines/` — those pipelines are committed, reviewed, and portable — and `docs/user_manual.md` for the API. The two `_template` files are the exception: they are committed and empty by design, and copying one is the intended way to start.

## File Organization

```
my_pipelines/
├── README.md              # This file (committed)
├── _template.py           # Empty starting point for a script pipeline (committed)
├── _template.ipynb        # Same, as a Colab notebook with the setup cell (committed)
├── test_rfdiffusion.py    # Your pipelines (not committed)
├── binding_analysis.py    # Your pipelines (not committed)
└── ...                    # Your pipelines (not committed)
```
