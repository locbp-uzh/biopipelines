# Vendored browser libraries

Third-party JavaScript committed here verbatim so that `renderers/pipeline_report.py` can inline it into a pipeline page. A page has to open from a `file://` URL on a laptop after being copied off a cluster, and page generation itself runs on compute nodes that often have no outbound network, so both rendering and viewing must work with no connectivity. That rules out fetch-and-cache and leaves a committed copy.

Nothing here is imported by Python. `renderers/pipeline_report.py` reads the file as bytes and embeds it once per page; if a file listed below is absent, the page falls back to the CDN URL (only when `render_page(..., allow_external=True)`) and then to a metadata table.

## 3Dmol-min.js

| | |
|---|---|
| library | [3Dmol.js](http://3dmol.org) |
| version | 2.5.2 |
| source URL | `https://cdn.jsdelivr.net/npm/3dmol@2.5.2/build/3Dmol-min.js` |
| size | 524,222 bytes |
| SHA-256 | `7b26bfd8170372ea78ea06df4a45e4a55fa5f538c2dfd716166fedf31ebfce9a` |
| licence | BSD-3-Clause — full text in `3Dmol-LICENSE.txt` |
| used by | `renderers/structures.py`, `renderers/grids.py` |

The version must stay in step with the `_LIB_3DMOL_URL` constant in `renderers/structures.py` and `renderers/grids.py`, which is what a CDN fallback would load.

### Updating it

```sh
curl -o renderers/vendor/3Dmol-min.js https://cdn.jsdelivr.net/npm/3dmol@<version>/build/3Dmol-min.js
curl -o renderers/vendor/3Dmol-LICENSE.txt https://cdn.jsdelivr.net/npm/3dmol@<version>/LICENSE
sha256sum renderers/vendor/3Dmol-min.js
```

Then update the table above and the `_LIB_3DMOL_URL` constants. The file is inlined into `<script>` verbatim, so a future version that contains the literal text `</script` would break page assembly; `renderers/pipeline_report.py` neutralizes that case, but check it anyway.
