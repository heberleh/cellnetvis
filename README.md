# Citation
Heberle, H., Carazzolle, M. F., Telles, G. P., Meirelles, G. V., & Minghim, R.  (2017). CellNetVis : a web tool for visualization of biological networks using force-directed layout constrained by cellular components. BMC Bioinformatics, 18(10), 395. 

https://doi.org/10.1186/s12859-017-1787-5

## What it is

CellNetVis is a web tool that draws biomolecular networks (XGMML files) inside a diagram of a cell. Each node is placed in
its cellular component (extracellular space, plasma membrane, cytoplasm, nucleus, mitochondrion, ...) by a force-directed
layout that is constrained by the shape of the component. Nodes can be searched, filtered, coloured by any attribute,
labelled, dragged, bundled, summarised in a distribution chart and exported as SVG/PNG. It is a static site: everything
runs in the browser and **your network never leaves your computer**.

This is **CellNetVis 1 (2017)**, kept on its original stack (D3 v3, jQuery 2.1, Bootstrap 3, nvd3, jsnetworkx,
d3.ForceBundle, FileSaver, cytoscape.js 2.2 / CerebralWeb). A new version is in development.

The hosted copy of the original publication (`http://bioinfo03.ibi.unicamp.br/lnbio/IIS2/cellnetvis/`) is still reachable
(2026-10-20, owner check), but its example-network links no longer load and it is served over plain HTTP.
A bundled example network (RIG-I-like receptor signaling pathway, from InnateDB) is available through the
**Load example** button, or with `index.html?file=examples/rig-i-like-receptor-innatedb.xgmml`.

![Images](https://media.springernature.com/full/springer-static/image/art%3A10.1186%2Fs12859-017-1787-5/MediaObjects/12859_2017_1787_Fig4_HTML.gif)

## Run it locally

Everything the site needs is in `web/`. Any static file server works; the repository has a tiny one without dependencies:

```bash
npm run dev          # http://localhost:8080  (PORT=3000 npm run dev for another port)
```

The server also applies `web/_headers`, so the Content-Security-Policy that is deployed is the one you test against.

Tests (Playwright, Chromium; the server is started by the test run):

```bash
npm install
npx playwright install chromium   # once, if you do not have a browser yet
npm test
# CHROMIUM_PATH=/path/to/chrome npm test   # use an already installed Chromium
```

The suite checks that no page requests another origin, that the console stays free of errors and warnings, the
cell diagram and the example network (107 nodes, 554 links after the 78 self-edges are removed), the force layout
(running, Stop, Start), network measures, search/filter/clear, labels, the donut chart, SVG/PNG downloads, the
attribute table, the CerebralWeb layout, the responsive layout (screenshots are written to `test-results/screenshots/`) and the notices
that replace `alert()`.

The `?file=` parameter loads a network from an address: `index.html?file=<address>`. It accepts a path on this site or an
`http(s)` address of a server that allows cross-origin requests (CORS); other schemes are refused.

## Deploy to Cloudflare Pages

It is a static site, there is nothing to build:

| Setting | Value |
|---|---|
| Framework preset | None |
| Build command | *(empty)* |
| Build output directory | `web` |

`web/_headers` is picked up automatically by Cloudflare Pages: `X-Content-Type-Options`, `Referrer-Policy`,
`Permissions-Policy` and a Content-Security-Policy that only allows the site's own files
(`'unsafe-inline'` for scripts and styles because the pages have inline blocks; `connect-src 'self' https:` is only
for `?file=` and the InnateDB localization call). Any other static host works the same way: publish the `web/` folder.

## What changed in 2026

Privacy and self-containment
- Removed Google Analytics (`gtag.js`, and the commented older snippet) from every page.
- Bootstrap 3.3.7 CSS, theme and JS (checked byte for byte against the official release) now come from `web/bootstrap/`,
  plus the Glyphicons fonts in `web/fonts/` (they were missing: the panel arrows were blank). The IE8 shim/Respond CDN
  scripts and the Google Fonts links were removed. The two YouTube iframes of the Help page are plain links, so nothing is
  requested from another server until you click.
- Added `web/_headers` (CSP and friends), `package.json`, `scripts/serve.mjs`, `.gitignore`.
- Removed unused or stale files: `WEB-INF/` (Java App Engine), `untitled.css`, `minify.js`, `all_scripts.min.js`
  (stale copy), `circular_cell_structure - Copia.svg`, `.vscode/`, the nvd3 source tree (tests and examples; its licence is
  kept in `web/nvd3/LICENSE.md`), the d3.ForceBundle example and images, unused jQuery copies and the unminified D3.

Bugs fixed in `all_scripts.js`
- `convertToCellNetVisCC`: `value.indexOf("cell_surface">=0)` compared the string with 0, so every unmatched value became
  "extracellular" (e.g. "mitochondrial inner membrane"). It is now `value.indexOf("cell_surface") >= 0`.
- `String.prototype.capitalize` was defined twice; the later one (every word) is kept.
- `integrate()` had two `error:` handlers in one `$.ajax` call; now one, with a clear message.
- Choosing a file looped over the file list and called `alert` on the second; only the first file is used.
- The InnateDB localization call used `http://` and `async: false`; it is now `https://`, asynchronous, and when the
  service is down (its status in 2026 is unknown) a notice explains that the network needs a "Selected CC"/"Localization"
  attribute instead of leaving the page half-loaded.
- Every `alert()` became a dismissible in-page notice (`role="status"`) at the top of the content column.
- `?file=` is decoded with `decodeURIComponent` (not `unescape`), may contain its own `=`/`?`, and must be http(s) or a
  path on this site.
- Attribute names and values from a file were written into the page as HTML (script injection); they are escaped now.
- Computing the network measures twice duplicated the options; loading a second network left the old options in the
  label menu; the "Submit CC" handler piled up once per loaded network; edges pointing to missing nodes broke the force
  layout; the `Selected CC` select was missing its closing tag; node counts were compared to numbers instead of
  lengths in `force_only_selected`; a failed load left the "loading" overlay on forever.
- The example link pointed to the dead IIS host: the bundled example is now `web/examples/rig-i-like-receptor-innatedb.xgmml`
  with a **Load example** button; the Help page links to it and marks the other examples (MAPK, HPA, yeast) as no longer
  available.

Layout
- The diagram is a responsive SVG (`viewBox="0 0 850 850"`, `width: 100%`, `max-width: 850px`): the force layout still
  works in 850 x 850 units, the browser only scales it. The columns stack on small screens (the old
  `min-width: 1704px` on `html` is gone) and the export keeps its 850 x 850 size.
- A quiet banner under the navbar of every page: "CellNetVis 1 (2017). A new version is in development ..." (dismissible for
  the visit, nothing is stored).

## Licence

GNU General Public License v3.0 or later, see [LICENSE](LICENSE). Third-party libraries keep their own licences
(D3, jQuery, Bootstrap, nvd3, jsnetworkx, d3.ForceBundle, FileSaver, cytoscape.js, bootstrap-slider, jQuery BlockUI).

## Citation

Please cite:

Heberle, H., Carazzolle, M. F., Telles, G. P., Meirelles, G. V., & Minghim, R. (2017). CellNetVis : a web tool for
visualization of biological networks using force-directed layout constrained by cellular components.
*BMC Bioinformatics*, 18(10), 395. https://doi.org/10.1186/s12859-017-1787-5
