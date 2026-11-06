# Atlas skeleton

Public repository: the one source of the atlas viewer every map site uses, and its demo at
https://atlas-skeleton.github.io. Names no private repository in files, pages, code comments or commit messages: it
says "the workbench", "the converter's repository", "a private repository".

- The viewer: `_includes/atlas-app.html` (list on the left, the map, its panel, search, downloads row). It never
  knows a map's source format.
- The renderers: `renderers/yrender/` (yEd GraphML, drawn as yEd draws it) and `renderers/converter/cdrender.py`
  (CellDesigner, drawn as CellDesigner draws it, with the converter's `mapconv/`). Where each was copied from:
  `VERSION`.
- The build: `python3 build.py` turns every file in `maps/` (`.graphml` or `.xml`) into the same output
  (`downloads/`, `atlas-data/`, `_maps/`); its header says what it writes. `maps/<id>.md` gives a map's title, group,
  order and text. The workflow `build.yml` runs it on every push that touches `maps/`, `renderers/` or `build.py`.
- Only what a map file holds is shown; nothing is fetched from databases (for now).
- Demo maps are maps already public (the Metabolism Regulation website's official files), never a collaborator's
  unpublished map.
- Every change goes to `main` and is checked on the live demo before a site takes the new version. A site gets the
  viewer as a copy with a VERSION file naming this repository and its commit.
- Hidden from search engines for now (`noindex` in `_config.yml`).
