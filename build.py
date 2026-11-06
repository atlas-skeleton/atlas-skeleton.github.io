# Builds the atlas from maps/: every map file drawn as its own editor draws it, with the data a click shows.
#   python3 build.py
# maps/<id>.graphml (yEd) is drawn by renderers/yrender, maps/<id>.xml (CellDesigner) by renderers/converter/cdrender.py.
# Both give the same output, so the viewer (_includes/atlas-app.html) never knows the source:
#   downloads/<id>.svg and a copy of the map file;
#   atlas-data/<id>.json: {"nodes": [...]}, each clickable node's box in SVG pixels, label, class and properties
#     ([name, value] pairs), a process with references also its reaction (in, out, by) and references;
#   atlas-data/search.json: the search index of all maps;
#   _maps/<id>.md: the map's entry (front matter for the list and the downloads row, then its text).
# maps/<id>.md gives the entry: front matter title, group, order (and topic), its body the text under the map.
# A node is clickable when it has annotation, a process when it has references (as on the Metabolism Regulation
# Atlas). The panel's names, function, location and cross-references come from UniProt, HGNC and ChEBI, and the
# citations from PubMed, fetched once per identifier and cached (databases.py); a node whose identifier field holds
# several IDs (a family drawn as one node in a CellDesigner file) shows only what the file holds.
# python3 build.py --refresh fetches every cached entry again. Needs lxml.
import collections, glob, json, os, re, shutil, sys
from urllib.parse import unquote
from lxml import etree

ROOT = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(ROOT, 'renderers', 'yrender'))
sys.path.insert(0, os.path.join(ROOT, 'renderers', 'converter'))
import yrender  # noqa: E402
import cdrender  # noqa: E402
from mapconv.celldesigner import plain_name  # noqa: E402
import databases as db  # noqa: E402

FIELDS = db.FIELDS
PROCESSES = ('Process', 'OmittedProcess', 'UncertainProcess', 'Association', 'Dissociation')



# ---------- yEd GraphML ----------
GN, YN = '{http://graphml.graphdrawing.org/xmlns}', '{http://www.yworks.com/xml/graphml}'


def realizer(n):
    gn = n.find(GN + 'data/' + YN + 'GenericNode')
    return gn if gn is not None else n.find(GN + 'data/' + YN + 'ProxyAutoBoundsNode/' + YN + 'Realizers/*')


def first_label(gn):
    labels = [] if gn is None else [l.text.strip() for l in gn.findall(YN + 'NodeLabel')
                                    if l.text and l.text.strip() and not l.get('iconData')]
    return ' '.join(labels[0].split()) if labels else ''


def graphml(src, out_svg):
    yrender.render(src, out_svg)
    tx, ty = (float(v) for v in re.search(r'translate\((-?[\d.]+),(-?[\d.]+)\)', open(out_svg).read()).groups())
    root = etree.parse(src).getroot()
    keys = {k.get('id'): k.get('attr.name') for k in root.findall(GN + 'key') if k.get('for') == 'node'}
    label = {n.get('id'): first_label(realizer(n)) for n in root.iter(GN + 'node')}
    arrow = lambda e: next((a.get('target', 'none') for a in e.iter(YN + 'Arrows')), 'none')
    arcs = [(e.get('source'), e.get('target'), arrow(e)) for e in root.iter(GN + 'edge')]
    nodes = []
    for n in root.iter(GN + 'node'):
        data = {keys.get(d.get('key')): (d.text or '').strip() for d in n.findall(GN + 'data') if d.get('key') in keys}
        gn = realizer(n)
        if gn is None:
            continue
        kind = (gn.get('configuration') or '').replace('com.yworks.sbgn.', '')
        geo = gn.find(YN + 'Geometry')
        x, y, w, h = (float(geo.get(a)) for a in ('x', 'y', 'width', 'height'))
        if kind in PROCESSES:
            if data.get('References'):
                nid = n.get('id')
                nodes.append(process(kind, [x + tx + w / 2, y + ty + h / 2],
                                     [label[s] for s, t, a in arcs if t == nid and a == 'none'],
                                     [label[t] for s, t, a in arcs if s == nid],
                                     [label[s] for s, t, a in arcs if t == nid and a != 'none'],
                                     db.references(data['References'])))
            continue
        own = {f: data[f] for f in FIELDS if data.get(f)}
        if not own:
            continue
        node = {'label': first_label(gn), 'class': kind, 'box': [round(x + tx, 1), round(y + ty, 1), round(w, 1), round(h, 1)],
                'props': record(own)}
        # a family's members, by their HGNC symbols in Members; a member's sites (PTM) from Member annotation
        old = {m['symbol']: dict(m['props']) for m in json.loads(data.get('Member annotation') or '[]')}
        symbols = [x.strip() for x in data.get('Members', '').split(';') if x.strip()] or list(old)
        members = [{'symbol': sym, 'props': record(dict({k: v for k, v in old.get(sym, {}).items() if k in FIELDS}, **{'HGNC symbol': sym}))}
                   for sym in symbols]
        if members:
            node['members'] = members
        nodes.append(node)
    return nodes


def record(own):
    # the panel's fields: fetched where every identifier is single, else only what the file holds
    if any(',' in own.get(k, '') for k in ('HGNC', 'HGNC symbol', 'UniProt', 'ChEBI')):
        return [[f, own[f]] for f in FIELDS if own.get(f)]
    return db.record(own)


def process(kind, at, ins, outs, by, refs):
    # a process glyph is small: a click a little outside it still finds it
    return {'label': '', 'class': kind, 'box': [round(at[0] - 12, 1), round(at[1] - 12, 1), 24, 24], 'props': [],
            'reaction': {'in': ins, 'out': outs, 'by': by}, 'references': refs}


# ---------- CellDesigner ----------
C = '{http://www.sbml.org/2001/ns/celldesigner}'
S = '{http://www.sbml.org/sbml/level2/version4}'
RDF = '{http://www.w3.org/1999/02/22-rdf-syntax-ns#}'
# MIRIAM namespaces to the viewer's fields (others are kept under their own name; the panel shows the ones it knows)
FIELD = {'uniprot': 'UniProt', 'hgnc': 'HGNC', 'hgnc.symbol': 'HGNC symbol', 'obo.chebi': 'ChEBI', 'chebi': 'ChEBI',
         'ncbigene': 'Entrez', 'entrez.gene': 'Entrez', 'ensembl': 'Ensembl', 'refseq': 'RefSeq', 'kegg.compound': 'KEGG',
         'hmdb': 'HMDB', 'lipidmaps': 'LIPID MAPS'}
CLASS = {'PROTEIN': 'Macromolecule', 'RECEPTOR': 'Macromolecule', 'ION_CHANNEL': 'Macromolecule',
         'TRUNCATED': 'Macromolecule', 'SIMPLE_MOLECULE': 'SimpleChemical', 'ION': 'SimpleChemical',
         'DRUG': 'SimpleChemical', 'GENE': 'NucleicAcidFeature', 'RNA': 'NucleicAcidFeature',
         'ANTISENSE_RNA': 'NucleicAcidFeature', 'COMPLEX': 'Complex', 'PHENOTYPE': 'Phenotype',
         'UNKNOWN': 'UnspecifiedEntity', 'DEGRADED': 'SourceAndSink'}


def miriam(e):
    out = []
    for li in e.iter(RDF + 'li'):
        m = re.match(r'urn:miriam:([^:]+):(.+)$', li.get(RDF + 'resource') or '')
        if m and (m.group(1), unquote(m.group(2)).strip()) not in out:
            out.append((m.group(1), unquote(m.group(2)).strip()))
    return out


def name(n):
    return ' '.join(plain_name(n or '').replace('_br_', ' ').split())


def celldesigner(src, out_svg):
    svg, hits = cdrender.render_hits(src)
    open(out_svg, 'w', encoding='utf-8').write(svg)
    model = etree.parse(src).getroot().find(S + 'model')
    species = {}
    for s in list(model.iter(S + 'species')) + list(model.iter(C + 'species')):
        cls = s.findtext('.//' + C + 'class') or ''
        props = collections.defaultdict(list)
        for ns, i in miriam(s):
            if ns in ('obo.chebi', 'chebi') and not i.upper().startswith('CHEBI:'):
                i = 'CHEBI:' + i
            if ns == 'hgnc' and not i.upper().startswith('HGNC:'):
                i = 'HGNC:' + i
            props[FIELD.get(ns, ns)].append(i)
        species[s.get('id')] = {'label': name(s.get('name')), 'class': CLASS.get(cls, cls),
                                'props': record({f: ', '.join(props[f]) for f in FIELDS if f in props}) if props else []}
    nodes = [{'label': species[a['species']]['label'], 'class': species[a['species']]['class'], 'box': a['box'],
              'props': species[a['species']]['props']}
             for a in hits['aliases'] if species.get(a['species'], {}).get('props')]
    alias_name = {a['id']: species.get(a['species'], {}).get('label', '') for a in hits['aliases']}
    sp_name = lambda e: species.get(e.get('species') or e.get('reactant') or e.get('product') or '', {}).get('label') \
        or alias_name.get(e.get('alias'), '')
    for rx in model.iter(S + 'reaction'):
        x = rx.find(S + 'annotation/' + C + 'extension')
        refs = db.references(' '.join(ns + ':' + i for ns, i in miriam(rx) if ns in ('pubmed', 'doi')))
        at = hits['reactions'].get(rx.get('id'), {}).get('at')
        if x is None or not refs or not at:
            continue
        ins = [sp_name(e) for e in x.iterfind(C + 'baseReactants/' + C + 'baseReactant')] + \
              [sp_name(e) for e in x.iterfind(C + 'listOfReactantLinks/' + C + 'reactantLink')]
        outs = [sp_name(e) for e in x.iterfind(C + 'baseProducts/' + C + 'baseProduct')] + \
               [sp_name(e) for e in x.iterfind(C + 'listOfProductLinks/' + C + 'productLink')]
        by = [' and '.join(alias_name.get(a, '') for a in (m.get('aliases') or '').split(','))
              for m in x.iterfind(C + 'listOfModification/' + C + 'modification')
              if not (m.get('type') or '').startswith('BOOLEAN_LOGIC_GATE')]
        nodes.append(process('Process', at, ins, outs, by, refs))
    return nodes


# ---------- the atlas ----------
def search_entries(mid, nodes):
    """The search index of one map: one entry per molecule (copies of a node
    are joined), with its node indices in atlas-data/<base>.json.
    a: label and symbols, b: other names, i: database IDs (exact match only)."""
    out = collections.OrderedDict()
    for idx, n in enumerate(nodes):
        if 'reaction' in n:  # processes are not searched
            continue
        p = dict(n['props'])
        members = n.get('members') or []
        if members:
            key = 'fam:' + ','.join(sorted(m['symbol'] for m in members))
        else:
            key = (p.get('UniProt') or p.get('ChEBI') or p.get('HGNC symbol') or '').strip()
            key = key or 'lab:%s:%s' % (n['class'], n['label'].lower())
        a, b, i = [n['label']], [], []
        for rec in [p] + [dict(m['props'], _sym=m['symbol']) for m in members]:
            for f in ('_sym', 'HGNC symbol'):
                if rec.get(f):
                    a.append(rec[f])
            for f in ('Full name', 'Family'):
                if rec.get(f):
                    b.append(rec[f])
            for f in ('Alternative names', 'Synonyms', 'Previous symbols', 'Alias symbols'):
                if rec.get(f):
                    # HGNC symbols are separated by commas, the names by semicolons
                    b += [x.strip() for x in rec[f].split(';' if f in ('Alternative names', 'Synonyms') else ',') if x.strip()]
            for f in ('UniProt', 'ChEBI', 'HGNC', 'KEGG', 'HMDB', 'Entrez'):
                if rec.get(f):
                    i += [x.strip() for x in rec[f].split(',') if x.strip()]
        e = out.get(key)
        if not e:
            e = out[key] = {'k': key, 'l': n['label'], 'c': n['class'], 'm': mid, 'n': [], 'a': [], 'b': [], 'i': []}
        e['n'].append(idx)
        for f, v in (('a', a), ('b', b), ('i', i)):
            for x in v:
                if x not in e[f]:
                    e[f].append(x)
    return list(out.values())


def entry(mid):
    # maps/<id>.md: front matter (simple "key: value" lines) and the text
    path = os.path.join(ROOT, 'maps', mid + '.md')
    meta, body = {}, ''
    if os.path.exists(path):
        _, head, body = open(path, encoding='utf-8').read().split('---\n', 2)
        meta = {k.strip(): v.strip().strip('"') for k, _, v in (l.partition(':') for l in head.splitlines() if ':' in l)}
    return meta, body.strip()


def main():
    for d in ('downloads', 'atlas-data', '_maps'):
        shutil.rmtree(os.path.join(ROOT, d), ignore_errors=True)
        os.makedirs(os.path.join(ROOT, d))
    search = []
    sources = sorted(glob.glob(os.path.join(ROOT, 'maps', '*.graphml')) + glob.glob(os.path.join(ROOT, 'maps', '*.xml')))
    for src in sources:
        mid, ext = os.path.splitext(os.path.basename(src))
        db.current = mid
        svg = os.path.join(ROOT, 'downloads', mid + '.svg')
        nodes = (graphml if ext == '.graphml' else celldesigner)(src, svg)
        shutil.copy(src, os.path.join(ROOT, 'downloads', mid + ext))
        json.dump({'nodes': nodes}, open(os.path.join(ROOT, 'atlas-data', mid + '.json'), 'w', encoding='utf-8'),
                  ensure_ascii=False, separators=(',', ':'))
        search += search_entries(mid, nodes)
        w, h = (int(float(v)) for v in re.search(r'<svg[^>]* width="([\d.]+)" height="([\d.]+)"', open(svg).read()).groups())
        meta, body = entry(mid)
        fm = {'title': meta.get('title', mid), 'topic': meta.get('topic', meta.get('title', mid)), 'tab': '', 'page': '',
              'group': meta.get('group', ''), 'order': int(meta.get('order', 0)), 'kind': 'map', 'draft': False, 'tag': '',
              'image': '/downloads/%s.svg' % mid, 'width': w, 'height': h, 'scale': 1,
              'annotation': '/atlas-data/%s.json' % mid, 'image_error': '',
              'downloads': [{'label': 'GraphML' if ext == '.graphml' else 'CellDesigner', 'url': '/downloads/%s%s' % (mid, ext)},
                            {'label': 'SVG', 'url': '/downloads/%s.svg' % mid}]}
        with open(os.path.join(ROOT, '_maps', mid + '.md'), 'w', encoding='utf-8') as f:
            f.write('---\n' + ''.join('%s: %s\n' % (k, json.dumps(v, ensure_ascii=False)) for k, v in fm.items()) + '---\n' + body + '\n')
        print('%s: %s, %d clickable' % (mid, 'yEd' if ext == '.graphml' else 'CellDesigner', len(nodes)))
    db.save()
    json.dump(search, open(os.path.join(ROOT, 'atlas-data', 'search.json'), 'w', encoding='utf-8'),
              ensure_ascii=False, separators=(',', ':'))


if __name__ == '__main__':
    main()
