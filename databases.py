# What the panel shows beyond a map file's identifiers: names, function, location, formula, cross-references and the
# enzymes of drawn sites from UniProt, HGNC and ChEBI, and citations from PubMed. Each identifier is fetched once and
# kept in cache/<source>.json (shared by all maps); "python3 build.py --refresh" fetches every entry again, keeping
# the old one when a fetch fails. An identifier that no longer resolves keeps its old entry, marked stale with the
# reason, and every build warns, naming the maps. cache/overrides.json, by hand, beats a fetched value and survives a
# refresh: {"uniprot": {"P47712": {"Location": "Cytoplasm"}}}, the same for "hgnc" (by HGNC ID) and "chebi"
# (CHEBI:...); "" hides a field. Copied from the Metabolism Regulation website's scripts/make_maps.py (VERSION).
import collections, json, os, re, sys, time, urllib.request, urllib.parse

ROOT = os.path.dirname(os.path.abspath(__file__))
FIELDS = ['PTM', 'Full name', 'Alternative names', 'Synonyms', 'Definition', 'Function', 'Location', 'Formula',
          'HGNC', 'HGNC symbol', 'Previous symbols', 'Alias symbols', 'Entrez', 'Ensembl', 'RefSeq', 'UniProt',
          'ChEBI', 'LIPID MAPS', 'KEGG', 'HMDB', 'InChI', 'InChIKey', 'Family', 'Members']
MODIFIED = {'P': 'Phosphorylated', 'Ac': 'Acetylated', 'Me': 'Methylated', 'Ub': 'Ubiquitinated',
            'Su': 'Succinylated', 'Glut': 'Glutarylated', 'Pa': 'Palmitoylated', 'My': 'Myristoylated'}
RESIDUE = {'S': 'Ser', 'T': 'Thr', 'Y': 'Tyr', 'K': 'Lys', 'C': 'Cys', 'G': 'Gly', 'R': 'Arg', 'H': 'His'}
REFRESH = '--refresh' in sys.argv
SOURCE = {'uniprot': 'https://rest.uniprot.org/uniprotkb/%s.json?fields=protein_name,cc_function,cc_subcellular_location,ft_mod_res',
          'hgnc': 'https://rest.genenames.org/fetch/hgnc_id/%s',
          'hgnc-symbol': 'https://rest.genenames.org/fetch/symbol/%s',
          'chebi': 'https://www.ebi.ac.uk/chebi/backend/api/public/compound/%s/',
          'pubmed': 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=pubmed&retmode=json&id=%s'}
CACHE = {n: os.path.join(ROOT, 'cache', n + '.json') for n in ('uniprot', 'hgnc', 'chebi', 'pubmed')}
cache = {n: json.load(open(f, encoding='utf-8')) if os.path.exists(f) else {} for n, f in CACHE.items()}
fetched = set()  # what this run has fetched (with --refresh, each entry once)
OVERRIDES = json.load(open(os.path.join(ROOT, 'cache', 'overrides.json'), encoding='utf-8')) if os.path.exists(os.path.join(ROOT, 'cache', 'overrides.json')) else {}
stale = collections.defaultdict(set)  # (source, identifier, reason) -> the maps that have it
current = ''  # the map being read, for the warnings
PREFIX = {'P': 'phospho', 'Ac': 'acetyl', 'Me': 'methyl', 'Su': 'succinyl', 'Glut': 'glutaryl'}


def fetch(url):
    req = urllib.request.Request(url, headers={'Accept': 'application/json', 'User-Agent': 'atlas-skeleton.github.io'})
    with urllib.request.urlopen(req, timeout=60) as r:
        return json.load(r)


def gone(src, key, d):
    """Why an identifier no longer resolves, from the source's answer ('' when it does)."""
    if src == 'uniprot':
        r = d.get('inactiveReason')
        if r:
            return r['inactiveReasonType'].lower() + (' into ' + ', '.join(r['mergeDemergeTo']) if r.get('mergeDemergeTo') else '')
        return 'now ' + d['primaryAccession'] if d['primaryAccession'] != key else ''
    if src == 'chebi':
        return 'now ' + d['chebi_accession'] if d['chebi_accession'] != key else ''
    if src == 'pubmed':
        return 'not found' if 'error' in d['result'].get(key, {'error': 1}) else ''
    docs = d['response']['docs']
    return 'not found' if not docs else docs[0]['status'].lower() if docs[0]['status'] != 'Approved' else ''


def cached(src, key, url, parse):
    if key and (key not in cache[src] or REFRESH and (src, key) not in fetched):
        fetched.add((src, key))
        why = ''
        try:
            d = fetch(url)
            why = gone(src, key, d)
            if not why:
                cache[src][key] = parse(d)
        except (OSError, ValueError, KeyError, IndexError) as e:
            if getattr(e, 'code', None) in (400, 404):  # the source has no such entry
                why = 'not found'
            else:
                print('%s %s: not fetched (%s)%s' % (src, key, e, ', the old entry kept' if key in cache[src] else ''))
        if why:
            cache[src][key] = dict(cache[src].get(key, {}), stale=why)
    got = cache[src].get(key, {})
    if got.get('stale'):
        stale[(src, key, got['stale'])].add(current)
    return dict(got, **OVERRIDES.get(src, {}).get(key, {}))


def first_sentences(text, n=2):
    # UniProt's function text without its evidence, in its first two sentences
    text = re.sub(r'\s*\((?:PubMed|By similarity|Ref\.)[^()]*\)', '', text)
    return ' '.join(re.split(r'(?<=\.)\s+(?=[A-Z0-9])', text)[:n]).strip()


def parse_uniprot(u):
    pd, com = u.get('proteinDescription', {}), u.get('comments', [])
    fn = [c['texts'][0]['value'] for c in com if c['commentType'] == 'FUNCTION' and c.get('texts')]
    return {'Full name': pd.get('recommendedName', {}).get('fullName', {}).get('value', ''),
            'Alternative names': '; '.join(a['fullName']['value'] for a in pd.get('alternativeNames', [])),
            'Function': first_sentences(fn[0]) if fn else '',
            # every location, the isoforms' too, in UniProt's order
            'Location': ', '.join(dict.fromkeys(l['location']['value'] for c in com if c['commentType'] == 'SUBCELLULAR LOCATION'
                                                for l in c.get('subcellularLocations', []) if 'location' in l)),
            # position -> [description, enzymes] for the sites UniProt gives an enzyme for
            'sites': {str(f['location']['start']['value']): [f['description'].split(';')[0], f['description'].split('; by ', 1)[1].split(';')[0]]
                      for f in u.get('features', []) if '; by ' in f.get('description', '')
                      and any(e.get('evidenceCode') == 'ECO:0000269' for e in f.get('evidences', []))}}


def parse_hgnc(h):
    d = h['response']['docs'][0]
    return {'HGNC': d['hgnc_id'], 'HGNC symbol': d['symbol'],
            'Previous symbols': ', '.join(d.get('prev_symbol', [])), 'Alias symbols': ', '.join(d.get('alias_symbol', [])),
            'Entrez': d.get('entrez_id', ''), 'Ensembl': d.get('ensembl_gene_id', ''),
            'RefSeq': ', '.join(d.get('refseq_accession', [])), 'uniprot': d.get('uniprot_ids', [])}


def parse_chebi(d):
    # names checked by ChEBI only (status C), not the brand names it takes from ChEMBL
    names, syn = sorted({x['ascii_name'] for l in (d.get('names') or {}).values() for x in l if x.get('status') == 'C'}), []
    seen = {d['name'].lower(), d['ascii_name'].lower()}
    for x in names:
        if x.lower() not in seen:
            seen.add(x.lower())
            syn.append(x)
    st, acc = d.get('default_structure') or {}, d.get('database_accessions') or {}
    xref = lambda pre: ', '.join(dict.fromkeys(a['accession_number'] for t, l in acc.items() if t != 'CAS'
                                               for a in l if (a.get('prefix') or '') == pre))
    return {'Full name': d['ascii_name'], 'Synonyms': '; '.join(syn),
            'Definition': re.sub(r'<[^>]+>', '', d.get('definition') or ''),
            'Formula': (d.get('chemical_data') or {}).get('formula') or '',
            'KEGG': xref('kegg.compound'), 'HMDB': xref('hmdb'),
            'InChI': st.get('standard_inchi') or '', 'InChIKey': st.get('standard_inchi_key') or ''}


# A process's references (user, 2026-10-27): the node property References, prefix:id separated by
# spaces (pubmed:32441256 doi:10.1016/...), as the converter reads and writes them. The panel lists each
# as a citation, first author, year, title and journal from PubMed (atlas-data/pubmed.json, fetched
# once per PMID as the other sources), with its PUBMED link; a DOI only as its link.
def parse_pubmed(d):
    r = next(v for k, v in d['result'].items() if k != 'uids')
    names = [a['name'] for a in r.get('authors', []) if a.get('authtype') == 'Author']
    return {'authors': (names[0] + (' et al.' if len(names) > 1 else '')) if names else '',
            'year': r.get('pubdate', '')[:4], 'title': r.get('title', ''), 'journal': r.get('source', '')}


def references(text):
    out = []
    for x in (text or '').split():
        db, _, i = x.partition(':') if ':' in x else ('pubmed', '', x)
        db = 'pubmed' if db.lower() in ('pubmed', 'pmid') else db.lower()
        if db == 'pubmed' and (i not in cache['pubmed'] or REFRESH):
            time.sleep(0.4)  # PubMed answers at most three requests a second
        out.append(dict(cached('pubmed', i, SOURCE['pubmed'] % i, parse_pubmed), db=db, id=i) if db == 'pubmed' else {'db': db, 'id': i})
    return out


def uniprot(acc):
    return cached('uniprot', acc, SOURCE['uniprot'] % acc, parse_uniprot)


def hgnc(hgnc_id='', symbol=''):
    if not hgnc_id and symbol:
        hgnc_id = next((k for k, v in cache['hgnc'].items() if v.get('HGNC symbol') == symbol), '')
        if not hgnc_id:  # a symbol not seen yet: found by the symbol, kept under its HGNC ID
            try:
                h = parse_hgnc(fetch(SOURCE['hgnc-symbol'] % urllib.parse.quote(symbol)))
                hgnc_id = h['HGNC']
                cache['hgnc'][hgnc_id] = h
                fetched.add(('hgnc', hgnc_id))
            except (OSError, ValueError, KeyError, IndexError) as e:
                if isinstance(e, IndexError):  # no approved symbol of that name
                    stale[('hgnc', symbol, 'not an approved symbol')].add(current)
                else:
                    print('hgnc %s: not fetched (%s)' % (symbol, e))
                return {}
    return cached('hgnc', hgnc_id, SOURCE['hgnc'] % hgnc_id, parse_hgnc)


def chebi(chebi_id):
    return cached('chebi', chebi_id, SOURCE['chebi'] % chebi_id.split(':')[-1], parse_chebi)


def ptm_text(ptm, acc=''):
    # by position in the protein (Gly2, Cys15, Cys26), not in the stored order around the node
    sites = sorted((s.strip().partition('@') for s in ptm.split('\n')),
                   key=lambda x: int(re.sub(r'\D', '', x[2]) or 1e9))
    known = uniprot(acc.split(',')[0].strip()).get('sites', {}) if acc else {}
    def by(sym, site):
        u = known.get(site[1:])
        return ', by ' + u[1] if u and u[0].lower().startswith(PREFIX.get(sym, '-')) else ''
    return '\n'.join(MODIFIED.get(sym, sym) + (' at ' + RESIDUE.get(site[:1], site[:1]) + site[1:] + by(sym, site) if site else '')
                     for sym, _, site in sites)


def record(data):
    """The panel's properties of a node or a family member, in FIELDS order: its identifiers
    and curated fields from the GraphML, the rest from the cache (the GraphML's own value where
    the cache has none)."""
    rec = dict(data)
    h = hgnc(data.get('HGNC', ''), data.get('HGNC symbol', ''))
    if data.get('HGNC symbol') and h.get('HGNC symbol', data['HGNC symbol']) != data['HGNC symbol']:
        stale[('hgnc', data['HGNC symbol'], 'now ' + h['HGNC symbol'])].add(current)
    acc = (data.get('UniProt') or '').split(',')[0].strip()
    if not acc and len(h.get('uniprot', [])) > 1:
        print('%s: HGNC gives several UniProt entries (%s), none taken' % (h['HGNC symbol'], ', '.join(h['uniprot'])))
    acc = acc or (h['uniprot'][0] if len(h.get('uniprot', [])) == 1 else '')
    if acc:
        rec['UniProt'] = data.get('UniProt') or acc
    for got in (h, uniprot(acc) if acc else {}, chebi(data['ChEBI']) if data.get('ChEBI') else {}):
        rec.update({k: v for k, v in got.items() if v and k in FIELDS})
    return [[f, ptm_text(rec[f], rec.get('UniProt', '')) if f == 'PTM' else rec[f]] for f in FIELDS if rec.get(f)]


def save():
    os.makedirs(os.path.join(ROOT, 'cache'), exist_ok=True)
    for n, f in CACHE.items():
        json.dump(cache[n], open(f, 'w', encoding='utf-8'), ensure_ascii=False, indent=0, sort_keys=True)
    for (src, key, why), maps in sorted(stale.items()):
        print('warning: %s %s %s (%s)' % (src, key, why, ', '.join(sorted(maps))))
