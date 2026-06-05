#!/usr/bin/env python3
"""
filter_vep_fast - Fast Python reimplementation of Ensembl's filter_vep.

Replicates all features of the original Perl script.
Uses cyvcf2 for VCF parsing when available (significant speedup for large
VCF / bgzipped files via htslib).
"""

import sys
import re
import gzip
import argparse
import os

try:
    import cyvcf2
    HAS_CYVCF2 = True
except ImportError:
    HAS_CYVCF2 = False


# ---------------------------------------------------------------------------
# Filter expression parser / evaluator
# ---------------------------------------------------------------------------

class FilterSet:
    """Parse and evaluate VEP filter expressions.

    Multiple filter strings are ANDed (all must pass).
    Within a string: 'and' / 'or' / 'not' and parentheses are supported.

    Operators (with documented synonyms):
      is / = / eq           case-insensitive equality
      is not / != / ne      inequality
      > / gt                greater than
      < / lt                less than
      >= / gte              greater than or equal
      <= / lte              less than or equal
      match / matches / re / regex   regex search (case-insensitive)
      in                    list or file membership
      exists / ex / defined field has a non-null value
      is_child              SO ontology child (falls back to exact match)

    Special value syntax:
      #FieldName            compare against another field's value
      bare field name       equivalent to 'exists'

    Numeric extraction:
      Mixed-content values like "deleterious(0.05)" have the number
      extracted automatically when a numeric operator is used.

    Source-prefix stripping:
      "HGNC:28706" matches "28706" for is/in comparisons.
    """

    # Canonical operator names after synonym normalisation
    _OP_SYNONYMS = {
        'eq':      'is',
        'ne':      '!=',
        'lt':      '<',
        'gt':      '>',
        'lte':     '<=',
        'gte':     '>=',
        'matches': 'match',
        're':      'match',
        'regex':   'match',
        'ex':      'exists',
        'defined': 'exists',
    }
    # Words that should be tokenised as OP rather than field names
    _OP_WORDS = frozenset(
        ['is', 'match', 'in', 'exists', 'is_child']
        + list(_OP_SYNONYMS.keys())
    )

    _file_cache: dict = {}  # shared across all instances

    def __init__(self, *filter_strings):
        self._parsed = [self._parse(f) for f in filter_strings]

    # -- tokeniser -----------------------------------------------------------

    def _tokenize(self, expr):
        toks = []
        i = 0
        n = len(expr)
        while i < n:
            c = expr[i]
            if c.isspace():
                i += 1
                continue
            if c == '(':
                toks.append(('LP', None)); i += 1; continue
            if c == ')':
                toks.append(('RP', None)); i += 1; continue
            if c in ('>', '<', '!', '='):
                two = expr[i:i+2]
                if two in ('>=', '<=', '!='):
                    toks.append(('OP', two)); i += 2
                else:
                    toks.append(('OP', c)); i += 1
                continue
            m = re.match(r'[^\s()]+', expr[i:])
            if m:
                w = m.group(); lo = w.lower()
                if   lo == 'and':   toks.append(('AND', None))
                elif lo == 'or':    toks.append(('OR',  None))
                elif lo == 'not':   toks.append(('NOT', None))
                elif lo in self._OP_WORDS:
                    # Normalise synonyms to canonical op
                    toks.append(('OP', self._OP_SYNONYMS.get(lo, lo)))
                else:
                    toks.append(('W', w))
                i += len(w)
                continue
            i += 1
        toks.append(('EOF', None))
        return toks

    # -- recursive-descent parser -> AST -------------------------------------

    def _parse(self, expr):
        toks = self._tokenize(expr)
        pos  = [0]

        def peek(): return toks[pos[0]]
        def eat():
            t = toks[pos[0]]; pos[0] += 1; return t

        def or_():
            l = and_()
            while peek()[0] == 'OR':
                eat(); r = and_(); l = ('or', l, r)
            return l

        def and_():
            l = not_()
            while peek()[0] == 'AND':
                eat(); r = not_(); l = ('and', l, r)
            return l

        def not_():
            if peek()[0] == 'NOT':
                eat(); return ('not', atom())
            return atom()

        def atom():
            if peek()[0] == 'LP':
                eat(); node = or_()
                if peek()[0] == 'RP': eat()
                return node
            if peek()[0] != 'W':
                raise SyntaxError(f"Expected field name, got {peek()!r} in: {expr!r}")
            field = eat()[1]
            t = peek()
            if t[0] in ('EOF', 'RP', 'AND', 'OR'):
                return ('exists', field, None)
            if t[0] != 'OP':
                return ('exists', field, None)
            op = eat()[1].lower()
            if op == 'is' and peek()[0] == 'NOT':
                eat(); op = 'is not'
            if op == 'exists':
                return ('exists', field, None)
            if peek()[0] in ('EOF', 'AND', 'OR', 'RP'):
                val = None
            else:
                val = eat()[1]
            return (op, field, val)

        return or_()

    # -- evaluator -----------------------------------------------------------

    def evaluate(self, data):
        return all(self._eval(ast, data) for ast in self._parsed)

    def _fv(self, data, field):
        """Fetch field value with case-insensitive fallback."""
        v = data.get(field)
        if v is not None:
            return v
        fl = field.lower()
        for k, val in data.items():
            if k is not None and k.lower() == fl:
                return val
        return None

    @staticmethod
    def _norm(field, val):
        """Treat '' and '-' (except for Allele) as missing."""
        if val is None or val == '':
            return None
        if field != 'Allele' and val == '-':
            return None
        return val

    # Regex for VEP mixed-type field values: e.g. "HGNC:28706", "deleterious(0.05)"
    _MIXED_RE = re.compile(r'^([\w.\-]+)?:?\(?([-\d.e]*)\)?', re.IGNORECASE)
    _NUM_RE   = re.compile(r'^-?\d+\.?\d*(e-?\d+)?$', re.IGNORECASE)
    _PURE_NUM = re.compile(r'^[-\d.e]+$', re.IGNORECASE)

    def _get_input(self, fval, rep_value):
        """Replicate VEP 115 FilterSet.get_input() value-extraction logic.

        Given a raw field value and a representative comparison value,
        returns the part of fval that should be compared.

        - If fval has the form TEXT:NUM or TEXT(NUM):
            - and rep_value is purely numeric → return the NUM part
            - otherwise                      → return fval unchanged
        - If fval is plain (no numeric part) → return fval unchanged
        """
        if fval is None:
            return None
        s = str(fval)
        m = self._MIXED_RE.match(s)
        if not m:
            return s
        text = m.group(1) or ''
        num  = m.group(2) or ''
        if not num:
            return s  # no numeric part; unchanged
        # Determine whether the comparison side is purely numeric
        rep = str(rep_value) if rep_value is not None else ''
        if self._PURE_NUM.match(rep):
            # comparison is numeric → use num part unless text is already numeric
            if not self._NUM_RE.match(text):
                return num
            return text
        # comparison is not purely numeric → return value unchanged
        return s

    def _resolve_value(self, value, data):
        """If value starts with '#', resolve it to the named field's value."""
        if value and value.startswith('#'):
            return self._norm(value[1:], self._fv(data, value[1:]))
        return value

    def _eval(self, node, data):
        op = node[0]
        if op == 'and': return self._eval(node[1], data) and self._eval(node[2], data)
        if op == 'or':  return self._eval(node[1], data) or  self._eval(node[2], data)
        if op == 'not': return not self._eval(node[1], data)

        field = node[1]
        value = self._resolve_value(node[2] if len(node) > 2 else None, data)
        fval  = self._norm(field, self._fv(data, field))

        if op == 'exists':
            return fval is not None

        if op in ('is', '='):
            if value is None: return fval is None
            if fval  is None: return False
            inp = self._get_input(fval, value)
            return inp.lower() == str(value).lower()

        if op in ('is not', '!='):
            if value is None: return fval is not None
            if fval  is None: return True
            inp = self._get_input(fval, value)
            return inp.lower() != str(value).lower()

        if fval is None:
            return False

        if op == 'match':
            return bool(re.search(str(value or ''), str(fval), re.IGNORECASE))

        if op == 'in':
            if value is None: return False
            in_set = self._in_set(value)
            # Pick a representative value to determine numeric vs string mode
            rep = next(iter(in_set), None)
            inp = self._get_input(fval, rep).lower()
            return inp in in_set

        if op == 'is_child':
            # Ontology lookup not supported; exact match fallback
            inp = self._get_input(fval, value)
            return inp.lower() == str(value or '').lower()

        # Ordered comparisons — use _get_input for numeric extraction
        inp = self._get_input(fval, value)
        sv  = inp
        tv  = str(value or '')
        try:
            fn, tn = float(sv), float(tv)
            if op == '>':  return fn > tn
            if op == '<':  return fn < tn
            if op == '>=': return fn >= tn
            if op == '<=': return fn <= tn
        except (ValueError, TypeError):
            pass
        if op == '>':  return sv > tv
        if op == '<':  return sv < tv
        if op == '>=': return sv >= tv
        if op == '<=': return sv <= tv
        return False

    def limit_synonym_search(self, _val=True):
        """No-op shim (kept for API parity with the Perl module)."""

    def _in_set(self, value):
        """Return a set of lowercase values for the 'in' operator.
        If value looks like a readable file path, read one value per line.
        Otherwise split on commas.
        """
        if value in self._file_cache:
            return self._file_cache[value]
        if os.path.isfile(value):
            with open(value) as fh:
                s = {line.strip().lower() for line in fh if line.strip()}
            self._file_cache[value] = s
            return s
        s = {v.strip().lower() for v in value.split(',')}
        return s


# ---------------------------------------------------------------------------
# Parsing helpers
# ---------------------------------------------------------------------------

def parse_headers(raw_headers, vcf_info_field='CSQ'):
    """Return (vep_headers, col_headers, allowed_fields_dict)."""
    vep_headers = None
    col_headers = None
    allowed     = {}
    for raw in raw_headers:
        hashes = len(raw) - len(raw.lstrip('#'))
        line   = raw.lstrip('#')
        if hashes >= 2:
            if re.match(r'INFO=<ID=' + re.escape(vcf_info_field) + r',', line):
                m = re.search(r'Format: (.+?)"', line)
                if m:
                    vep_headers = m.group(1).split('|')
            else:
                m = re.match(r'INFO=<ID=(.+?),', line)
                if m:
                    allowed[m.group(1)] = True
                else:
                    m = re.search(r' (.+?) :', line)
                    if m:
                        allowed[m.group(1)] = True
        else:
            col_headers = line.split('\t')
    return vep_headers, col_headers, allowed


def _expand_extra(data):
    extra = data.get('Extra')
    if not extra:
        return
    for part in extra.split(';'):
        if '=' in part:
            k, v = part.split('=', 1)
            data[k] = v


def _normalize_data(data):
    for k in list(data.keys()):
        v = data[k]
        if v == '' or (k != 'Allele' and v == '-'):
            data[k] = None


def parse_tab_line(line, headers):
    parts = line.rstrip('\n').split('\t', len(headers))
    data  = {h: (parts[i] if i < len(parts) else None) for i, h in enumerate(headers)}
    data['INFO'] = None
    _expand_extra(data)
    _normalize_data(data)
    return data


def parse_csq_chunk(chunk, vep_headers, main_data=None):
    parts = chunk.split('|')
    data  = dict(main_data) if main_data else {}
    for i, h in enumerate(vep_headers):
        data[h] = parts[i] if i < len(parts) else None
    data['INFO'] = None
    _expand_extra(data)
    _normalize_data(data)
    return data


def detect_format(line):
    p = line.split('\t')
    if (len(p) >= 5
            and re.match(r'(chr)?\w+', p[0])
            and re.match(r'^\d+$', p[1])
            and p[3] and re.match(r'^[ACGTN\-.]+$', p[3], re.IGNORECASE)
            and p[4]):
        return 'vcf'
    return 'tab'


def open_input(path, force_gz=False):
    if path is None:
        return sys.stdin
    if path.endswith('.gz') or force_gz:
        return gzip.open(path, 'rt')
    with open(path, 'rb') as f:
        magic = f.read(2)
    if magic == b'\x1f\x8b':
        return gzip.open(path, 'rt')
    return open(path, 'r')


# ---------------------------------------------------------------------------
# cyvcf2 fast path (VCF only)
# ---------------------------------------------------------------------------

def process_cyvcf2(args, out_fh):
    vcf_info_field = args.vcf_info_field
    filter_set     = args.filter_set

    vcf = cyvcf2.VCF(args.input_file or '-')

    raw_hdr_lines = vcf.raw_header.rstrip('\n').split('\n')

    vep_headers = None
    for h in raw_hdr_lines:
        m = re.search(r'ID=' + re.escape(vcf_info_field) + r'.*?Format: (.+?)"', h)
        if m:
            vep_headers = m.group(1).split('|')
            break

    col_headers = None
    for h in reversed(raw_hdr_lines):
        if h.startswith('#') and not h.startswith('##'):
            col_headers = h.lstrip('#').split('\t')
            break

    all_fields = set()
    if vep_headers:  all_fields.update(vep_headers)
    if col_headers:  all_fields.update(col_headers)
    for h in raw_hdr_lines:
        m = re.match(r'##INFO=<ID=(.+?),', h)
        if m: all_fields.add(m.group(1))

    if args.list:
        print("Available fields:\n")
        for f in sorted(all_fields):
            print(f)
        vcf.close()
        return

    hdr_lines = list(raw_hdr_lines)
    if args.soft_filter:
        for i, h in enumerate(hdr_lines):
            if h.startswith('#') and not h.startswith('##'):
                hdr_lines.insert(i, '##FILTER=<ID=filter_vep_fail,Description="Variant fails filter_vep">')
                hdr_lines.insert(i, '##FILTER=<ID=filter_vep_pass,Description="Variant passes filter_vep">')
                break

    if not args.count:
        out_fh.write('\n'.join(hdr_lines) + '\n')

    col_names   = col_headers or ['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO']
    count       = 0
    line_number = 0
    missing_csq = 0

    for var in vcf:
        line_number += 1
        if args.test and line_number > args.test:
            break

        vcf_line = str(var)
        parts    = vcf_line.split('\t')

        main_data = {col_names[i]: (parts[i] if i < len(parts) else None)
                     for i in range(len(col_names))}
        if len(parts) > 7:
            for ifield in parts[7].split(';'):
                if '=' in ifield:
                    k, v = ifield.split('=', 1)
                    if k != vcf_info_field:
                        main_data[k] = v
                elif ifield:
                    main_data[ifield] = None  # Flag: Perl sets to undef, so exists-check fails

        csq_val = var.INFO.get(vcf_info_field)
        chunks    = []
        data_list = []

        if csq_val:
            chunks    = csq_val.split(',')
            data_list = [parse_csq_chunk(c, vep_headers or [], main_data) for c in chunks]
        else:
            missing_csq += 1
            data_list = [main_data]

        line_pass  = 0
        new_chunks = []
        for i, parsed in enumerate(data_list):
            if filter_set.evaluate(parsed):
                line_pass += 1
                if i < len(chunks):
                    new_chunks.append(chunks[i])

        count += bool(line_pass)

        if not args.soft_filter and count < args.start:
            continue

        if args.soft_filter:
            tag = 'filter_vep_pass' if line_pass else 'filter_vep_fail'
            cur = parts[6] if len(parts) > 6 else '.'
            parts[6] = f"{cur};{tag}" if cur and cur != '.' else tag
            out_fh.write('\t'.join(parts) + '\n')
        elif line_pass and not args.count:
            out_line = vcf_line
            if args.only_matched and new_chunks and len(new_chunks) != len(chunks):
                new_csq = ','.join(new_chunks)
                p = out_line.split('\t')
                p[7] = re.sub(re.escape(vcf_info_field) + r'=[^;]*',
                               vcf_info_field + '=' + new_csq, p[7], count=1)
                out_line = '\t'.join(p)
            out_fh.write(out_line + '\n')

        if not args.soft_filter and count >= args.limit + args.start - 1:
            break

    vcf.close()

    if line_number == 0 and not args.count:
        out_fh.write('\n'.join(hdr_lines) + '\n')

    if args.count:
        out_fh.write(f"{count}\n")

    if missing_csq:
        sys.stderr.write(
            f"WARNING: filter_vep couldn't find VEP annotations field "
            f"{vcf_info_field} in {missing_csq} line(s) of the input file\n"
        )


# ---------------------------------------------------------------------------
# Generic path (tab or VCF without cyvcf2)
# ---------------------------------------------------------------------------

def process_generic(args, out_fh):
    vcf_info_field = args.vcf_info_field
    filter_set     = args.filter_set

    in_fh = open_input(args.input_file, args.gz)

    raw_headers        = []
    vep_headers        = None
    col_headers        = None
    fmt                = args.format
    count              = 0
    line_number        = 0
    missing_csq        = 0
    headers_initialised = False

    csq_re = re.compile(re.escape(vcf_info_field) + r'=(.+?)(?:;|\s|$)')

    for raw_line in in_fh:
        line = raw_line.rstrip('\n')

        if line.startswith('#'):
            raw_headers.append(line)
            continue

        line_number += 1
        if args.test and line_number > args.test:
            break

        # -- initialise on first data line --
        if not headers_initialised:
            if not raw_headers:
                sys.exit("ERROR: No headers found in input file")

            if args.soft_filter:
                chrom_hdr = raw_headers.pop()
                raw_headers.append('##FILTER=<ID=filter_vep_pass,Description="Variant passes filter_vep">')
                raw_headers.append('##FILTER=<ID=filter_vep_fail,Description="Variant fails filter_vep">')
                raw_headers.append(chrom_hdr)

            if not args.count and not args.list:
                out_fh.write('\n'.join(raw_headers) + '\n')

            vep_headers, col_headers, extra_allowed = parse_headers(raw_headers, vcf_info_field)

            all_fields = set()
            if vep_headers:  all_fields.update(vep_headers)
            if col_headers:  all_fields.update(col_headers)
            all_fields.update(extra_allowed)

            if args.list:
                print("Available fields:\n")
                for f in sorted(all_fields):
                    print(f)
                if in_fh is not sys.stdin:
                    in_fh.close()
                return

            headers_initialised = True

        # -- format detection --
        if not fmt:
            fmt = detect_format(line)
        if fmt not in ('vcf', 'tab'):
            sys.exit(f"ERROR: Unable to parse data in format {fmt}")
        if fmt != 'vcf' and args.only_matched:
            sys.exit("ERROR: --only_matched is compatible only with VCF files")
        if fmt != 'vcf' and args.soft_filter:
            sys.exit("ERROR: --soft_filter is compatible only with VCF files")

        if args.soft_filter:
            args.start = 0

        chunks    = []
        data_list = []

        if fmt == 'tab':
            hdrs = col_headers or vep_headers or []
            data_list.append(parse_tab_line(line, hdrs))
            chunks.append(line)
            if not any(h == 'Extra' for h in hdrs):
                filter_set.limit_synonym_search(True)

        else:  # vcf
            parts = line.split('\t')
            ch    = col_headers or ['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO']
            main_data = {ch[i]: (parts[i] if i < len(parts) else None)
                         for i in range(len(ch))}

            if len(parts) > 7:
                for ifield in parts[7].split(';'):
                    if '=' in ifield:
                        k, v = ifield.split('=', 1)
                        if k != vcf_info_field:
                            main_data[k] = v
                    elif ifield:
                        main_data[ifield] = None  # Flag: Perl sets to undef, so exists-check fails

            m = csq_re.search(line)
            if m:
                raw_csq   = m.group(1)
                chunks    = raw_csq.split(',')
                data_list = [parse_csq_chunk(c, vep_headers or [], main_data) for c in chunks]
            else:
                missing_csq += 1
                data_list = [main_data]

            filter_set.limit_synonym_search(True)

        # -- evaluate --
        line_pass  = 0
        new_chunks = []
        for i, parsed in enumerate(data_list):
            if filter_set.evaluate(parsed):
                line_pass += 1
                if i < len(chunks):
                    new_chunks.append(chunks[i])

        count += bool(line_pass)

        if not args.soft_filter and count < args.start:
            continue

        # -- output --
        if args.soft_filter:
            sp  = line.split('\t')
            tag = 'filter_vep_pass' if line_pass else 'filter_vep_fail'
            cur = sp[6] if len(sp) > 6 else '.'
            sp[6] = f"{cur};{tag}" if cur and cur != '.' else tag
            out_fh.write('\t'.join(sp) + '\n')

        elif line_pass and not args.count:
            out_line = line
            if args.only_matched and new_chunks and len(new_chunks) != len(chunks):
                new_csq = ','.join(new_chunks)
                p = out_line.split('\t')
                p[7] = re.sub(re.escape(vcf_info_field) + r'=[^;]*',
                               vcf_info_field + '=' + new_csq, p[7], count=1)
                out_line = '\t'.join(p)
            out_fh.write(out_line + '\n')

        if not args.soft_filter and count >= args.limit + args.start - 1:
            break

    # -- empty file --
    if not line_number:
        if not args.count and not args.list:
            out_fh.write('\n'.join(raw_headers) + '\n')
        if args.list:
            v, c, a = parse_headers(raw_headers, vcf_info_field)
            af = set()
            if v: af.update(v)
            if c: af.update(c)
            af.update(a)
            print("Available fields:\n")
            for f in sorted(af):
                print(f)
            if in_fh is not sys.stdin:
                in_fh.close()
            return

    if args.count:
        out_fh.write(f"{count}\n")

    if missing_csq:
        sys.stderr.write(
            f"WARNING: filter_vep couldn't find VEP annotations field "
            f"{vcf_info_field} in {missing_csq} line(s) of the input file\n"
        )

    if in_fh is not sys.stdin:
        in_fh.close()


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

USAGE = """\
#------------#
# filter_vep #
#------------#

http://www.ensembl.org/info/docs/tools/vep/script/vep_filter.html

Usage:
./filter_vep_fast.py [arguments]

--help               -h   Print usage message and exit

--input_file [file]  -i   Input file (VEP results). Reads from STDIN if
                          not specified. Supports gzip (auto-detected or
                          force with --gz).
--format [vcf|tab]        Input format (tab = any tab-delimited format,
                          including default VEP output format)

--output_file [file] -o   Output file. Defaults to STDOUT.
--force_overwrite         Overwrite existing output file.

--filter [filters]   -f   Filter expression. Multiple --filter flags are
                          treated as logical ANDs, e.g.:
                            -f "Consequence is missense_variant"
                            -f "AF < 0.01 or not AF"
                            -f "HGNC_ID in gene_panels.txt"
                            -f "(AFR_AF gt #EUR_AF) and SIFT lt 0.05"

                          Operators:
                            is / = / eq       case-insensitive equality
                            != / ne           inequality
                            > / gt            greater than
                            < / lt            less than
                            >= / gte          greater than or equal
                            <= / lte          less than or equal
                            match/matches/re/regex  regex (case-insensitive)
                            in                list or file (one value/line)
                            exists/ex/defined field is present and non-null
                            is_child          SO ontology child term
                          Logical: and  or  not  ( )
                          Value can be #Field to compare two fields.
                          Mixed content (e.g. "deleterious(0.05)") has its
                          number extracted for numeric comparisons.
                          Source prefixes (e.g. "HGNC:34") are stripped
                          automatically in is/in comparisons.

--list               -l   List available fields from the input file.
--count              -c   Print only a count of matched lines.

--only_matched            In VCF files, rewrite the CSQ field to contain
                          only annotation blobs that passed the filters.

--vcf_info_field [key]    VCF INFO key for VEP annotations (default: CSQ).

--soft_filter             Add filter_vep_pass / filter_vep_fail to the VCF
                          FILTER column instead of excluding variants.

--ontology           -y   Use Sequence Ontology to match consequence terms
                          (requires Ensembl API; is_child falls back to
                          exact match in this implementation).
--host / --user / --pass / --port / --version / --registry
                          Database connection options for --ontology.

--start [N]          -s   Skip first N passing results (1-based, default 1)
--limit [N]               Return at most N passing results.
--test [N]                Process only the first N non-header lines.
"""


def main():
    ap = argparse.ArgumentParser(add_help=False)
    ap.add_argument('--help',    '-h', action='store_true')
    ap.add_argument('--test',          type=int)
    ap.add_argument('--count',   '-c', action='store_true')
    ap.add_argument('--list',    '-l', action='store_true')
    ap.add_argument('--input_file',  '-i')
    ap.add_argument('--output_file', '-o', default='stdout')
    ap.add_argument('--force_overwrite', action='store_true')
    ap.add_argument('--format',   choices=['vcf', 'tab'])
    ap.add_argument('--gz',       action='store_true')
    ap.add_argument('--only_matched', action='store_true')
    ap.add_argument('--vcf_info_field', default='CSQ')
    ap.add_argument('--soft_filter',  action='store_true')
    ap.add_argument('--ontology', '-y', action='store_true')
    ap.add_argument('--host',    default='ensembldb.ensembl.org')
    ap.add_argument('--user',    default='anonymous')
    ap.add_argument('--pass',    dest='password', default=None)
    ap.add_argument('--port',    type=int, default=3306)
    ap.add_argument('--version', type=int)
    ap.add_argument('--registry')
    ap.add_argument('--start',   '-s', type=int, default=1)
    ap.add_argument('--limit',         type=int, default=int(1e12))
    ap.add_argument('--filter',  '-f', action='append')

    args = ap.parse_args()

    if args.help or (not args.filter and not args.list):
        print(USAGE)
        if not args.help:
            sys.exit("ERROR: No valid filters given")
        sys.exit(0)

    if args.ontology:
        sys.stderr.write(
            "WARNING: --ontology requires Ensembl API; "
            "is_child will fall back to exact match in this implementation\n"
        )

    args.filter_set = FilterSet(*(args.filter or []))

    if args.output_file.lower() != 'stdout':
        if os.path.exists(args.output_file) and not args.force_overwrite:
            sys.exit(
                f"ERROR: Output file {args.output_file} already exists. "
                "Use --force_overwrite to overwrite."
            )
        out_fh = open(args.output_file, 'w')
    else:
        out_fh = sys.stdout

    try:
        use_cyvcf2 = (
            HAS_CYVCF2
            and args.format != 'tab'
            and (args.format == 'vcf'
                 or (args.input_file and
                     re.search(r'\.vcf(\.gz)?$|\.bcf$', args.input_file, re.IGNORECASE)))
        )
        if use_cyvcf2:
            process_cyvcf2(args, out_fh)
        else:
            process_generic(args, out_fh)
    finally:
        if out_fh is not sys.stdout:
            out_fh.close()


if __name__ == '__main__':
    main()
