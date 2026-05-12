#!/usr/bin/env python
"""
Convert a GitHub Issue body (produced by an Issue Form in
.github/ISSUE_TEMPLATE/) into a manifest YAML.

GitHub renders Issue Forms into a markdown body where every form
field becomes a `### <label>` heading followed by a blank line and
the value. Dropdowns render the selected text; textareas render
multi-line content. This script parses that format and emits a
manifest matching the schema in Biochemistry/Updates/README.md.

The manifest type is read from a label on the issue, of the form
`manifest:<type>`. The labels are set automatically by the form
template's `labels:` field.

Usage:
    python Issue_To_Manifest.py \\
        --type structure_update \\
        --body /path/to/issue_body.md \\
        --author <github-handle> \\
        --output /path/to/manifest.yaml

Designed to run from inside .github/workflows/manifest-from-issue.yml.
"""
import argparse
import datetime
import os
import re
import sys


SUPPORTED_TYPES = {
    'structure_update',
    'override_add',
    'alias_add',
    'ignore_add',
}


def parse_form_body(body):
    """Parse a GitHub Issue Form body into {label: value}.

    Each form field is rendered as:
        ### <Label>
        \\n
        <value, possibly multi-line>
        \\n

    Empty fields render as '_No response_'.
    """
    out = {}
    lines = body.splitlines()
    i = 0
    while i < len(lines):
        line = lines[i].rstrip()
        m = re.match(r'^###\s+(.+)$', line)
        if not m:
            i += 1
            continue
        label = m.group(1).strip()
        # Skip the blank line after the heading
        i += 1
        while i < len(lines) and not lines[i].strip():
            i += 1
        # Read value until the next ### heading or EOF
        value_lines = []
        while i < len(lines) and not lines[i].startswith('### '):
            value_lines.append(lines[i])
            i += 1
        # Trim trailing blank lines
        while value_lines and not value_lines[-1].strip():
            value_lines.pop()
        value = '\n'.join(value_lines).strip()
        if value == '_No response_':
            value = ''
        out[label] = value
    return out


def required(form, label):
    val = form.get(label, '').strip()
    if not val:
        sys.exit(f'ERROR: required form field missing or empty: {label!r}')
    return val


def to_yaml_block(s, indent=2):
    """Render a possibly multi-line string as a YAML literal block."""
    pad = ' ' * indent
    return '|\n' + '\n'.join(pad + ln for ln in s.splitlines())


def yaml_scalar(s):
    """Quote a scalar for safe YAML embedding."""
    s = str(s)
    if any(c in s for c in ':#"\\\n') or s in ('', 'null', 'true', 'false'):
        return '"' + s.replace('\\', '\\\\').replace('"', '\\"') + '"'
    if re.match(r'^[+-]?\d', s):
        return '"' + s + '"'   # keep numeric-looking strings as strings
    return s


def build_structure_update(form, author, date):
    source = required(form, 'Source database')
    ext_id = required(form, 'External ID')
    field  = required(form, 'Field to update')
    val    = required(form, 'New structure value')
    reason = required(form, 'Reason')
    return f"""type: structure_update
title: "{ext_id} {field} update via Issue Form"
author: {yaml_scalar(author)}
date: {date}
reason: {to_yaml_block(reason)}
target:
  source: {source}
  external_id: {yaml_scalar(ext_id)}
  field: {field}
change:
  new_value: {yaml_scalar(val)}
"""


def build_override_add(form, author, date):
    cpd     = required(form, 'ModelSEED compound ID')
    name    = required(form, 'Compound name')
    formula = required(form, 'Formula')
    charge  = required(form, 'Charge')
    reason  = required(form, 'Reason')
    return f"""type: override_add
title: "Add {cpd} ({name}) to ACP override table"
author: {yaml_scalar(author)}
date: {date}
reason: {to_yaml_block(reason)}
target:
  file: Curation/overrides/acps_formula_charge.tsv
change:
  add_row:
    ID: {cpd}
    name: {yaml_scalar(name)}
    formula: {yaml_scalar(formula)}
    charge: {yaml_scalar(charge)}
"""


def build_alias_add(form, author, date):
    cpd     = required(form, 'ModelSEED compound ID')
    src1    = required(form, 'First source')
    ext1    = required(form, 'First external ID')
    src2    = form.get('Second source', '').strip()
    ext2    = form.get('Second external ID', '').strip()
    reason  = required(form, 'Reason')
    adds = [f"    - source: {src1}\n      external_id: {yaml_scalar(ext1)}"]
    if src2 and ext2:
        adds.append(f"    - source: {src2}\n      external_id: {yaml_scalar(ext2)}")
    return f"""type: alias_add
title: "Add alias(es) for {cpd}"
author: {yaml_scalar(author)}
date: {date}
reason: {to_yaml_block(reason)}
target:
  modelseed_id: {cpd}
change:
  add:
{chr(10).join(adds)}
"""


def build_ignore_add(form, author, date):
    ignore_file = required(form, 'Ignore file')
    ext_id      = required(form, 'External ID to ignore')
    accepted    = required(form, "Replacement structure ID (or 'None')")
    notes       = form.get('Notes', '').strip()
    reason      = required(form, 'Reason')
    return f"""type: ignore_add
title: "Ignore {ext_id}"
author: {yaml_scalar(author)}
date: {date}
reason: {to_yaml_block(reason)}
target:
  file: {ignore_file}
  create_if_missing: true
change:
  add_row:
    external_id: {yaml_scalar(ext_id)}
    accepted: {yaml_scalar(accepted)}
    notes: {yaml_scalar(notes)}
"""


BUILDERS = {
    'structure_update': build_structure_update,
    'override_add':     build_override_add,
    'alias_add':        build_alias_add,
    'ignore_add':       build_ignore_add,
}


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--type',   required=True, help='manifest type (must match a manifest:* label)')
    ap.add_argument('--body',   required=True, help='path to a file containing the issue body markdown')
    ap.add_argument('--author', required=True, help='GitHub handle of the issue author')
    ap.add_argument('--date',   default=datetime.date.today().isoformat(),
                    help='ISO date (default: today)')
    ap.add_argument('--output', required=True, help='where to write the manifest YAML')
    args = ap.parse_args()

    if args.type not in SUPPORTED_TYPES:
        sys.exit(f'ERROR: unknown manifest type {args.type!r}; supported: {sorted(SUPPORTED_TYPES)}')

    with open(args.body) as fh:
        body = fh.read()
    form = parse_form_body(body)

    builder = BUILDERS[args.type]
    yaml_str = builder(form, args.author, args.date)

    os.makedirs(os.path.dirname(args.output) or '.', exist_ok=True)
    with open(args.output, 'w') as fh:
        fh.write(yaml_str)
    print(f'wrote manifest: {args.output}')


if __name__ == '__main__':
    main()
