"""Validate explicitly declared experimental units without inferring donors."""


def validate_design(rows, design='unspecified'):
    if design not in {'paired', 'unpaired', 'unspecified'}:
        raise ValueError('Design must be paired, unpaired or unspecified')
    if not rows:
        raise ValueError('No samples in experimental design')
    groups = sorted({row['group'] for row in rows})
    donors = {}
    units = set()
    for row in rows:
        donor = str(row.get('donor') or '').strip()
        if design == 'paired' and not donor:
            raise ValueError(f"Paired design requires donor for {row['sample']}")
        if not donor:
            continue
        if (row['group'], donor) in units:
            raise ValueError(f'Duplicate donor within condition: {donor}/{row["group"]}; technical libraries are unsupported')
        units.add((row['group'], donor))
        donors.setdefault(donor, set()).add(row['group'])
    if design == 'paired':
        if len(groups) < 2:
            raise ValueError('Paired design requires at least two conditions')
        for donor, present in donors.items():
            if present != set(groups):
                raise ValueError(f'Incomplete pairing for {donor}: missing {", ".join(sorted(set(groups) - present))}')
    if design == 'unpaired':
        reused = [donor for donor, present in donors.items() if len(present) > 1]
        if reused:
            raise ValueError('Unpaired design reuses donors across conditions: ' + ', '.join(sorted(reused)))
    return {'design': design, 'conditions': groups, 'samples': len(rows),
            'declared_donors': len(donors),
            'model_formula': '~ donor + condition' if design == 'paired' else '~ condition' if design == 'unpaired' else None,
            'status': 'metadata validation only; no differential model fitted',
            'notes': 'Replicate labels are biological replicate identifiers within a condition, not inferred donor identities. Paired validation requires complete donor coverage across declared conditions. Unspecified design must be resolved before inference.',
            'sample_metadata': [{key: row.get(key, '') for key in ('sample', 'group', 'replicate', 'donor')} for row in rows]}
