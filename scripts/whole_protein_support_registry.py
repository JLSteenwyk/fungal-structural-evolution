"""Read a fully verified support registry and require exact fit-input binding."""
import json
from pathlib import Path
from screen_duplication_alignment_reuse import sha


class SupportRegistry:
    def __init__(self, root, proof):
        root = Path(root)
        proof = Path(proof)
        validation = json.loads(proof.read_text())
        if validation['status'] != 'passed_full_resolved_whole_protein_support_readback':
            raise ValueError('Resolved support has not passed full saved-output verification')
        receipt_path = root/'receipt.json'
        if sha(receipt_path) != validation['source_receipt_sha256']:
            raise ValueError('Support receipt changed after verification')
        receipt = json.loads(receipt_path.read_text())
        if receipt['status'] != 'complete_resolved_whole_protein_support_pending_readback':
            raise ValueError('Unexpected support source')
        for field in ['inputs', 'geometries', 'settings', 'changed_inputs', 'changed_settings']:
            if receipt[field] != validation[field]:
                raise ValueError('Verification scope differs: '+field)
        for path, digest in validation['source_hashes'].items():
            if sha(path) != digest:
                raise ValueError('Verified support dependency changed: '+path)
        for name, digest in receipt['artifacts'].items():
            if sha(root/name) != digest:
                raise ValueError('Support artifact changed: '+name)
        self.rows = {}
        with (root/'resolved_input_support.jsonl').open() as stream:
            for line in stream:
                row = json.loads(line)
                identifier = row['fit_input_id']
                if identifier in self.rows:
                    raise ValueError('Repeated support input: '+identifier)
                if validation['source_hashes'].get(row['certificate_file']) != row['certificate_file_sha256']:
                    raise ValueError('Unverified certificate source: '+identifier)
                self.rows[identifier] = row
        if len(self.rows) != receipt['inputs']:
            raise ValueError('Incomplete support input registry')
        self.receipt_sha256 = sha(receipt_path)
        self.proof_sha256 = sha(proof)

    def resolve(self, fit_input_id, input_sha256):
        """Return support provenance, not statistical/scientific eligibility."""
        row = self.rows[fit_input_id]
        if row['input_sha256'] != input_sha256:
            raise ValueError('Fit input hash does not match verified support: '+fit_input_id)
        return dict(**row, zero_reference_supported=row['classification'] == 'zero_supported_to_numeric_tolerance')
