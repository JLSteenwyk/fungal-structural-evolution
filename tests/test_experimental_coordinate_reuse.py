import gzip,json,sys,tempfile,unittest
from pathlib import Path
from unittest.mock import patch
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from retrieve_experimental_coordinates import main
from audit_busco_gene_copies import sha

class ReuseTests(unittest.TestCase):
    def test_pinned_reuse_and_integrative_deferral(self):
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory);screen=root/'screen';reuse=root/'reuse';meta=root/'meta'
            for path in (screen,reuse,meta):path.mkdir()
            def dump(path,value):path.write_text(json.dumps(value))
            table=screen/'sequence_correspondence.tsv'
            table.write_text('entity_id\tsequence_class\n1ABC_1\texact_full_sequence\n2ABC_1\texact_full_sequence\n')
            dump(screen/'receipt.json',{'artifacts':{table.name:sha(table)}})
            entries=meta/'entries.tsv';entries.write_text('entry_id\tmethodology\n1ABC\texperimental\n2ABC\tintegrative\n')
            dump(meta/'receipt.json',{'screen_receipt_sha256':sha(screen/'receipt.json'),'artifacts':{entries.name:sha(entries)}})
            dump(reuse/'config.json',{'test':True});archive=reuse/'1ABC.cif.gz';raw=b'data_1ABC\n#\n';archive.write_bytes(gzip.compress(raw))
            rp=reuse/'1ABC.receipt.json';dump(rp,{'entry_id':'1ABC','config_sha256':sha(reuse/'config.json'),'gzip_sha256':sha(archive),'compressed_bytes':archive.stat().st_size,'uncompressed_bytes':len(raw)})
            dump(reuse/'receipt.json',{'status':'complete_frozen_experimental_coordinate_download','config_sha256':sha(reuse/'config.json'),'results':[{'entry_id':'1ABC','receipt_sha256':sha(rp),'gzip_sha256':sha(archive)}]})
            def run(output):
                args=['download','--screen',str(screen),'--output',str(output),'--reuse',str(reuse),'--reference-metadata',str(meta)]
                with patch('sys.argv',args),patch('urllib.request.urlopen',side_effect=AssertionError('Unexpected network request')):main()
            output=root/'out';run(output)
            self.assertEqual((output/archive.name).read_bytes(),archive.read_bytes())
            self.assertEqual(json.loads((output/'config.json').read_text())['excluded_entries'],{'2ABC':'integrative'})
            self.assertEqual(json.loads((output/'receipt.json').read_text())['entries'],1)
            archive.write_bytes(gzip.compress(b'data_1ABC\nchanged\n'))
            with self.assertRaises(ValueError):run(root/'tampered')

if __name__=='__main__':unittest.main()
