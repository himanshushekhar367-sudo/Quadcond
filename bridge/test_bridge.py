"""Software fixtures only: no API access or trained-model claims."""
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

import anndata
import numpy as np
import pandas as pd
from alphagenome.data import genome
from export_atlas import convert, write_csv
from run_quadcond import DEFAULT_ROOT, load_quadcond, validate_export, flatten, model_path, prompt_value, start_value, chromosome_value, sequence_value


class BridgeTests(unittest.TestCase):
    def test_blank_prompts_retry(self):
        for parser, entries, expected in [
            (chromosome_value, ['', 'bad', 'chr22'], 'chr22'),
            (start_value, ['', 'text', '0', '36201698'], 36201698),
            (sequence_value, ['', 'NN', 'acgt'], 'ACGT')]:
            with patch('builtins.input', side_effect=entries), patch('builtins.print'):
                self.assertEqual(prompt_value('test', parser, 'retry'), expected)

    def test_missing_export_has_actionable_message(self):
        with self.assertRaisesRegex(ValueError, 'Run --export-only'):
            validate_export(self.folder / 'missing')

    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.folder = Path(self.tmp.name)
        self.variant = genome.Variant('chr22', 101, 'G', 'A')
        self.adata = anndata.AnnData(X=np.array([[-2., 1.], [3., np.nan]]),
            obs=pd.DataFrame({'variant': [self.variant, self.variant],
                              'gene_id': ['gene1', 'gene2']}, index=['0', '1']),
            var=pd.DataFrame({'biosample_name': ['liver', 'heart']}, index=['0', '1']))
        self.details, self.summary = convert({'AVI_SCORE': self.adata}, 'chr22', 100, 'AGG')
        write_csv(self.folder / 'quadcond_scores.csv', self.summary)
        self.meta = dict(assembly='GRCh38', chromosome='chr22', start_1based=100,
            sequence='AGG', interval_0based=[99, 102], scorer='AVI_SCORE', returned_scored_snv_count=1)
        (self.folder / 'provenance.json').write_text(json.dumps(self.meta))

    def tearDown(self):
        self.tmp.cleanup()

    def test_sdk_conversion_and_export_validation(self):
        self.assertEqual(self.details[0]['raw_score'], -2.)
        self.assertEqual(len(self.details), 4)
        self.assertEqual(self.summary[0]['score'], 3.)
        self.assertEqual(validate_export(self.folder), self.meta)

    def test_coordinate_and_reference_errors(self):
        for chrom, start, seq in [('chr21', 100, 'AGG'), ('chr22', 101, 'AGG'), ('chr22', 100, 'AAA')]:
            with self.assertRaises(ValueError):
                convert({'AVI_SCORE': self.adata}, chrom, start, seq)

    def test_duplicate_score_refused(self):
        write_csv(self.folder / 'quadcond_scores.csv', self.summary * 2)
        with self.assertRaisesRegex(ValueError, 'Duplicate'):
            validate_export(self.folder)

    def test_mixed_scorer_refused(self):
        self.summary[0]['scorer'] = 'RNA_SEQ'
        write_csv(self.folder / 'quadcond_scores.csv', self.summary)
        with self.assertRaisesRegex(ValueError, 'Mixed scorers'):
            validate_export(self.folder)

    def test_real_quadcond_join_and_missing_variant(self):
        _, variants = load_quadcond(DEFAULT_ROOT)
        from quadcond.alphagenome import TableAtlas
        # Only the structural model is replaced; the actual table reader,
        # genomic join, ranking and output conversion run unchanged.
        subs = []
        for alt in ['A', 'T']:
            subs.append(dict(position=1, wild_type_base='G', mutant_base=alt,
                label='fixture', heads={'g4_tm': dict(delta=-4., strand='+',
                    motif={'state': 'retained'}, mutant={'in_domain': True})}))
        inner = dict(heads={'g4_tm': {}}, substitutions=subs, run_record={},
            wild_type={}, delta_note='fixture', ranking_note='fixture')
        with patch.object(variants.scans, 'mutation_scan', return_value=inner):
            result = variants.variant_scan(None, 'AGG', 'chr22', 100,
                atlas_source=TableAtlas.from_path(self.folder / 'quadcond_scores.csv'),
                structural_head='g4_tm', regulatory_scorer='AVI_SCORE')
        rows = flatten(result)
        scored = next(r for r in rows if r['alternate'] == 'A')
        missing = next(r for r in rows if r['alternate'] == 'T')
        self.assertEqual(scored['regulatory_score'], 3.)
        self.assertEqual(scored['structural_delta'], -4.)
        self.assertEqual(scored['position'], 101)
        self.assertIsNone(missing['regulatory_score'])
        self.assertEqual(missing['quadrant'], 'unclassified')

    def test_bad_model_fails_before_loading(self):
        assets, _ = load_quadcond(DEFAULT_ROOT)
        bad = self.folder / 'bad.joblib'
        bad.write_bytes(b'not a model')
        with self.assertRaises(assets.AssetError):
            model_path(assets, str(bad), DEFAULT_ROOT)


if __name__ == '__main__':
    unittest.main()
