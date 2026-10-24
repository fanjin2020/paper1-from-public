"""Synthetic regression tests for publication-specific statistical definitions."""
import importlib.util
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch
import numpy as np
import pandas as pd

from pipeline import publication_tables as P


class PublicationTests(unittest.TestCase):
    def test_native_and_shared_factors_are_distinct(self):
        raw=pd.DataFrame(dict(spec=['s','s'],gene_symbol=['g','g'],branch=['baseline','cage_fixed'],
            beta_direct=[1.,2.],SE_direct_ordinary=[2.,6.],unscaled_direct=[1.,2.]))
        shared=pd.DataFrame(dict(spec=['s'],gene_symbol=['g'],SE_ratio=[3.],design_factor=[1.5],
            posterior_factor=[1.],reweight_factor=[2.],unweighted_design_prediction=[2.5]))
        z=P.native_precision(raw,shared).iloc[0]
        self.assertEqual(z.unscaled_ratio,2.)
        self.assertEqual(z.posterior_variance_ratio,2.25)
        self.assertEqual(z.design_factor,1.5)
        self.assertEqual(z.identity_ratio,z.common_weight_identity_ratio)
        bad=shared.assign(SE_ratio=4.)
        with self.assertRaisesRegex(ValueError,'disagree'):
            P.native_precision(raw,bad)

    @staticmethod
    def icc_data():
        v=pd.DataFrame(dict(spec=['s']*2,gene_symbol=['g']*2,weight_policy=['shared','native'],
            beta=[1.,2.],ICC=[.05,.1],status='available',converged=True,singular=False,at_zero=False,
            mean_within_cage_correlation=[.05,.1]))
        f=pd.DataFrame(dict(spec=['s']*4,gene_symbol=['g']*4,
            model=['baseline','dream_shared','dream_native','cage_fixed'],SE=[2.,1.,1.6,6.],beta=[1.,1.,2.,4.],status='available'))
        return v,f

    def test_icc_joins_dream_by_weight_not_fixed(self):
        v,f=self.icc_data();z=P.icc_precision(v,f).set_index('weight_policy')
        self.assertAlmostEqual(z.loc['shared','SE_ratio'],.5)
        self.assertAlmostEqual(z.loc['native','SE_ratio'],.8)
        self.assertEqual(str(z.loc['shared','ICC_bin']),'0.01-0.10')
        self.assertEqual(str(z.loc['native','ICC_bin']),'0.10-0.20')
        summary,strata=P.icc_summaries(z.reset_index())
        self.assertEqual(len(strata),8)
        self.assertIn('singular_fraction',summary)
        self.assertIn('ICC_q95',summary)

    def test_icc_missing_or_duplicate_key_stops(self):
        v,f=self.icc_data()
        with self.assertRaises(ValueError):P.icc_precision(v,f[f.model.ne('dream_native')])
        with self.assertRaisesRegex(ValueError,'Duplicate'):
            P.icc_precision(pd.concat([v,v]),f)

    def test_icc_join_difference_is_absolute_with_signed_diagnostic_separate(self):
        v,f=self.icc_data()
        v['beta']=[1.25,1.5]
        z=P.icc_precision(v,f).set_index('weight_policy')
        self.assertEqual(z.loc['shared','coefficient_join_difference'],.25)
        self.assertEqual(z.loc['shared','coefficient_join_signed_difference'],-.25)
        self.assertEqual(z.loc['native','coefficient_join_difference'],.5)
        self.assertEqual(z.loc['native','coefficient_join_signed_difference'],.5)

    def test_clinical_summary_includes_nonsignificant_and_zero_selected_contexts(self):
        frame=pd.DataFrame(dict(context=['A','A','B','B'], beta_base=1.,
            original_q_base=[.01,.7,.8,.9],log_abs_beta_ratio=[-2.,4.,1.,3.],
            minus_log_SE_ratio=[-1.,3.,2.,4.]))
        z=P.clinical_change_summary(frame).set_index('context')
        self.assertEqual(list(z.n_selected),[1,0])
        self.assertEqual(list(z.n_decomposable),[2,2])
        self.assertEqual(z.loc['A','beta_log_median'],1.)
        self.assertEqual(z.loc['A','beta_log_q25'],-.5)
        self.assertEqual(z.loc['B','beta_log_median'],2.)
        self.assertEqual(z.loc['B','se_log_q75'],3.5)

    def test_clinical_summary_excludes_nonfinite_log_pairs(self):
        frame=pd.DataFrame(dict(context='A',beta_base=[1.,0.,1.,np.nan],
            original_q_base=[.8,.2,.9,np.nan],log_abs_beta_ratio=[2.,np.inf,5.,np.nan],
            minus_log_SE_ratio=[3.,4.,np.nan,np.nan]))
        z=P.clinical_change_summary(frame).iloc[0]
        self.assertEqual(z.n_decomposable,1)
        self.assertEqual(z.n_selected,0)
        self.assertEqual(z.beta_log_median,2.)
        self.assertEqual(z.se_log_median,3.)

    def test_layout_membership_uses_results_not_feasibility(self):
        layout=pd.DataFrame(dict(layout=['both','three','unrun'],feasible=True))
        z=P.layout_membership(layout,pd.Series(['both']*5),pd.Series(['both','three']*3))
        self.assertEqual(list(z.simulated_five_model),[True,False,False])
        self.assertEqual(list(z.simulated_three_model),[True,True,False])
        pd.testing.assert_series_equal(z.simulated_here,z.simulated_five_model,check_names=False)
        with self.assertRaisesRegex(ValueError,'Unregistered'):
            P.layout_membership(layout,pd.Series(['unknown']),pd.Series(['three']))

    def test_external_interactions_are_not_disease_direction_comparisons(self):
        from pipeline.supplementary.common import discovery_direction_comparison
        from pipeline import reporting
        from pipeline.supplementary import driver
        self.assertIs(reporting.discovery_direction_comparison,discovery_direction_comparison)
        self.assertIs(driver.discovery_direction_comparison,discovery_direction_comparison)
        z=pd.DataFrame(dict(coefficient=['case','case','case:treat','case:treat','case','case','case','case'],
            discovery_direction=['Up','Down','Up','Down','Up',None,'Up','Up'],
            logFC=[1.,-1.,-1.,1.,2.,-1.,0.,np.nan],P_Value=[.01]*4+[np.nan,.1,.2,.2],
            same_as_discovery_direction=False))
        before=z.copy(deep=True)
        out=discovery_direction_comparison(z)
        pd.testing.assert_frame_equal(z,before)
        pd.testing.assert_frame_equal(out[z.columns.drop('same_as_discovery_direction')],z.drop(columns='same_as_discovery_direction'))
        self.assertEqual(list(out.same_as_discovery_direction.iloc[:2]),[True,True])
        self.assertTrue(out.same_as_discovery_direction.iloc[2:6].isna().all())
        self.assertFalse(out.same_as_discovery_direction.iloc[6])
        self.assertTrue(pd.isna(out.same_as_discovery_direction.iloc[7]))
        self.assertEqual(list(out.direction_comparison_status.iloc[2:4]),['not_applicable_interaction']*2)
        self.assertEqual(out.direction_comparison_status.iloc[4],'unavailable_estimate')
        self.assertEqual(out.direction_comparison_status.iloc[5],'unknown_discovery_direction')
        pd.testing.assert_frame_equal(discovery_direction_comparison(out),out)

    def test_review_output_is_short_and_does_not_repeat_run_name(self):
        from datetime import datetime
        from pipeline.result_collection import default_review_output
        root=Path('/results')/('long_experiment_name_'*8)
        out=default_review_output(root,datetime(2026,10,24,15,30,12))
        self.assertEqual(out.parent,root.parent)
        self.assertEqual(out.name,'p1_review_20261024_153012')
        self.assertLess(len(out.name),30)

    def test_figure_labels_use_family_positions_not_discovery_counts(self):
        from pipeline.figures import family_scope_labels, highlighted_paired_design
        scopes=pd.DataFrame(dict(model=['base','base','base'],scope=['all','context','context'],
            denominator='planned',planned=[17064,711,711]))
        labels=family_scope_labels(scopes)
        self.assertIn('17,064 / family',labels['all'])
        self.assertIn('711 / family',labels['context'])
        designs=pd.DataFrame(dict(n=[80,80,40,80],leverage='ordinary',weights='uniform',
            partial_correlation=[.3,.7,.3,.0]))
        self.assertEqual(list(highlighted_paired_design(designs)),[True,False,False,False])

    def test_bins_include_boundary_fits(self):
        bins=pd.cut([0,.009,.01,.099,.1,.199,.2,1],P.ICC_EDGES,right=False,labels=P.ICC_LABELS)
        self.assertEqual(list(bins),[P.ICC_LABELS[i] for i in [0,0,1,1,2,2,3,3]])

    def test_complete_bh_collections_and_missing_denominators(self):
        h=[];m=[]
        for arm in ['base','quality']:
            for gene,p in [('G',.01),('H',np.nan)]:
                h.append(dict(spec='A::clinical',context='A',contrast='clinical',model=arm,
                    gene_symbol=gene,beta=1.,SE=.2,p=p))
        for arm in ['baseline','cage_fixed','cage_correlation']:
            for gene,p in [('G',.01),('H',np.nan)]:
                m.append(dict(spec='whole_linear',model=arm,gene_symbol=gene,mouse_gene=gene.title(),beta=1.,SE=.2,p=p))
        t=P.correction_tables(pd.DataFrame(h),pd.DataFrame(m));z=t['values']
        self.assertEqual(set(z.collection),{'human_17064','mouse_2844','mouse_4266'})
        self.assertEqual(len(z),4+4+6)
        a=z[z.collection.eq('human_17064')&z.branch.eq('base')&z.gene_symbol.eq('G')].iloc[0]
        self.assertAlmostEqual(a.q__original__planned,.02)
        self.assertAlmostEqual(a.q__original__available_only,.01)
        self.assertEqual(len(t['continuous']),2+2+4)
        self.assertFalse(z.loc[z.p.isna(),'q__original__planned'].notna().any())

    def test_geometry_retains_full_factor_columns_for_subset(self):
        u=pd.DataFrame(dict(cage=[1,1,2,2,3,3,4,4],genotype=['5X','WT','5X','WT','5X','5X','WT','WT'],
            age=[4,4,8,8,4,8,4,8],sex=['M','F']*4,slide=['A','B']*4))
        z=P.mouse_geometry(u);s=z[z.spec.eq('whole_linear')].set_index('branch')
        self.assertEqual(s.loc['mixed_cages_geometry_only','columns'],s.loc['cage_fixed','columns'])
        self.assertEqual(s.loc['mixed_cages_geometry_only','n'],4)
        a=s.loc['baseline'];self.assertAlmostEqual(a.unit_variance_beta_factor,1/a.genotype_residual_SS)

    def test_reexport_environment_does_not_import_excel_or_R(self):
        import paper1_reexport_results as R
        from types import SimpleNamespace
        record=dict(python=R.platform.python_version(),packages={n:'1' for n in R.PACKAGES})
        with patch.object(R.importlib,'import_module',return_value=SimpleNamespace(__version__='1')) as p:
            R.environment(record)
        self.assertEqual([x.args[0] for x in p.call_args_list],list(R.PACKAGES))

    def test_reexport_refuses_inplace_or_existing_output(self):
        from paper1_reexport_results import safe_output
        with tempfile.TemporaryDirectory() as d:
            root=Path(d)/'run';project=Path(d)/'code';root.mkdir();project.mkdir()
            with self.assertRaises(ValueError):safe_output(root,root/'export',project)
            with self.assertRaises(ValueError):safe_output(root,project,project)
            with self.assertRaises(ValueError):safe_output(root,Path(d),project)
            safe_output(root,Path(d)/'new_export',project)

    def test_collection_handles_large_protocol_fields_when_imported(self):
        import csv
        from pipeline.result_collection import table_info
        with tempfile.TemporaryDirectory() as d:
            path=Path(d)/'protocol.csv'
            pd.DataFrame({'value':['x'*300000,'quoted\nline']}).to_csv(path,index=False)
            previous=csv.field_size_limit(131072)
            try:
                self.assertEqual(table_info(path),(2,1))
                self.assertEqual(csv.field_size_limit(),131072)
            finally:csv.field_size_limit(previous)


if __name__=='__main__':unittest.main()
