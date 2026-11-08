"""Conditional rule development uses training designs ONLY; test results stay separate."""
import csv
import html
import json
import math
from pathlib import Path
from collections import defaultdict
import numpy as np
from scipy.stats import beta, t
from .study import write_csv

COUNT_MODELS = ('baseline','cage_fixed','cage_correlation','dream_native','mundlak_native')
FEATURES = ('conditional_ewi', 'effective_information_groups', 'groups', 'largest_group_information_share')


def csv_rows(path):
    with Path(path).open(newline='', encoding='utf-8-sig') as f:
        return list(csv.DictReader(f))


def number(value):
    try:
        x = float(value)
        return x if math.isfinite(x) else None
    except (ValueError, TypeError):
        return None


def interval(rows, metric, global_null=False):
    values = np.array([v for row in rows if (v := number(row.get(metric))) is not None])
    n=len(values)
    if not n:
        return dict(mean=None,upper=None,lower_two=None,upper_two=None,mcse=None,n=0,method='unavailable')
    mean = float(values.mean())
    se=float(values.std(ddof=1)/np.sqrt(n)) if n>1 else None
    bounded=metric in ('FDP','power','coverage','availability','available')
    if metric=='FDP' and global_null:
        if not np.all((values==0)|(values==1)):
            raise ValueError('Global-null FDP must be binary')
        k=int(values.sum())
        upper=float(beta.ppf(.95,k+1,n-k)) if k<n else 1.
        lower=float(beta.ppf(.025,k,n-k+1)) if k else 0.
        two_upper=float(beta.ppf(.975,k+1,n-k)) if k<n else 1.
        method='exact_binomial_global_null'
    elif metric=='FDP' and not np.any(values):
        # For general FDP, E(FDP) <= P(FDP>0). This is an upper bound, not a
        # binomial model for fractional FDP or a claim that all datasets are null.
        upper=1.-.05**(1./n);lower=0.;two_upper=1.-.025**(1./n)
        method='zero_event_conservative_bound'
    elif n>1:
        upper=mean+float(t.ppf(.95,n-1))*se
        lower=mean-float(t.ppf(.975,n-1))*se
        two_upper=mean+float(t.ppf(.975,n-1))*se
        if bounded:upper=min(1.,upper);lower=max(0.,lower);two_upper=min(1.,two_upper)
        method='approximate_t_monte_carlo'
    else:
        upper=1. if bounded else None;lower=0. if bounded else None;two_upper=upper
        method='insufficient_repetitions'
    return dict(mean=mean,upper=upper,lower_two=lower,upper_two=two_upper,mcse=se,n=n,method=method)


def summarize(rows, metric, global_null=False):
    z=interval(rows,metric,global_null)
    return z['mean'],z['upper'],z['n']


def aggregate(rows, phase, scenarios=None):
    scenario_by_id={int(s['scenario']):s for s in scenarios or []}
    grouped = defaultdict(list)
    for r in rows:
        grouped[(r['design'], r['split'], int(r['scenario']), r['model'])].append(r)
    result = []
    for (design, split, scenario, model), items in sorted(grouped.items()):
        if scenario in scenario_by_id:
            global_null=float(scenario_by_id[scenario]['de_fraction'])==0
        elif all(number(r.get('de_fraction')) is not None for r in items):
            fractions={float(r['de_fraction']) for r in items}
            if len(fractions)!=1:raise ValueError('Inconsistent nonnull fraction')
            global_null=next(iter(fractions))==0
        else:
            raise ValueError('Scenario metadata is required to identify the global null')
        uncertainty=interval(items,'FDP',global_null)
        fdr,upper,n=uncertainty['mean'],uncertainty['upper'],uncertainty['n']
        power, _, _ = summarize(items, 'power')
        coverage, _, _ = summarize(items, 'coverage')
        availability, _, _ = summarize(items, 'availability' if phase=='counts' else 'available')
        bias, _, _ = summarize(items, 'bias_at_positive_gamma')
        prediction, _, _ = summarize(items, 'predicted_bias_at_positive_gamma')
        result.append(dict(phase=phase,design=design,split=split,scenario=scenario,model=model,
            planned_datasets=len(items),usable_datasets=n,empirical_fdr=fdr,
            fdr_one_sided_95_mc_upper=upper,
            fdr_two_sided_95_mc_lower=uncertainty['lower_two'],
            fdr_two_sided_95_mc_upper=uncertainty['upper_two'],fdr_mcse=uncertainty['mcse'],
            fdr_interval_method=uncertainty['method'],global_null=global_null,
            power=power,coverage=coverage,availability=availability,
            empirical_positive_gamma_bias=bias,analytic_positive_gamma_bias=prediction))
    return result


def vector(design):
    return np.array([float(design.get(k) or 0.) for k in FEATURES])


def fit_rule(summary, designs, tolerance=.075, minimum_reps=20, validation_only=False):
    # No holdout metrics are serialized into the frozen rule.
    training = [dict(r) for r in summary if r['split']=='training' and r['phase']=='counts']
    ds = {d['design']:d for d in designs if d['split']=='training'}
    values = np.array([vector(d) for d in ds.values()])
    scale = np.maximum(values.std(axis=0), 1e-8)
    return dict(training_metrics=training,training_designs=ds,features=FEATURES,scale=scale.tolist(),
                fdr_mc_upper_screen=tolerance,minimum_usable_reps=minimum_reps,
                minimum_availability=.99,neighbors=2,validation_only=validation_only,
                interpretation='Conditional empirical screening, not a guarantee of FDR control',
                requires_generating_context=True,unknown_context_policy='abstain',
                mc_interval='Global null: exact binomial; other FDP: approximate t, with a conservative zero-event bound',
                formula='Nearest training designs in scaled design geometry; maximize training power among eligible models')


def choose(rule, design, scenario):
    if rule.get('validation_only'):
        return dict(model=None,status='pilot_not_for_selection',neighbors='')
    if scenario is None:
        return dict(model=None,status='generating_context_required',neighbors='')
    if scenario not in {r['scenario'] for r in rule['training_metrics']}:
        return dict(model=None,status='unsupported_generating_context',neighbors='')
    if not design['within_estimable']:
        return dict(model=None,status='within_effect_not_estimable',neighbors='')
    if design['fixed_residual_df'] < 5:
        return dict(model=None,status='insufficient_residual_df',neighbors='')
    ds = rule['training_designs']
    distances = sorted((float(np.linalg.norm((vector(design)-vector(d))/rule['scale'])), name)
                       for name,d in ds.items() if d['within_estimable'])
    names = [name for _,name in distances[:rule['neighbors']]]
    candidates = []
    for model in COUNT_MODELS:
        rows = [r for r in rule['training_metrics'] if r['design'] in names and
                r['scenario']==scenario and r['model']==model]
        if len(rows)!=len(names) or not rows:
            continue
        eligible = all(r['usable_datasets']>=rule['minimum_usable_reps'] and
                       r['availability'] is not None and r['availability']>=rule['minimum_availability'] and
                       r['fdr_one_sided_95_mc_upper'] is not None and
                       r['fdr_one_sided_95_mc_upper']<=rule['fdr_mc_upper_screen'] for r in rows)
        if eligible:
            powers = [r['power'] for r in rows if r['power'] is not None]
            # Global-null cases use error rate rather than an undefined power.
            score = float(np.mean(powers)) if powers else -float(np.mean([r['empirical_fdr'] for r in rows]))
            candidates.append((score, -max(r['fdr_one_sided_95_mc_upper'] for r in rows), model))
    return dict(model=max(candidates)[2] if candidates else None,
                status='conditionally_eligible' if candidates else 'no_supported_selection',
                neighbors=';'.join(names))


def paired_mean_shift_contrasts(rows, scenarios):
    lookup={(r['design'],int(r['scenario']),r['model'],int(r['replicate'])):r for r in rows}
    results=[]
    for s in scenarios:
        if s['cage_mean_amplitude']==0:continue
        reference=next(c for c in scenarios if c['cage_mean_amplitude']==0 and
                       c['cage_log2_sd']==s['cage_log2_sd'] and c['de_fraction']==s['de_fraction'])
        grouped=defaultdict(list)
        for r in rows:
            if int(r['scenario'])!=s['scenario']:continue
            ref=lookup[(r['design'],reference['scenario'],r['model'],int(r['replicate']))]
            if number(r['FDP']) is not None and number(ref['FDP']) is not None:
                grouped[(r['design'],r['split'],r['model'])].append(float(r['FDP'])-float(ref['FDP']))
        for (design,split,model),differences in sorted(grouped.items()):
            n=len(differences);mean=float(np.mean(differences))
            mcse=float(np.std(differences,ddof=1)/np.sqrt(n)) if n>1 else None
            results.append(dict(design=design,split=split,model=model,scenario=s['scenario'],
                reference_scenario=reference['scenario'],paired_usable_datasets=n,
                mean_FDP_difference=mean,paired_mcse=mcse,
                interpretation='Same design and paired streams; changing the generating cage mean'))
    return results


def build_report(out, protocol, supplied_rows=None, sensitivity_rows=None, include_dataset_metrics=True):
    summaries=[]
    for phase in ('gaussian','counts'):
        rows=list(supplied_rows.get(phase,[])) if supplied_rows is not None else []
        jobs=[j for j in protocol['jobs'] if j['phase']==phase]
        if supplied_rows is None:
            for job in jobs:
                rows.extend(csv_rows(out/'tasks'/job['id']/'metrics.csv'))
        if rows:
            table=aggregate(rows,phase,protocol['scenarios']);summaries.extend(table)
            write_csv(out/'review'/f'{phase}_summary.csv',table)
            # Dataset-level metrics remain compact and allow independent MC calculations.
            if include_dataset_metrics:write_csv(out/'review'/f'{phase}_dataset_metrics.csv',rows)
            write_csv(out/'review'/f'{phase}_paired_mean_shift.csv',paired_mean_shift_contrasts(rows,protocol['scenarios']))
    counts=[r for r in summaries if r['phase']=='counts']
    decisions=[]
    if counts:
        rule=fit_rule(counts,protocol['designs'],protocol['fdr_screen'],validation_only=protocol['profile']=='pilot')
        (out/'review'/'frozen_training_rule.json').write_text(json.dumps(rule,indent=2,allow_nan=False),encoding='utf-8')
        for d in protocol['designs']:
            if d['split']!='holdout':continue
            for scenario in protocol['scenarios']:
                choice=choose(rule,d,scenario['scenario'])
                actual=[r for r in counts if r['design']==d['design'] and r['scenario']==scenario['scenario']
                        and r['model']==choice['model']]
                row=dict(design=d['design'],scenario=scenario['scenario'],**choice,
                         generating_gamma=scenario['cage_mean_amplitude'],generating_sd=scenario['cage_log2_sd'],
                         generating_nonnull_fraction=scenario['de_fraction'],nominal_fdr=.05,
                         screening_mc_upper_limit=protocol['fdr_screen'],requires_generating_context=True)
                if actual:
                    a=actual[0]
                    row.update(heldout_fdr=actual[0]['empirical_fdr'],heldout_power=actual[0]['power'],
                               heldout_coverage=actual[0]['coverage'],heldout_availability=actual[0]['availability'],
                               heldout_fdr_mc_upper=a['fdr_one_sided_95_mc_upper'],
                               heldout_fdr_ci_lower=a['fdr_two_sided_95_mc_lower'],
                               heldout_fdr_ci_upper=a['fdr_two_sided_95_mc_upper'],
                               heldout_interval_method=a['fdr_interval_method'],
                               heldout_usable_datasets=a['usable_datasets'],
                               heldout_point_exceeds_nominal=(a['empirical_fdr']>.05) if a['empirical_fdr'] is not None else None,
                               heldout_point_exceeds_screen=(a['empirical_fdr']>protocol['fdr_screen']) if a['empirical_fdr'] is not None else None,
                               heldout_mc_upper_exceeds_screen=(a['fdr_one_sided_95_mc_upper']>protocol['fdr_screen']) if a['fdr_one_sided_95_mc_upper'] is not None else None)
                decisions.append(row)
        write_csv(out/'review'/'heldout_rule_evaluation.csv',decisions)
    write_csv(out/'review'/'design_diagnostics.csv',protocol['designs'])
    if protocol.get('empirical_design'):
        write_csv(out/'review'/'empirical_design_diagnostics.csv',[protocol['empirical_design']])
    sensitivity=list(sensitivity_rows or [])
    if supplied_rows is None:
        for job in protocol['jobs']:
            if job['phase']=='gaussian':
                sensitivity.extend(csv_rows(out/'tasks'/job['id']/'sensitivity.csv'))
    unique={tuple(sorted(r.items())):r for r in sensitivity}
    write_csv(out/'review'/'analytic_sensitivity.csv',unique.values())
    if supplied_rows is None:
        for name in ('diagnostics.csv','retry_attempts.csv','messages.csv','unavailable_genes.csv','failure_counts.csv'):
            records=[]
            for job in protocol['jobs']:
                path=out/'tasks'/job['id']/name
                if job['phase']=='counts' and path.is_file():records.extend(csv_rows(path))
            write_csv(out/'review'/name,records)
    issues=[]
    design_by_name={d['design']:d for d in protocol['designs']}
    for r in counts:
        expected = r['model'] in ('baseline','cage_correlation','dream_native') or (
            design_by_name[r['design']]['within_estimable'] and design_by_name[r['design']]['fixed_residual_df']>=5)
        if expected and r['usable_datasets']==0:
            issues.append(f"Entire arm unavailable: {r['design']} / {r['scenario']} / {r['model']}")
    selection_count=sum(r['model'] is not None for r in decisions)
    lines=['# AdjScope-Tx design extension', '', f"Profile: {protocol['profile']}",
           'Scope: design diagnostics and assumption-conditional sensitivity.', '',
           f"Completed planned jobs: {len(protocol['jobs'])}",
           f"Gaussian datasets: {sum(j['end']-j['start'] for j in protocol['jobs'] if j['phase']=='gaussian')}",
           f"Count datasets: {sum(j['end']-j['start'] for j in protocol['jobs'] if j['phase']=='counts')}",
           f"Held-out conditional selections: {selection_count}/{len(decisions)}",'',
           '## Interpretation',
           '- EWI is a within/total information ratio, not an unconditional FDR predictor.',
           '- Gaussian known-covariance GLS is a mechanism check, not dream.',
           '- The count benchmark uses filtering, TMM, voom/dream weights, fitting, eBayes and BH.',
           '- Complete held-out design families are excluded from rule development.',
           '- The selection rule is conditioned on the stated generating gamma, SD and nonnull fraction.',
           '- These generating values are unknown in empirical data; no empirical FDR estimate is produced.',
           f"- The MC upper screen is {protocol['fdr_screen']}; nominal BH level is 0.05.",
           '- Global-null MC bounds use exact binomial intervals; fractional FDP uses approximate t intervals.',
           '- All-zero fractional FDP uses a conservative zero-event upper bound.',
           '- Unknown generating context leads to abstention; this is not a design-only FDR selector.',
           '- The empirical screen provides no unconditional FDR guarantee.',
           '- Pilot results use reduced genes/repetitions and do not support model recommendations.',
           '- No Shapley FDR attribution, Hausman-triggered selection or new machine-learning estimator was run.',
           '', '## Whole-arm failures', *(issues or ['None detected in the structurally estimable count arms.'])]
    (out/'review'/'REPORT.md').write_text('\n'.join(lines)+'\n',encoding='utf-8')
    header=''.join(f'<th>{html.escape(k)}</th>' for k in ['Design','Split','EWI','Conditional EWI','Effective cages','Status'])
    body=''.join('<tr>'+''.join(f'<td>{html.escape(str(v))}</td>' for v in
         (d['design'],d['split'],round(d['ewi'],4) if d['ewi'] is not None else '',
          round(d['conditional_ewi'],4) if d['conditional_ewi'] is not None else '',
          round(d['effective_information_groups'],2),d['design_status']))+'</tr>' for d in protocol['designs'])
    page=f'''<!doctype html><meta charset="utf-8"><title>AdjScope-Tx design diagnostics</title>
<style>body{{font:16px Arial;margin:36px;color:#222;max-width:1200px}}table{{border-collapse:collapse;width:100%}}td,th{{padding:9px;border-bottom:1px solid #ddd;text-align:left}}th{{background:#e6f2ef}}pre{{white-space:pre-wrap;line-height:1.5}}</style>
<h1>AdjScope-Tx design diagnostics</h1><table><tr>{header}</tr>{body}</table><pre>{html.escape(chr(10).join(lines))}</pre>'''
    (out/'review'/'report.html').write_text(page,encoding='utf-8')
    return issues
