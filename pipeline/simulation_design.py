"""四种代表布局的条件模拟设计；真实布局保留年龄和性别。"""
import itertools
import numpy as np
import pandas as pd
from .io import require


def synthetic(size, mixed, fraction, seed):
    """精确满足动物数、混合笼数和基因型人数，不用近似布局。"""
    n, rng = 72, np.random.default_rng(seed)
    k, m, target = n // size, round(n // size * mixed), round(n * fraction)
    options = [p for p in range(k - m + 1) if m <= target - p * size <= m * (size - 1)]
    require(bool(options), '给定布局无法满足精确配额')
    pure = options[len(options) // 2]
    nums, left = [1] * m, target - pure * size - m
    while left:
        candidates = [j for j, v in enumerate(nums) if v < size - 1]
        low = min(nums[j] for j in candidates)
        j = int(rng.choice([j for j in candidates if nums[j] == low]))
        nums[j] += 1
        left -= 1
    nums += [size] * pure + [0] * (k - m - pure)
    rng.shuffle(nums)
    rows = []
    for j, num in enumerate(nums):
        gs = [1] * num + [0] * (size - num)
        rng.shuffle(gs)
        rows.extend(dict(cage=j, genotype=v) for v in gs)
    return pd.DataFrame(rows)


def canonical_actual(units):
    """按来源SAMPLE排序并显式编码笼位，使同一种子不随清单行序改变。"""
    require('SAMPLE' in units, '真实布局缺少来源SAMPLE，不能猜测随机数分配顺序')
    u = units.copy()
    require(u.SAMPLE.notna().all(), '来源SAMPLE为空或重复')
    u['SAMPLE'] = u['SAMPLE'].astype(str)
    require(u.SAMPLE.str.strip().ne('').all() and not u.SAMPLE.duplicated().any(), '来源SAMPLE为空或重复')
    require(u.cage.notna().all(), '真实布局存在缺失笼位')
    u = u.sort_values('SAMPLE', kind='stable').reset_index(drop=True)
    u['source_cage'] = u.cage.astype(str)
    u['cage'] = pd.Categorical(u.source_cage, categories=sorted(u.source_cage.unique())).codes
    return u


def matrices(u):
    """固定效应保留完整干扰空间；SVD只移除干扰项内部的重复列。"""
    x = pd.DataFrame(dict(Intercept=np.ones(len(u)), genotype=u.genotype.to_numpy(float)))
    if 'age' in u:
        x['age_centered'] = u.age.to_numpy(float) - u.age.mean()
        x['male'] = u.sex.eq('M').astype(float).to_numpy()
    require(np.linalg.matrix_rank(x) == x.shape[1], '基础设计秩不足')
    z = x.drop(columns='genotype').to_numpy(float)
    cage = pd.get_dummies(u.cage, drop_first=True, dtype=float).to_numpy()
    nuisance = np.column_stack([z, cage])
    q, _, _ = np.linalg.svd(nuisance, full_matrices=False)
    rank = np.linalg.matrix_rank(nuisance)
    basis = q[:, :rank]
    g = x.genotype.to_numpy()
    r0 = g - z @ np.linalg.lstsq(z, g, rcond=None)[0]
    r1 = g - basis @ (basis.T @ g)
    ss0, ss1 = float(r0 @ r0), float(r1 @ r1)
    fixed = pd.DataFrame(np.column_stack([g, basis]), columns=['genotype'] + [f'z{i}' for i in range(rank)])
    info = dict(n=len(u), cages=u.cage.nunique(),
                mixed=int(u.groupby('cage').genotype.nunique().eq(2).sum()),
                fixed_estimable=ss1 > 1e-9 and len(u) - rank - 1 > 1,
                baseline_rank=x.shape[1], fixed_rank=np.linalg.matrix_rank(fixed),
                information_fraction=ss1 / ss0)
    return x, fixed, info


def prepare(work, config):
    """生成四布局和32个情境，逐情境种子不随续跑次序变化。"""
    require(config['simulation_grid'] == 'targeted', '本版只实现已说明的四布局模拟')
    u = pd.read_csv(work / 'prepared' / 'mouse_units.csv', dtype={'cage': str, 'SAMPLE': str})
    require(u.genotype.isin(['WT', '5X']).all(), '模拟基因型映射不明确')
    u['genotype'] = u.genotype.eq('5X').astype(int)
    layouts = {'actual_linear': canonical_actual(u)}
    for index, (size, mixed, fraction) in enumerate(itertools.product((2, 4, 6), (0, .5, 1), (1 / 3, .5))):
        if size == 6 and fraction == .5:
            name = f'n72_s6_m{int(mixed * 100):03d}_g36'
            layouts[name] = synthetic(size, mixed, fraction, 3252026 + index)
    folder = work / 'prepared' / 'simulation'
    folder.mkdir(parents=True, exist_ok=True)
    geometry, tasks = [], []
    for name, units in layouts.items():
        x, fixed, info = matrices(units)
        units.to_csv(folder / f'{name}_units.csv', index=False)
        x.to_csv(folder / f'{name}_base.csv', index=False)
        fixed.to_csv(folder / f'{name}_fixed.csv', index=False)
        geometry.append(dict(layout=name, **info))
    for index, (name, rho, shift, truth) in enumerate(itertools.product(layouts, (0., .2), (0., .5), (0., .3))):
        tasks.append(dict(scenario=f's{index:03d}', layout=name, rho_truth=rho,
                          shift=shift, truth=truth, seed=config['simulation_seed'] + index * 100000,
                          planned=config['simulation_repeats']))
    pd.DataFrame(tasks).to_csv(work / 'results' / 'simulation_scenarios.csv', index=False)
    pd.DataFrame(geometry).to_csv(work / 'results' / 'simulation_layouts.csv', index=False)
