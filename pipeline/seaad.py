"""SEA-AD 原始 UMIs 的供体级计数准备。"""

import csv
import hashlib
import json
import re
import shutil
import urllib.parse
import xml.etree.ElementTree as ET
from pathlib import Path
import numpy as np
import pandas as pd
from .io import download

VERSION = "2.0.0-public-clinical"
S3 = "https://sea-ad-single-cell-profiling.s3.amazonaws.com/"
RELEASE = "2026-06-22"
EXPECTED_DONORS = 80
SYMBOL = re.compile("[A-Za-z0-9][A-Za-z0-9_.-]*")


def need(ok, why):
    """检查必需条件。条件不成立时抛出异常，保留调用者给出的具体原因。"""
    if not ok:
        raise ValueError(why)


def digest(p):
    """分块计算文件的 SHA-256，用于确认缓存内容是否发生变化；不把文件全部读入内存。"""
    h = hashlib.sha256()
    with Path(p).open("rb") as f:
        for b in iter(lambda: f.read(4 * 1024 * 1024), b""):
            h.update(b)
    return h.hexdigest()


def dump(p, value):
    """将字典写成 UTF-8 JSON。拒绝 NaN，避免把无效数值保存成看似有效的元数据。"""
    Path(p).write_text(
        json.dumps(value, ensure_ascii=False, indent=2, allow_nan=False), encoding="utf-8"
    )


def records(p):
    """按原列名读取 CSV，返回字典列表；兼容带 BOM 的 UTF-8 文件。"""
    with Path(p).open(encoding="utf-8-sig", newline="") as f:
        return list(csv.DictReader(f))


def write_csv(p, rows):
    """把记录列表写成 CSV，保留输入行序；空的必需记录表会报错。"""
    need(len(rows) > 0, "No records for " + str(p))
    pd.DataFrame(rows).to_csv(p, index=False)


def finite_number(value):
    """解析有限数值；缺失、非数值和无穷大返回 None，由上层决定是否纳入。"""
    try:
        x = float(value)
    except (ValueError, TypeError):
        return None
    return x if np.isfinite(x) else None


def braak(value):
    """把数字或罗马数字 Braak 分期转成 0 至 6；其他值保留为缺失。"""
    s = str(value).strip().removeprefix("Braak ").upper()
    roman = {"0": 0, "I": 1, "II": 2, "III": 3, "IV": 4, "V": 5, "VI": 6}
    x = roman.get(s, finite_number(s))
    return float(x) if x is not None and x in range(7) else None


def unit(
    uid, donor, region, context, age, sex, pathology, platform, batch="not_recorded", n_cells=0
):
    """建立一条供体与细胞亚类的伪批量记录。

    uid 是表达矩阵列名；donor 是生物学供体；n_cells 是实际回收核数。
    脑区、年龄、性别、病理和批次保留来源定义，不据此推断未知批次。"""
    return dict(
        unit_id=uid,
        donor_id=str(donor),
        region=region,
        context=context,
        age=float(age),
        sex=sex.lower(),
        pathology=float(pathology),
        platform=platform,
        batch=batch,
        replicate_role="biological",
        n_cells=int(n_cells),
    )


def gene_symbols(values):
    """检查来源基因符号，排除空标识和未转换的 Ensembl 编号。

    人类符号统一为大写，返回顺序不变；无效行用空字符串表示。
    这里只处理大小写，不把旧别名猜成当前官方符号，也不在线查询别名。"""
    return np.array(
        [
            str(v).upper() if SYMBOL.fullmatch(str(v)) and (not str(v).startswith("ENSG")) else ""
            for v in values
        ]
    )


def collapse_genes(symbols, values, mean=False):
    """按有效基因符号合并重复行。

    默认相加，适用于原始计数；返回合并后矩阵与符号。mean 模式仅用于强度数据，
    计数准备流程不使用该模式。缺少足够可识别基因时停止。"""
    symbols = gene_symbols(symbols)
    valid = symbols != ""
    need(valid.sum() > 100, "No valid source gene symbols; do not guess Ensembl conversion")
    names, inverse = np.unique(symbols[valid], return_inverse=True)
    result = np.zeros((len(names), values.shape[1]), dtype=float)
    np.add.at(result, inverse, values[valid])
    if mean:
        result /= np.bincount(inverse)[:, None]
    return (names, result)


def read_strings(dataset):
    """读取 HDF5 字符串数据集，将字节串解码为文本并保留原顺序。"""
    return np.array([x.decode() if isinstance(x, bytes) else str(x) for x in dataset[:]])


def h5_column(node, name):
    """读取 H5AD 的普通列或分类编码列。分类代码为缺失时返回空字符串。"""
    import h5py

    x = node[name]
    if isinstance(x, h5py.Dataset):
        data = x[:]
        return np.array([v.decode() if isinstance(v, bytes) else v for v in data])
    if "categories" in x and "codes" in x:
        categories = h5_column(x, "categories")
        codes = x["codes"][:]
        return np.array([categories[c] if c >= 0 else "" for c in codes])
    raise ValueError("Unsupported H5AD column: " + name)


def count_sum_sparse(group, selected, labels, n_features, cell_axis):
    """按选中细胞的分组汇总稀疏原始计数。

    selected 为细胞索引，labels 为对应分组。函数检查计数非负且为整数，
    按细胞指针分块读取；输出基因乘分组矩阵，不先构造完整稠密细胞矩阵。"""
    need(cell_axis == "rows", "本流程只接受细胞在行的 CSR 计数")
    need(len(selected) == len(labels), "分组索引与标签数量不一致")
    ptr = group["indptr"][:]
    need(len(ptr) > 0 and ptr[0] == 0 and np.all(np.diff(ptr) >= 0), "CSR 指针无效")
    need(ptr[-1] == len(group["data"]) == len(group["indices"]), "CSR 数据长度不一致")
    flat = np.asarray([i for rows in selected for i in rows], dtype=int)
    need(
        np.all((flat >= 0) & (flat < len(ptr) - 1)) and len(np.unique(flat)) == len(flat),
        "细胞索引重复或越界",
    )
    result = np.zeros((n_features, len(labels)), dtype=float)
    membership = np.full(len(ptr)-1, -1, dtype=int)
    for j, rows in enumerate(selected):
        membership[rows] = j
    # 连续读取一段CSR，避免逐核解压HDF5块造成大量重复I/O。
    for first in range(0, len(membership), 512):
        last = min(first+512, len(membership))
        if np.all(membership[first:last] < 0):
            continue
        lo, hi = int(ptr[first]), int(ptr[last])
        groups = np.repeat(membership[first:last], np.diff(ptr[first:last+1]))
        keep = groups >= 0
        data = group['data'][lo:hi].astype(float)[keep]
        ix = group['indices'][lo:hi][keep]
        need(np.all((ix >= 0) & (ix < n_features)), '基因索引超出矩阵范围')
        need(np.isfinite(data).all() and (data >= 0).all()
             and np.allclose(data, np.round(data), rtol=0, atol=1e-8),
             '选中层不是非负整数计数，不能取整或反变换后替代')
        np.add.at(result, (ix, groups[keep]), data)
    return result


def atomic_npz(path, **data):
    """先写临时 NPZ，再替换目标缓存，避免中断留下半个可复用文件。"""
    tmp = path.with_name(path.name + ".part")
    with tmp.open("wb") as f:
        np.savez_compressed(f, **data)
    tmp.replace(path)
    dump(path.with_name(path.name + ".sha256.json"), dict(sha256=digest(path)))


def verified_npz(path):
    """核对缓存键和数组信息；只复用与当前来源对象匹配的供体聚合结果。"""
    receipt = path.with_name(path.name + ".sha256.json")
    need(
        receipt.exists() and json.loads(receipt.read_text())["sha256"] == digest(path),
        "Pseudobulk cache checksum mismatch",
    )
    return np.load(path, allow_pickle=False)


def sea_inventory(cache):
    """列出指定 SEA-AD 发布日期的额叶供体对象。

    检查清单唯一性和预期供体数量；不自动改用最新发布版本。返回远端对象记录。"""
    token = None
    result = []
    page = 0
    ns = {"s": "http://s3.amazonaws.com/doc/2006-03-01/"}
    while True:
        query = {"list-type": "2", "prefix": "PFC/RNAseq/donors_objects/", "max-keys": "1000"}
        if token:
            query["continuation-token"] = token
        path = download(
            S3 + "?" + urllib.parse.urlencode(query), cache / f"SEAAD_inventory_{page}.xml"
        )
        root = ET.parse(path).getroot()
        for item in root.findall("s:Contents", ns):
            key = item.find("s:Key", ns).text
            if key.endswith(f"_SEAAD_DFC_RNAseq_final-nuclei.{RELEASE}.h5ad"):
                result.append(
                    dict(
                        key=key,
                        size=int(item.find("s:Size", ns).text),
                        etag=item.find("s:ETag", ns).text,
                    )
                )
        if root.find("s:IsTruncated", ns).text == "false":
            break
        token = root.find("s:NextContinuationToken", ns).text
        page += 1
        need(page < 100, "Unexpected inventory pagination")
    need(
        len(result) == EXPECTED_DONORS and len({r["key"] for r in result}) == EXPECTED_DONORS,
        f"所选发行版应有{EXPECTED_DONORS}个唯一供体对象",
    )
    return sorted(result, key=lambda x: x["key"])


def prepare_seaad(cache, out, run_id=None):
    """读取 SEA-AD 原始 UMIs 层并汇总到供体与作者细胞亚类。

    cache 保存下载和供体级 NPZ；out 保存计数准备所需的小型记录表。
    保留作者纳入标记，核对供体、脑区与人口学字段，返回聚合矩阵和单位记录。"""
    import h5py

    need(run_id, "Raw-input run identity is required for donor preparation")
    donor_cache = out / "donor_cache"
    donor_cache.mkdir(exist_ok=True)
    inventory = sea_inventory(cache)
    write_csv(out / "source_object_inventory.csv", inventory)
    blocks = []
    manifests = []
    excluded = []
    schemas = []
    for number, entry in enumerate(inventory, 1):
        name = Path(entry["key"]).name
        donor_file = name.split("_SEAAD_")[0]
        path = cache / name
        checkpoint = donor_cache / (donor_file + "_pseudobulk.npz")
        signature = entry["key"] + entry["etag"] + str(entry["size"]) + VERSION + run_id
        if checkpoint.exists():
            z = verified_npz(checkpoint)
            need(str(z["fingerprint"]) == signature, "Changed donor checkpoint")
            genes = z["genes"]
            data = z["values"]
            rows = json.loads(str(z["rows"]))
            audit = json.loads(str(z["audit"]))
            if "schema" in z.files:
                schemas.append(json.loads(str(z["schema"])))
        else:
            need(
                shutil.disk_usage(cache).free > (0 if path.is_file() else entry["size"]) + 2 * 2**30,
                "Insufficient disk for next donor file; free space, do not select fewer donors",
            )
            path = download(S3 + entry["key"], path)
            need(path.stat().st_size == entry["size"], "Inventory/file length mismatch")
            with h5py.File(path, "r") as f:
                obs = f["obs"]
                required = [
                    "Donor ID",
                    "Brain Region",
                    "Age at Death",
                    "Sex",
                    "Braak",
                    "Cognitive Status",
                    "Subclass",
                    "Used in analysis",
                    "Neurotypical reference",
                ]
                obs_index = h5_column(obs, obs.attrs.get("_index", "_index"))
                need(len(obs_index) == len(set(obs_index)), "Duplicate H5AD cell identifier")
                need(
                    all((k in obs for k in required)), "Changed H5AD obs schema: " + str(list(obs))
                )
                columns = {k: h5_column(obs, k) for k in required}
                n = len(columns["Donor ID"])
                need(
                    set(columns["Donor ID"].astype(str)) == {donor_file},
                    "File donor does not match observed donor",
                )
                need(set(columns["Brain Region"].astype(str)) == {"DFC"}, "Unexpected anatomy")
                identity = {}
                for k in ["Age at Death", "Sex", "Braak", "Cognitive Status", "Neurotypical reference"]:
                    vals = set(map(str, columns[k]))
                    need(len(vals) == 1, "Inconsistent donor phenotype " + k)
                    identity[k] = next(iter(vals))
                age = finite_number(identity["Age at Death"])
                pathology = braak(identity["Braak"])
                sex = identity["Sex"].lower()
                reason = (
                    "included"
                    if age is not None
                    and pathology is not None
                    and (sex in ("male", "female"))
                    and (identity["Neurotypical reference"].lower() == "false")
                    else "invalid_metadata_or_reference"
                )
                audit = dict(
                    donor=donor_file,
                    age_raw=identity["Age at Death"],
                    braak_raw=identity["Braak"],
                    status=reason,
                )
                groups = sorted(set(map(str, columns["Subclass"])) - {"", "nan", "None"})
                selected = [
                    (
                        np.flatnonzero(
                            (columns["Subclass"].astype(str) == g)
                            & (np.char.lower(columns["Used in analysis"].astype(str)) == "true")
                        )
                        if reason == "included"
                        else np.array([], dtype=int)
                    )
                    for g in groups
                ]
                matrix = f["layers/UMIs"] if "layers/UMIs" in f else None
                need(
                    matrix is not None,
                    "Fixed SEA-AD release lacks layers/UMIs; refuse to use normalized X",
                )
                encoding = matrix.attrs.get("encoding-type", "")
                encoding = encoding.decode() if isinstance(encoding, bytes) else encoding
                need(
                    encoding == "csr_matrix",
                    "Expected cell-by-gene CSR raw counts, got " + str(encoding),
                )
                var = f["var"]
                symbol_key = next(
                    (k for k in ("gene_symbol", "gene_symbols", "feature_name") if k in var), None
                )
                if symbol_key:
                    symbols = h5_column(var, symbol_key)
                else:
                    index = var.attrs.get("_index", "_index")
                    symbols = h5_column(var, index)
                need(len(matrix["indptr"]) == n + 1, "Counts/obs row mismatch")
                data = count_sum_sparse(matrix, selected, groups, len(symbols), "rows")
                genes, data = collapse_genes(symbols, data)
                rows = (
                    [
                        unit(
                            donor_file + "__" + str(j),
                            donor_file,
                            "dorsolateral_frontal_cortex",
                            g,
                            age if age is not None else 0,
                            sex if sex in ("male", "female") else "male",
                            pathology if pathology is not None else 0,
                            "10X_UMI",
                            n_cells=len(selected[j]),
                        )
                        for j, g in enumerate(groups)
                    ]
                    if reason == "included"
                    else []
                )
                if not rows:
                    data = data[:, :0]
                for row in rows:
                    row['cognitive_status_source'] = identity['Cognitive Status']
                audit["count_object"] = matrix.name
                audit["gene_symbol_field"] = symbol_key or str(var.attrs.get("_index", "_index"))
                schemas.append(
                    dict(donor=donor_file, obs_columns=list(obs), count_object=matrix.name)
                )
            atomic_npz(
                checkpoint,
                genes=genes,
                values=data,
                rows=json.dumps(rows),
                audit=json.dumps(audit),
                fingerprint=signature,
                schema=json.dumps(schemas[-1]),
            )
        excluded.append(audit)
        if rows:
            blocks.append((genes, data, [r["unit_id"] for r in rows]))
            manifests.extend(rows)
        print("SEA-AD donor", number, "/", len(inventory), donor_file, "pseudobulk complete", flush=True)
    write_csv(out / "donor_inclusion.csv", excluded)
    dump(out / "H5AD_schema_records.json", schemas)
    return (
        blocks,
        manifests,
        dict(
            scope=f"All {len(inventory)} DFC donor objects in {RELEASE} release; all source subclasses",
            phenotype="Clinical status and Braak read directly from official donor metadata",
            model="Model definitions are supplied by analysis/models.R",
            limitations="Post-QC recovered nuclei, not tissue composition. No EC, spatial mechanism or independent donor-overlap certification.",
        ),
    )


def combine(blocks, rows, out, kind):
    """合并各供体的基因集合与计数列。

    不同对象缺少的基因记录为缺失，不冒充零表达；输出 expression.csv.gz 和
    metadata.csv，并检查表达列与单位记录一一对应。"""
    need(blocks and rows, "No analyzable units prepared")
    allgenes = sorted(set().union(*(set(b[0]) for b in blocks)))
    lookup = {g: i for i, g in enumerate(allgenes)}
    ids = [x for _, _, cols in blocks for x in cols]
    need(ids == [r["unit_id"] for r in rows] and len(ids) == len(set(ids)), "Unit mismatch")
    data = np.full((len(allgenes), len(ids)), np.nan)
    start = 0
    for genes, values, cols in blocks:
        ix = np.array([lookup[g] for g in genes])
        data[np.ix_(ix, np.arange(start, start + len(cols)))] = values
        start += len(cols)
    pd.DataFrame(data, index=pd.Index(allgenes, name="gene_symbol"), columns=ids).to_csv(
        out / "expression.csv.gz", compression="gzip"
    )
    write_csv(out / "metadata.csv", rows)
