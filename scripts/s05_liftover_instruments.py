#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""s05_liftover_instruments.py —— 通用工具变量位置 liftover GRCh37 -> GRCh38

输入表须含 variant_id 列（格式 chr:pos，GRCh37）。
输出：variant_id 被替换为 GRCh38；新增 variant_id_grch37 / liftover_status / pos_grch37 / pos_grch38。

★ 坐标系背景（见 00b_执行期更正_2026-09-22.md §1）：
  OneK1K bim 与 raw/top 汇总的位置均为 GRCh37；
  GWAS Catalog harmonised 结局的 base_pair_location 为 GRCh38。
  两者直接按位置连接会系统性错配。
"""
import argparse, os, sys, traceback
import pandas as pd
from pyliftover import LiftOver


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--in', dest='inp', required=True)
    ap.add_argument('--out', dest='outp', required=True)
    ap.add_argument('--log', default=None)
    a = ap.parse_args()

    L = []
    try:
        d = pd.read_csv(a.inp, dtype={'variant_id': str})
        if 'variant_id_grch37' in d.columns:
            d['variant_id_grch37'] = d['variant_id_grch37'].astype(str)
            d['variant_id'] = d['variant_id_grch37']
        L.append('输入: %d 行 / %d 基因 / %d 细胞类型' % (len(d), d.gene.nunique(), d.cell_type.nunique()))
        ch = d['variant_id'].str.split(':').str[0]
        ps = d['variant_id'].str.split(':').str[1].astype(int)
        d['pos_grch37'] = ps
        # ★ 勿用 ch.min()/ch.max()：对 object 型字符串是字典序比较（'9' > '22'），会给出误导性区间
        _uks = sorted(set(ch.dropna().astype(str)), key=lambda x: (0, int(x)) if x.isdigit() else (1, 0))
        L.append('GRCh37 位置域: 染色体数=%d %s ; pos %d..%d' % (len(_uks), ','.join(_uks), ps.min(), ps.max()))

        lo = LiftOver('hg19', 'hg38')
        uniq = pd.DataFrame({'c': ch, 'p': ps}).drop_duplicates()
        L.append('唯一位点 = %d' % len(uniq))
        mp = {}
        nok = nmul = nmis = 0
        for c, p in zip(uniq['c'], uniq['p']):
            r = lo.convert_coordinate('chr' + str(c), int(p) - 1)
            if not r:
                mp[(c, p)] = (None, None, 'unmapped'); nmis += 1
            elif len(r) > 1:
                mp[(c, p)] = (r[0][0].replace('chr', ''), int(r[0][1]) + 1, 'multi'); nmul += 1
            else:
                mp[(c, p)] = (r[0][0].replace('chr', ''), int(r[0][1]) + 1, 'ok'); nok += 1
        L.append('唯一映射 = %d ; 多重 = %d ; 无映射 = %d' % (nok, nmul, nmis))

        d['chr_grch38'] = [mp[(c, p)][0] for c, p in zip(ch, ps)]
        d['pos_grch38'] = [mp[(c, p)][1] for c, p in zip(ch, ps)]
        d['liftover_status'] = [mp[(c, p)][2] for c, p in zip(ch, ps)]
        d['variant_id_grch37'] = d['variant_id']
        okm = d['pos_grch38'] == d['pos_grch38']
        d.loc[okm, 'variant_id'] = (d.loc[okm, 'chr_grch38'].astype(str) + ':' +
                                    d.loc[okm, 'pos_grch38'].astype(int).astype(str))
        d.loc[~okm, 'variant_id'] = None

        ok = d[d['liftover_status'] == 'ok'].copy()
        ok['shift'] = ok['pos_grch38'] - ok['pos_grch37']
        L.append('位移 min/median/max = %d / %.0f / %d bp' % (ok['shift'].min(), ok['shift'].median(), ok['shift'].max()))
        L.append('位移==0 比例 = %.1f%%' % (100.0 * (ok['shift'] == 0).mean()))
        L.append('|位移|>1Mb 比例 = %.2f%%' % (100.0 * (ok['shift'].abs() > 1e6).mean()))
        L.append('状态计数: %s' % d['liftover_status'].value_counts().to_dict())

        d.to_csv(a.outp, index=False, encoding='utf-8-sig')
        L.append('已写出: %s (%d 行)' % (a.outp, len(d)))
    except Exception:
        L.append('!!! EXCEPTION !!!'); L.append(traceback.format_exc())
    txt = '\n'.join(L)
    print(txt)
    if a.log:
        open(a.log, 'w', encoding='utf-8').write(txt)


if __name__ == '__main__':
    main()
