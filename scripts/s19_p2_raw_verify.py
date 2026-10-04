#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""s19_p2_raw_verify.py —— P2 OneK1K raw tensorQTL summary tar.gz 格式核验（暂停点 2）

流式读取，不落盘解压。产出 JSON 证据文件。
阶段：
  A) gzip magic + tarfile 可开 + 成员清单（前 40）+ 成员总数（限时遍历）
  B) 抽样成员首行表头/首行数据
"""
import argparse, json, os, sys, time, io

import tarfile


def natural_chr_key(s):
    try:
        return (0, int(s))
    except Exception:
        return (1, 0)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--tar', required=True)
    ap.add_argument('--out', required=True)
    ap.add_argument('--time-limit', type=float, default=480.0)
    ap.add_argument('--max-names', type=int, default=40)
    ap.add_argument('--n-samples', type=int, default=3)
    ap.add_argument('--sample-lines', type=int, default=3)
    ap.add_argument('--dump-names', default=None, help='把所有成员名/大小写入该 TSV')
    ap.add_argument('--no-early-exit', action='store_true', help='禁用提前收工，穷尽遍历')
    a = ap.parse_args()

    dump = None
    if a.dump_names:
        dump = open(a.dump_names, 'w', encoding='utf-8', newline='')
        dump.write('idx\tname\tsize\tisdir\n')

    R = {'tar': a.tar, 'goal_bytes': 10344009571}
    t0 = time.time()

    # ---------- 阶段 0：gzip magic ----------
    with open(a.tar, 'rb') as f:
        head = f.read(4)
    R['size_bytes_handle'] = os.path.getsize(a.tar)
    R['head_hex'] = head.hex()
    R['gzip_magic_ok'] = (head[:2] == b'\x1f\x8b')
    R['gzip_cm_flg'] = {'cm': head[2], 'flg': head[3]} if len(head) >= 4 else None
    R['mtime_field'] = int.from_bytes(head[4:8], 'little') if len(head) >= 8 else None
    R['mtime_iso'] = time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(R['mtime_field'])) if R.get('mtime_field') else None

    # ---------- 阶段 A/B：流式遍历 ----------
    names = []
    samples = {}
    n_members = 0
    n_dir = 0
    n_file = 0
    decomp_bytes = 0
    truncated = False
    sizes_top = []
    err = None
    t_open = time.time()
    try:
        tf = tarfile.open(a.tar, mode='r|gz')
        R['tarfile_open_ok'] = True
        R['tar_format'] = getattr(tf, 'format', None)
        for m in tf:
            n_members += 1
            if dump is not None:
                dump.write('%d\t%s\t%d\t%d\n' % (n_members, m.name, m.size, 1 if m.isdir() else 0))
            if m.isdir():
                n_dir += 1
            else:
                n_file += 1
                decomp_bytes += max(m.size, 0)
            if len(names) < a.max_names:
                names.append({'i': n_members, 'name': m.name, 'size': m.size,
                              'isdir': bool(m.isdir())})
            if len(sizes_top) < 5 and m.isfile():
                sizes_top.append([m.name, m.size])
            if m.isfile() and len(samples) < a.n_samples and m.size > 0:
                fh = tf.extractfile(m)
                lines = []
                for _ in range(a.sample_lines):
                    ln = fh.readline()
                    if not ln:
                        break
                    lines.append(ln.decode('utf-8', 'replace').rstrip('\r\n'))
                try:
                    fh.close()
                except Exception:
                    pass
                samples[m.name] = {'size': m.size, 'head_lines': lines}
            elif (not a.no_early_exit) and m.isfile() and len(samples) >= a.n_samples and len(names) >= a.max_names and n_members > 200:
                # 已采集够信息且成员数已远超 40，可提前收工（仍受限时保护）
                if time.time() - t0 > 120:
                    truncated = True
                    break
            if time.time() - t0 > a.time_limit:
                truncated = True
                break
        try:
            tf.close()
        except Exception:
            pass
    except Exception as e:
        R['tarfile_open_ok'] = R.get('tarfile_open_ok', False)
        err = repr(e)
    if dump is not None:
        dump.close()
    R['traverse_seconds'] = round(time.time() - t_open, 2)
    R['error'] = err
    R['n_members_seen'] = n_members
    R['n_dirs_seen'] = n_dir
    R['n_files_seen'] = n_file
    R['decomp_bytes_seen'] = decomp_bytes
    R['traverse_truncated'] = truncated
    R['first_members'] = names
    R['first_file_sizes'] = sizes_top
    R['samples'] = samples

    # 成员名规律推断
    exts = {}
    for x in names:
        if x['isdir']:
            continue
        b = os.path.basename(x['name'])
        if '.' in b:
            e = b.rsplit('.', 1)[1]
        else:
            e = '<none>'
        exts[e] = exts.get(e, 0) + 1
    R['ext_hist_first40'] = exts
    # 目录第一段
    tops = {}
    for x in names:
        p = x['name'].split('/')
        if len(p) > 1:
            tops[p[0]] = tops.get(p[0], 0) + 1
    R['top_dirs_first40'] = tops

    with open(a.out, 'w', encoding='utf-8') as f:
        json.dump(R, f, ensure_ascii=False, indent=2)
    print('WROTE ' + a.out)
    print('members_seen=%d files=%d dirs=%d truncated=%s secs=%.1f' %
          (n_members, n_file, n_dir, truncated, R['traverse_seconds']))


if __name__ == '__main__':
    main()
