# -*- coding: utf-8 -*-
"""
s70 · 把 63_figure_inventory.csv 中 5 张 D 系列诊断图的「状态」由笼统的
`探索性诊断` 改写为「已被主图/补充图承接」的分级措辞。

依据（2026-09-29）：
  - 用户的提问「诊断图是不需要在投稿系统上传吗？也不用在正文中提及？」
    需要 inventory 能自证「D 系列已 superseded、不上传、不引用」，
    而不是只写一个含义不明的 `探索性诊断`。
  - 技能 publication-figure-export-qa 已确立的分级措辞：
    `superseded by Fig. X<panel>; retained for provenance`。

映射（逐图核过图注，非推测）：
  D1 volcano            -> Fig. S4        （差异：D1 含 MHC 着色，S4 依正文口径去除）
  D2 Manhattan 14 小倍数 -> Fig. S3        （14 面板逐细胞类型 Manhattan，同构）
  D3 forest 15 对        -> Fig. 2b        （15 对森林图；Fig. 2b 另加 IVW 菱形）
  D4 细胞类型负荷        -> Fig. S1(b–d)   （名义数/FDR 数/最小 q 三指标逐项对应）
  D5 chr1 LD 结构        -> Fig. 4a        （7 变异两两 r² 热图；Fig. 4b 另给基因轨道）

硬纪律：写文件用 with io.open(..., newline='') 保 LF；改后必须 getsize + 计数复核。
"""
import io
import os
import shutil
import sys

ROOT = r"D:/endometriosis_project/11_sc_eqtl_mr_project"
CSV = os.path.join(ROOT, "tables", "63_figure_inventory.csv")
BACKUP_DIR = r"D:/_transfer_logs/_ms_backup_20260929_s1s4"
LOG = os.path.join(ROOT, "logs", "_s70_patch_inventory.log")

NEW = {
    "FigD1": "superseded by Fig. S4；仅留存溯源，不投稿、不引用",
    "FigD2": "superseded by Fig. S3；仅留存溯源，不投稿、不引用",
    "FigD3": "superseded by Fig. 2b；仅留存溯源，不投稿、不引用",
    "FigD4": "superseded by Fig. S1(b-d)；仅留存溯源，不投稿、不引用",
    "FigD5": "superseded by Fig. 4a；仅留存溯源，不投稿、不引用",
}


def main():
    out = []
    with open(CSV, "rb") as fh:
        raw = fh.read()
    nl = b"\r\n" if (raw.count(b"\r\n") > 0 and raw.count(b"\r\n") == raw.count(b"\n")) else b"\n"
    out.append("file bytes before = %d" % len(raw))
    out.append("newline = %s" % ("CRLF" if nl == b"\r\n" else "LF"))
    out.append("old token count = %d" % raw.count("探索性诊断".encode("utf-8")))

    text = raw.decode("utf-8-sig")
    lines = text.split("\n")
    out.append("lines = %d" % len(lines))
    out.append("header = %r" % lines[0])

    patched = 0
    for i, ln in enumerate(lines):
        if not ln.startswith("诊断图,"):
            continue
        cells = ln.split(",")
        tag = cells[1]
        if tag not in NEW:
            out.append("  !! unexpected diagnostic row: %r" % ln[:80])
            continue
        if cells[-1] != "探索性诊断":
            out.append("  -- %s already patched? last=%r" % (tag, cells[-1]))
            continue
        cells[-1] = NEW[tag]
        lines[i] = ",".join(cells)
        patched += 1
        out.append("  ++ %s : %s" % (tag, NEW[tag]))

    out.append("patched rows = %d" % patched)

    if patched != 5:
        out.append("ABORT: expected 5 rows, patched %d" % patched)
        write_log(out)
        sys.exit(2)

    new_text = "\n".join(lines)
    new_raw = new_text.encode("utf-8")
    # 还原 BOM
    if raw.startswith(b"\xef\xbb\xbf"):
        new_raw = b"\xef\xbb\xbf" + new_raw

    # 先备份
    os.makedirs(BACKUP_DIR, exist_ok=True)
    shutil.copy2(CSV, os.path.join(BACKUP_DIR, "63_figure_inventory.csv.bak"))
    with io.open(CSV, "wb") as fh:
        fh.write(new_raw)

    # 复核
    chk = open(CSV, "rb").read()
    out.append("file bytes after = %d (delta %+d)" % (len(chk), len(chk) - len(raw)))
    out.append("old token count after = %d" % chk.count("探索性诊断".encode("utf-8")))
    for tag, s in NEW.items():
        out.append("new[%s] count = %d" % (tag, chk.count(s.encode("utf-8"))))
    out.append("BOM preserved = %s" % chk.startswith(b"\xef\xbb\xbf"))
    out.append("CRLF in new file = %d" % chk.count(b"\r\n"))
    out.append("VERDICT = %s" % ("PASS" if (chk.count("探索性诊断".encode("utf-8")) == 0
                                            and all(chk.count(s.encode("utf-8")) == 1 for s in NEW.values()))
                                 else "FAIL"))
    write_log(out)
    print("done rc=0")


def write_log(lines):
    os.makedirs(os.path.dirname(LOG), exist_ok=True)
    with io.open(LOG, "w", encoding="utf-8", newline="\n") as fh:
        fh.write("\n".join(lines) + "\n")


if __name__ == "__main__":
    main()
