# -*- coding: utf-8 -*-
"""s64d_rename.py —— 裁定 6 的档务重命名（★纯重命名，被重命名文件的内容零改动）

范围（仅「同前缀 + 同暂停点标签」的硬冲突组，以及 28 组的三重前缀）：
  28b_任务A_三位点比较_暂停点12.md            -> 28b_任务A_三位点比较_暂停点12.md
  28c_正文框架_审稿意见修正版.md              -> 28c_正文框架_审稿意见修正版.md
  31b_裁定1GTEx补充与任务B_MVMR_暂停点14.md   -> 31b_裁定1GTEx补充与任务B_MVMR_暂停点14.md
保留为规范名（引用最多 / 最贴合暂停点主交付）：
  28_论文正文初稿与图注定稿_暂停点12.md、31_任务三3.2条件PheWAS_暂停点14.md

不做的事：
  · 不改任何被重命名文件的**内容**（sha256 前后必须一致）
  · 不改项目内任何**报告 .md/.docx** 的内容（悬空引用只**报告**，不擅改，守「只做重命名不改内容」）
  · 只更新**我自己的**记忆文件与脚本中的指针
  · 其他同前缀组（29/36/37/38/39/40/41/42/43/39b）无暂停点标签，属「同一产物多格式/多版本家族」，**有意保留**

备份：D:/_transfer_logs/_rename_backup_20260929/
日志：logs/s64d_rename.log
"""
import io, os, re, shutil, hashlib

ROOT = r"D:/endometriosis_project/11_sc_eqtl_mr_project"
MEM = r"D:/endometriosis_project/.workbuddy/memory"
BAK = r"D:/_transfer_logs/_rename_backup_20260929"
LOG = os.path.join(ROOT, "logs", "s64d_rename.log")

RENAMES = [
    ("28b_任务A_三位点比较_暂停点12.md", "28b_任务A_三位点比较_暂停点12.md"),
    ("28c_正文框架_审稿意见修正版.md", "28c_正文框架_审稿意见修正版.md"),
    ("31b_裁定1GTEx补充与任务B_MVMR_暂停点14.md", "31b_裁定1GTEx补充与任务B_MVMR_暂停点14.md"),
]
SELF_SCRIPTS = ["s64a_recon.py", "s64c_refscan.py", "s64d_rename.py"]

L = []
def p(*a):
    s = " ".join(str(x) for x in a); L.append(s)
    try: print(s)
    except UnicodeEncodeError: print(s.encode("ascii", "backslashreplace").decode("ascii"))


def sha(pth):
    return hashlib.sha256(open(pth, "rb").read()).hexdigest()


def main():
    import atexit
    atexit.register(lambda: io.open(LOG, "w", encoding="utf-8", newline="\n").write("\n".join(L) + "\n"))
    os.makedirs(BAK, exist_ok=True)
    p("== 档务重命名（纯重命名）==")
    p("备份目录 %s" % BAK)

    # 0) 前置校验
    for old, new in RENAMES:
        op = os.path.join(ROOT, old); np_ = os.path.join(ROOT, new)
        assert os.path.exists(op), "源不存在：%s" % old
        assert not os.path.exists(np_), "目标已存在，拒绝覆盖：%s" % new

    # 1) 备份 + 重命名
    p("-- 重命名 --")
    for old, new in RENAMES:
        op = os.path.join(ROOT, old); np_ = os.path.join(ROOT, new)
        h_before, sz = sha(op), os.path.getsize(op)
        shutil.copy2(op, os.path.join(BAK, old))
        os.rename(op, np_)
        h_after = sha(np_)
        assert h_before == h_after, "★内容被改动：%s" % old
        assert not os.path.exists(op)
        p("  %s" % old)
        p("      -> %s" % new)
        p("      size=%d B  sha256=%s（前后一致，内容零改动）" % (sz, h_after[:16]))

    # 2) 更新我自己的记忆文件中的指针
    p("-- 更新记忆文件指针（我自己的 bookkeeping）--")
    n_mem = 0
    for f in sorted(os.listdir(MEM)):
        if not f.endswith(".md"):
            continue
        fp = os.path.join(MEM, f)
        txt = io.open(fp, encoding="utf-8").read()
        orig = txt
        for old, new in RENAMES:
            o_stem, n_stem = old[:-3], new[:-3]
            txt = txt.replace(old, new).replace(o_stem, n_stem)
        if txt != orig:
            io.open(fp, "w", encoding="utf-8", newline="\n").write(txt)
            n_mem += 1
            p("  updated %s" % f)
    p("  记忆文件更新数 = %d" % n_mem)

    # 2b) 更新我自己的侦察脚本中的字面量
    p("-- 更新本次侦察脚本中的字面量 --")
    for s in SELF_SCRIPTS:
        fp = os.path.join(ROOT, "scripts", s)
        if not os.path.exists(fp):
            continue
        txt = io.open(fp, encoding="utf-8").read(); orig = txt
        for old, new in RENAMES:
            txt = txt.replace(old, new)
        if txt != orig:
            io.open(fp, "w", encoding="utf-8", newline="\n").write(txt)
            p("  updated scripts/%s" % s)

    # 3) 扫描：项目内（非记忆）仍指向旧名的悬空引用（只报告，不修改）
    p("-- 项目内悬空引用扫描（只报告，不擅改）--")
    SKIP = {"00_data_raw", "02_data_processed", "figures", "_onek1k_plink",
            "_chr1_ld", "_panel", "_tools", "__pycache__"}
    stale = []
    for dp, dn, fn in os.walk(ROOT):
        dn[:] = [d for d in dn if d not in SKIP and not d.startswith("_fig")]
        for f in fn:
            if not f.endswith((".md", ".txt")):
                continue
            if os.path.basename(dp) == "scripts":
                continue
            fp = os.path.join(dp, f)
            if os.path.basename(fp) in [o for o, _ in RENAMES]:
                continue
            try:
                txt = io.open(fp, encoding="utf-8", errors="replace").read()
            except Exception:
                continue
            for old, new in RENAMES:
                c = txt.count(old)
                if c:
                    stale.append((os.path.relpath(fp, r"D:/endometriosis_project"), old, c))
    if not stale:
        p("  无悬空引用")
    else:
        for rel, old, c in stale:
            p("  STALE  %s  ×%d  -> %s" % (rel, c, old))
    p("  悬空引用条目数 = %d" % len(stale))

    p("")
    p("RENAME_DONE（被重命名文件内容 sha256 前后一致；未改动任何报告 .md）")


if __name__ == "__main__":
    main()
