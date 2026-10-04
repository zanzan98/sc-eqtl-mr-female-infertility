# -*- coding: utf-8 -*-
"""s75a · 给 s37_final_figures.py 增加 FIG_NOTES 开关（默认 =1，保持原行为）

做法：4 处 fig.text(0.5, ...) 脚注块前置 if SHOW_NOTES: 并整体缩进 4 空格。
不动 fig.subplots_adjust（保留脚注曾占的底部空间 -> tight 裁剪后自动变矮，
面板的绝对位置与尺寸不变）。
"""
import io, os, sys, hashlib

P = r"D:/endometriosis_project/11_sc_eqtl_mr_project/scripts/s37_final_figures.py"
BAK = r"D:/_transfer_logs/_figs_backup_20260929_nonotes/s37_final_figures.py.pre_switch"
LOG = r"D:/endometriosis_project/11_sc_eqtl_mr_project/logs/_s75a_add_notes_switch.txt"

out = []
def w(s=""):
    out.append(str(s)); print(s)

raw = open(P, "rb").read()
w("before bytes=%d md5=%s" % (len(raw), hashlib.md5(raw).hexdigest()))
txt = raw.decode("utf-8")
lines = txt.splitlines(keepends=True)

# 幂等守卫
if "SHOW_NOTES" in txt:
    w("!! 已含 SHOW_NOTES，判定为已打过补丁 -> 中止（幂等）")
    with io.open(LOG, "w", encoding="utf-8", newline="\n") as fh:
        fh.write("\n".join(out) + "\n")
    sys.exit(0)

ENDS = 'color="#555555")'
starts = [i for i, l in enumerate(lines) if l.lstrip().startswith("fig.text(0.5,")]
w("找到 fig.text(0.5, ...) 起始行 %d 处：%s" % (len(starts), [i + 1 for i in starts]))
assert len(starts) == 4, "预期 4 处脚注，实测 %d" % len(starts)

# 定位每处结束行
spans = []
for i in starts:
    j = None
    for k in range(i, min(i + 40, len(lines))):
        if ENDS in lines[k]:
            j = k
            break
    assert j is not None, "第 %d 行起的脚注块未找到结束行" % (i + 1)
    spans.append((i, j))
for i, j in spans:
    w("  块 %d..%d  (缩进 %d 空格)  首行: %s" % (i + 1, j + 1, len(lines[i]) - len(lines[i].lstrip()), lines[i].strip()[:70]))

# 从后往前：整块 +4 缩进，并在块前插入 if
for i, j in reversed(spans):
    ind = len(lines[i]) - len(lines[i].lstrip())
    assert ind == 4, "缩进异常 %d" % ind
    for k in range(i, j + 1):
        if lines[k].strip():
            lines[k] = "    " + lines[k]
    lines.insert(i, "    if SHOW_NOTES:\n")

txt2 = "".join(lines)

# 插入开关定义
ANCHOR = "from matplotlib.patches import FancyBboxPatch, FancyArrowPatch, Rectangle"
assert txt2.count(ANCHOR) == 1
NEW = ANCHOR + (
    "\n\n# ★图内脚注（panel notes）开关：FIG_NOTES=0 时不画底部说明文字（2026-09-29 用户裁定）\n"
    'SHOW_NOTES = os.environ.get("FIG_NOTES", "1") != "0"'
)
txt2 = txt2.replace(ANCHOR, NEW, 1)

with io.open(P, "w", encoding="utf-8", newline="") as fh:
    fh.write(txt2)

raw2 = open(P, "rb").read()
w("")
w("after  bytes=%d md5=%s" % (len(raw2), hashlib.md5(raw2).hexdigest()))
w("SHOW_NOTES 出现次数 = %d" % raw2.decode("utf-8").count("SHOW_NOTES"))
w("if SHOW_NOTES: 出现次数 = %d" % raw2.decode("utf-8").count("if SHOW_NOTES:"))
w("fig.text(0.5, 仍在 = %d" % raw2.decode("utf-8").count("fig.text(0.5,"))

# 语法编译检查
import py_compile
try:
    py_compile.compile(P, cfile=os.path.join(os.environ.get("TEMP", "."), "_s37_check.pyc"), doraise=True)
    w("py_compile = OK")
    OK = True
except Exception as e:
    w("py_compile = FAILED: %s" % e)
    OK = False

with io.open(BAK, "w", encoding="utf-8", newline="") as fh:
    fh.write(raw.decode("utf-8"))
w("原始副本 -> %s" % BAK)
w("")
w("VERDICT = %s" % ("PASS" if OK and raw2.decode("utf-8").count("if SHOW_NOTES:") == 4 else "FAIL"))

with io.open(LOG, "w", encoding="utf-8", newline="\n") as fh:
    fh.write("\n".join(out) + "\n")
