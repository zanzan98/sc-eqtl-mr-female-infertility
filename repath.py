#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""repath.py -- 把本仓库内脚本的硬编码绝对路径重定位到你的机器。

为什么需要
----------
本仓库的分析脚本由作者在 Windows 本机上编写，脚本顶部直接写了绝对路径。
脚本内容与产出结果的脚本**逐字节一致**（保证可复现性），因此仓库内不预先改写，
而是提供本工具：克隆后运行一次即可把两处根目录前缀整体换到你机器的位置。

需要替换的两处根目录
--------------------
  A. 工程根  --project
     原值：D:/endometriosis_project/11_sc_eqtl_mr_project
  B. 外置数据根  --external
     原值：D:/endometriosis_project
     其下含 _panel / _onek1k_plink / _chr1_ld / _tools（见 docs/DATA_SOURCES.md）

用法
----
    # 1) 先看会发生什么（默认 dry-run，不写盘）
    python repath.py --project "D:/my/11_sc_eqtl_mr_project" --external "D:/my/endometriosis_project"

    # 2) 确认后落盘（自动备份）
    python repath.py --project "D:/my/11_sc_eqtl_mr_project" --external "D:/my/endometriosis_project" --apply

    # 3) 连文档一起改（可选）
    python repath.py ... --apply --include-docs

安全设计
--------
1. **字节级替换**：按 bytes 读写，不编解码、不改变换行符 —— 除被替换的前缀外逐字节不变。
2. **先长后短**：先替换较长的「工程根」，再替换较短的「外置根」，避免后者误伤前者。
3. **整字面量归一**：命中旧前缀的字符串字面量，其**全部分隔符统一为 `/`**
   （避免 `X:/a/b\\c\\d` 这类混合分隔符在非 Windows 上失效）。
4. **两遍处理**：先按字符串字面量重写（保引号与 `r`/`f` 前缀不动），
   再对注释、docstring 等非单行字面量场合做残余前缀替换。
5. **写盘前自证**（任一失败即整体中止，不写任何文件）：
   (a) 逐文件字节增量 == 已记录替换的实测增量之和（证明没有越界改动）；
   (b) 输出中不得再出现任何形态的旧根前缀。
6. `--apply` 时先备份到 `_repath_backup_<时间戳>/`（保留相对路径），再写入。
7. 结束后列出**仍未处理**的绝对路径及出处文件，供人工确认（如 C:/Users/... 之类）。
"""

from __future__ import annotations

import argparse
import re
import sys
import time
from pathlib import Path

REPO = Path(__file__).resolve().parent
SCRIPTS_DIR = REPO / "scripts"

#: 原始根目录（作者本机）。仅作匹配用，勿改。
OLD_PROJECT = "D:/endometriosis_project/11_sc_eqtl_mr_project"
OLD_EXTERNAL = "D:/endometriosis_project"

#: 默认处理范围
SCRIPT_GLOBS = ("**/*.py", "**/*.R")
DOC_GLOBS = ("docs/**/*.md", "*.md", "*.txt", "*.yml", "*.cff")

_SEP = rb"[\\/]+"
# 字符串字面量（同一行内，不含跨行）：可选 r/b/u/f 前缀 + 成对引号
_RX_LIT = re.compile(rb"""(?P<pre>[rRbBuUfF]{0,3})(?P<q>["'])(?P<body>[^"'\n]*)(?P=q)""")
# 绝对路径探测：盘符不得紧跟在字母/数字/下划线之后（排除 "MODEL:\n" 这类转义假阳性），
# 且首段至少 2 字符或以 '_' 开头（排除 \n \t \\ 等单字母转义）。
_RX_ABS = re.compile(
    rb"(?<![A-Za-z0-9_])[A-Za-z]:[\\/]{1,2}(?:[A-Za-z0-9_.\-]{2,}|_[A-Za-z0-9_.\-]+)[^\s\"'`)\]]*"
)


def _prefix_rx(fwd_path: str) -> re.Pattern:
    """把正斜杠形式的根目录编译成「容忍 / 与 \\\\ 混用」的字节正则。"""
    parts = fwd_path.replace("\\", "/").strip("/").split("/")
    body = _SEP.join(re.escape(p).encode("ascii") for p in parts)
    return re.compile(body, re.IGNORECASE)


def _norm_root(p: str) -> str:
    return p.replace("\\", "/").strip().rstrip("/")


def _norm(b: bytes) -> bytes:
    return b.replace(b"\\", b"/")


def _rewrite_literals(text: bytes, old_p: bytes, old_e: bytes | None,
                      new_p: bytes, new_e: bytes | None) -> tuple[bytes, list[int], int]:
    """第一遍：命中旧前缀的字符串字面量，整段分隔符归一为 '/'。

    返回 (新文本, [规则1命中, 规则2命中], 实测字节增量)。
    """
    hits = [0, 0]
    delta = 0

    def repl(m: re.Match) -> bytes:
        nonlocal delta
        pre, q, body = m.group("pre"), m.group("q"), m.group("body")
        norm = _norm(body)
        if old_p in norm:
            norm = norm.replace(old_p, new_p, 1)
            hits[0] += 1
        elif old_e is not None and new_e is not None and old_e in norm:
            norm = norm.replace(old_e, new_e, 1)
            hits[1] += 1
        else:
            return m.group(0)
        out = pre + q + norm + q
        delta += len(out) - len(m.group(0))
        return out

    return _RX_LIT.sub(repl, text), hits, delta


def _rewrite_residual(text: bytes, rules: list[tuple[re.Pattern, bytes]]) -> tuple[bytes, list[int], int]:
    """第二遍：注释 / docstring 等非单行字面量场合的残余前缀替换。

    同时把紧随新根之后的多余分隔符归一为单个 '/'。返回 (新文本, 命中数, 实测增量)。
    """
    hits: list[int] = []
    delta = 0
    for rx, new_b in rules:
        n = 0

        def repl(m: re.Match, _new=new_b) -> bytes:
            nonlocal n, delta
            n += 1
            delta += len(_new) - len(m.group(0))
            return _new

        text = rx.sub(repl, text)
        hits.append(n)

        def collapse(m: re.Match, _new=new_b) -> bytes:
            nonlocal delta
            out = _new + b"/"
            delta += len(out) - len(m.group(0))
            return out

        text = re.sub(re.escape(new_b) + _SEP, collapse, text)
    while len(hits) < 2:      # 与 _rewrite_literals 的 [规则1, 规则2] 对齐
        hits.append(0)
    return text, hits, delta


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(
        description="重定位本仓库脚本内的硬编码绝对路径（默认 dry-run）",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    ap.add_argument("--project", required=True,
                    help="新的工程根，例如 D:/my/11_sc_eqtl_mr_project")
    ap.add_argument("--external",
                    help="新的外置数据根，例如 D:/my/endometriosis_project（省略则跳过该规则）")
    ap.add_argument("--apply", action="store_true",
                    help="真正写盘（默认只报告）")
    ap.add_argument("--include-docs", action="store_true",
                    help="连同 docs/ 与仓库根级 *.md/*.txt/*.yml/*.cff 一起处理")
    ap.add_argument("--scripts-dir", default=None,
                    help="覆盖 scripts/ 目录（自测用）")
    args = ap.parse_args(argv)

    scripts_dir = Path(args.scripts_dir) if args.scripts_dir else SCRIPTS_DIR
    if not scripts_dir.is_dir():
        print(f"[ERR] 找不到脚本目录：{scripts_dir}", file=sys.stderr)
        return 2

    new_project = _norm_root(args.project)
    new_external = _norm_root(args.external) if args.external else None

    old_p_b, old_e_b = OLD_PROJECT.encode(), OLD_EXTERNAL.encode()
    new_p_b = new_project.encode()
    new_e_b = new_external.encode() if new_external else None

    rules: list[tuple[re.Pattern, bytes]] = [(_prefix_rx(OLD_PROJECT), new_p_b)]
    if new_external:
        rules.append((_prefix_rx(OLD_EXTERNAL), new_e_b))
    # 断言 (b) 用：旧根在任意分隔符形态下的探测
    old_any = re.compile(b"|".join(rx.pattern for rx, _ in rules), re.IGNORECASE)

    # ---- 收集目标文件 ----
    files: list[Path] = []
    for g in SCRIPT_GLOBS:
        files.extend(sorted(scripts_dir.glob(g)))
    if args.include_docs:
        for g in DOC_GLOBS:
            files.extend(sorted(REPO.glob(g)))
    files = sorted({f for f in files if f.is_file()})

    print(f"[i] 仓库根：{REPO}")
    print(f"[i] 规则 1（工程根）：{OLD_PROJECT}  ->  {new_project}")
    print(f"[i] 规则 2（外置根）：{OLD_EXTERNAL}  ->  {new_external}" if new_external
          else "[i] 规则 2（外置根）：已跳过（未提供 --external）")
    print(f"[i] 目标文件数：{len(files)}｜模式：{'APPLY' if args.apply else 'DRY-RUN'}")
    print("-" * 72)

    # ================= 阶段 A：全部计算 + 校验（不写盘） =================
    planned: list[tuple[Path, bytes, bytes, list[int]]] = []
    total_hits = [0, 0]
    total_delta = 0
    remaining: dict[str, set[str]] = {}

    for f in files:
        raw = f.read_bytes()
        new, h_lit, d_lit = _rewrite_literals(raw, old_p_b, old_e_b, new_p_b, new_e_b)
        new, h_res, d_res = _rewrite_residual(new, rules)
        per_rule = [h_lit[i] + h_res[i] for i in range(2)]
        for i in range(2):
            total_hits[i] += per_rule[i]

        rel = str(f.relative_to(REPO)) if f.is_relative_to(REPO) else str(f)

        # (a) 字节增量自证：整体增量必须等于已记录替换的实测增量之和
        if len(new) - len(raw) != d_lit + d_res:
            print(f"[FATAL] {rel}: 字节增量自证失败 —— 整体 {len(new) - len(raw)}，"
                  f"替换累计 {d_lit + d_res}。已中止，未写任何文件。", file=sys.stderr)
            return 3
        # (b) 旧根残留自证（先屏蔽新根，避免「新根合法包含旧根串」时误报）
        masked = new
        for nb in (new_p_b, new_e_b):
            if nb:
                masked = masked.replace(nb, b"\x00" * len(nb))
        if old_any.search(masked):
            print(f"[FATAL] {rel}: 替换后仍存在旧根前缀，已中止，未写任何文件。", file=sys.stderr)
            return 3

        total_delta += len(new) - len(raw)
        if new != raw:
            planned.append((f, raw, new, per_rule))

        for m in _RX_ABS.finditer(new):
            s = m.group(0).decode("ascii", "replace")
            if s.startswith(new_project) or (new_external and s.startswith(new_external)):
                continue
            if "\ufffd" in s:
                s += "  [含非 ASCII，请人工确认]"
            remaining.setdefault(s, set()).add(rel)

    # ================= 阶段 B：落盘（含备份） =================
    backup_dir = None
    if args.apply and planned:
        backup_dir = REPO / f"_repath_backup_{time.strftime('%Y%m%d_%H%M%S')}"
        backup_dir.mkdir(parents=True, exist_ok=True)
        for f, raw, new, _ in planned:
            dst = backup_dir / f.relative_to(REPO)
            dst.parent.mkdir(parents=True, exist_ok=True)
            dst.write_bytes(raw)
            f.write_bytes(new)

    for f, _, _, hits in planned:
        rel = str(f.relative_to(REPO)) if f.is_relative_to(REPO) else str(f)
        print(f"  {rel}   命中 = {hits}")
    print("-" * 72)
    print(f"[i] 变更文件：{len(planned)} / {len(files)}")
    print(f"[i] 各规则命中合计：{total_hits}（顺序 = 工程根 -> 外置根）")
    print(f"[i] 总字节增量：{total_delta}（全部文件自证通过）")
    if args.apply and backup_dir is not None:
        print(f"[i] 备份目录：{backup_dir}")

    if remaining:
        print(f"[i] 仍未处理的绝对路径（{len(remaining)} 个，请人工确认）：")
        for i, (s, where) in enumerate(sorted(remaining.items())):
            if i >= 25:
                print(f"      ...（另有 {len(remaining) - 25} 个）")
                break
            w = ", ".join(sorted(where)[:2])
            more = f" 等 {len(where)} 处" if len(where) > 2 else ""
            print(f"      {s}   <- {w}{more}")
    else:
        print("[i] 未发现其他绝对路径。")

    if not args.apply:
        print("[i] 以上为 dry-run，未写盘。确认无误后加 --apply。")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
