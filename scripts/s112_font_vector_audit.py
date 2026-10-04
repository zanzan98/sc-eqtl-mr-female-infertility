# -*- coding: utf-8 -*-
"""
s112_font_vector_audit.py  ——  只读审计：figures/ 各矢量 PDF 的
  ① 字体族与是否嵌入（FontFile2 = TrueType/Type42，FontFile3 = CFF，FontFile = Type1）
  ② 文字是否为"真文字层"（可提取）还是已被描边化（outline，提取为空）
  ③ 是否含位图（内嵌 raster）
Nature 硬要求：字体须嵌入且为 TrueType 2 或 42；文字不得 outline；主图须矢量可编辑。

★ 只读，不改任何 figures/ 正本。
★ 结果写文件再读（本机 stdout 为 GBK，禁 print 非 ASCII）。
"""
import io
import os
import pymupdf  # import pymupdf as fitz 亦可

ROOT = r"D:/endometriosis_project/11_sc_eqtl_mr_project"
FIGD = os.path.join(ROOT, "figures")
LOG = os.path.join(ROOT, "logs", "_s112_font_vector_audit.txt")

IN_PAPER = [
    "Fig1_study_design",
    "Fig2_discovery_replication",
    "Fig3_coloc_sensitivity",
    "Fig4_chr1p36_12_finemap",
    "Fig5_opentargets_phewas_safety",
    "FigS1_instrument_power_gating",
    "FigS2_sensitivity_outcome_concordance",
    "FigS3_manhattan_by_celltype",
    "FigS4_discovery_volcano",
]
NOT_IN_PAPER = [
    "FigD1_discovery_volcano",
    "FigD2_discovery_manhattan",
    "FigD3_discovery_forest",
    "FigD4_discovery_celltype",
    "FigD5_chr1_LD_structure",
    "FigOverview_all_figures",
]

L = []


def w(s=""):
    L.append(str(s))  # ★ 只落盘，不 print


def fontfile_kind(doc, xref):
    """
    xref = 字体资源 xref。返回嵌入类型字符串。
    在 FontDescriptor 里找 FontFile2 / FontFile3 / FontFile。
    """
    seen = set()
    stack = [xref]
    while stack:
        x = stack.pop()
        if x in seen:
            continue
        seen.add(x)
        try:
            obj = doc.xref_object(x, compressed=True)
        except Exception:
            obj = ""
        if "/FontFile2" in obj:
            return "FontFile2(TrueType=2/42) OK"
        if "/FontFile3" in obj:
            return "FontFile3(CFF/OpenType) ~"
        if "/FontFile" in obj:
            return "FontFile(Type1) BAD"
        if "/FontDescriptor" in obj:
            # 顺着 FontDescriptor 找
            import re
            m = re.search(r"/FontDescriptor\s+(\d+)\s+\d+\s+R", obj)
            if m:
                stack.append(int(m.group(1)))
            # 也可能 FontDescriptor 是内联字典
        if "/Font" in obj and "/FontDescriptor" not in obj:
            import re
            for m in re.finditer(r"/DescendantFonts\s*\[\s*(\d+)\s+\d+\s+R", obj):
                stack.append(int(m.group(1)))
    return "未找到嵌入字体文件(?)"


def inspect(pdf):
    d = pymupdf.open(pdf)
    rec = {"pages": d.page_count, "fonts": [], "txt_chars": 0, "imgs": 0,
           "fontfile": set()}
    for pno in range(d.page_count):
        pg = d[pno]
        rec["txt_chars"] += len(pg.get_text().strip())
        rec["imgs"] += len(pg.get_images(full=True))
        try:
            fl = d.get_page_fonts(pno)
        except Exception:
            fl = []
        for t in fl:
            xref, ext, typ, basefont, name, enc = t[:6]
            rec["fonts"].append((typ, basefont, str(enc)))
            rec["fontfile"].add(fontfile_kind(d, xref))
    d.close()
    return rec


def main():
    w("=" * 100)
    w("figures/ 矢量 PDF · 字体嵌入 / 文字可编辑性 / 位图 审计（Nature 硬要求对照）")
    w("=" * 100)
    for tag, names in (("入论文", IN_PAPER), ("未入论文", NOT_IN_PAPER)):
        w("")
        w("### %s" % tag)
        w("%-40s %-6s %-8s %-7s %-9s %s" %
          ("file", "页", "真文字字符", "位图", "字体类型", "嵌入类型"))
        w("-" * 100)
        for n in names:
            p = os.path.join(FIGD, n + ".pdf")
            if not os.path.exists(p):
                w("%-40s  !! 缺文件" % n)
                continue
            r = inspect(p)
            typs = sorted(set(t[0] for t in r["fonts"]))
            bases = sorted(set(t[1] for t in r["fonts"]))
            ff = sorted(r["fontfile"])
            w("%-40s %-6d %-8d %-7d %-9s %s" %
              (n, r["pages"], r["txt_chars"], r["imgs"],
               "/".join(typs), " | ".join(ff)))
            w("%-40s     basefont: %s" % ("", ", ".join(bases)))

    # ---- 判定 ----
    w("")
    w("=" * 100)
    w("判定")
    w("=" * 100)
    bad_outline = []
    bad_embed = []
    nonstd_font = []
    has_raster = []
    for n in IN_PAPER:
        p = os.path.join(FIGD, n + ".pdf")
        if not os.path.exists(p):
            continue
        r = inspect(p)
        if r["txt_chars"] == 0:
            bad_outline.append(n)          # 提不出文字 ⇒ 疑似 outline
        if not any("FontFile2" in s for s in r["fontfile"]):
            bad_embed.append((n, " | ".join(sorted(r["fontfile"]))))
        if any(("Arial" not in b and "Helvetica" not in b)
               for _, b, _ in r["fonts"]):
            nonstd_font.append((n, ", ".join(sorted(set(b for _, b, _ in r["fonts"])))))
        if r["imgs"] > 0:
            has_raster.append((n, r["imgs"]))

    w("· 提不出真文字（疑 outline，Nature 禁）：%s" %
      (bad_outline if bad_outline else "无"))
    w("· 未见 TrueType(FontFile2) 嵌入（Nature 要求 TrueType 2/42）：%s" %
      (bad_embed if bad_embed else "无"))
    w("· 含非 Arial/Helvetica 族字体：%s" % (nonstd_font if nonstd_font else "无"))
    w("· 含内嵌位图（raster）的面板：%s" % (has_raster if has_raster else "无"))
    w("")
    w("注：FontFile2 说明字体以 TrueType 嵌入（matplotlib pdf.fonttype=42 即为 Type42）；")
    w("    能提取出文字 = 文字层可编辑、未被描边化。")

    with io.open(LOG, "w", encoding="utf-8", newline="\n") as fh:
        fh.write("\n".join(L) + "\n")
    print("done -> %s" % LOG)


if __name__ == "__main__":
    main()
