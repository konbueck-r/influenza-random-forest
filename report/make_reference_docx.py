"""Build report/reference.docx, the Word style template for the Quarto report.

Starts from pandoc's default reference.docx and changes fonts, sizes, colours, paragraph
spacing, page size and margins, and adds a centred page number in the footer.

Usage (from the project root):
    quarto pandoc -o report/reference_default.docx --print-default-data-file reference.docx
    python report/make_reference_docx.py
"""

import re
import zipfile
from pathlib import Path

SOURCE = Path("report/reference_default.docx")
TARGET = Path("report/reference.docx")

BODY_FONT = "Calibri"
HEADING_FONT = "Calibri"
ACCENT = "1F5A96"  # the "wt infected" blue of the figures


def fonts(name):
    return (f'<w:rFonts w:ascii="{name}" w:hAnsi="{name}" w:eastAsia="{name}" '
            f'w:cs="{name}" />')


def replace_style(styles, style_id, ppr, rpr):
    """Replace the pPr and rPr of one paragraph style, keeping its name and links."""
    pattern = re.compile(r'(<w:style [^>]*w:styleId="%s">)(.*?)(</w:style>)' % style_id, re.S)
    match = pattern.search(styles)
    if match is None:
        raise ValueError(f"style {style_id} not found")
    body = re.sub(r'<w:pPr>.*?</w:pPr>|<w:pPr ?/>|<w:rPr>.*?</w:rPr>|<w:rPr ?/>', '',
                  match.group(2), flags=re.S)
    new = match.group(1) + body + f'<w:pPr>{ppr}</w:pPr><w:rPr>{rpr}</w:rPr>' + match.group(3)
    return styles[:match.start()] + new + styles[match.end():]


def size(half_points):
    return f'<w:sz w:val="{half_points}" /><w:szCs w:val="{half_points}" />'


with zipfile.ZipFile(SOURCE) as archive:
    parts = {name: archive.read(name) for name in archive.namelist()}

styles = parts["word/styles.xml"].decode("utf-8")

# Document defaults: 11 pt body font, 1.15 line spacing, 6 pt after paragraphs.
styles = re.sub(
    r'<w:docDefaults>.*?</w:docDefaults>',
    '<w:docDefaults><w:rPrDefault><w:rPr>' + fonts(BODY_FONT) + size(22) +
    '<w:lang w:val="en-GB" /></w:rPr></w:rPrDefault><w:pPrDefault><w:pPr>'
    '<w:spacing w:after="120" w:line="276" w:lineRule="auto" /></w:pPr></w:pPrDefault>'
    '</w:docDefaults>',
    styles, flags=re.S)

styles = replace_style(styles, "BodyText",
                       '<w:spacing w:before="0" w:after="120" /><w:jc w:val="both" />', '')
styles = replace_style(styles, "Title",
                       '<w:spacing w:before="0" w:after="160" w:line="240" w:lineRule="auto" />'
                       '<w:jc w:val="left" />',
                       fonts(HEADING_FONT) + f'<w:b /><w:color w:val="{ACCENT}" />' + size(36))
styles = replace_style(styles, "Subtitle",
                       '<w:spacing w:before="0" w:after="240" /><w:jc w:val="left" />',
                       fonts(HEADING_FONT) + '<w:color w:val="404040" />' + size(24))
styles = replace_style(styles, "Author", '<w:spacing w:before="0" w:after="60" />'
                       '<w:jc w:val="left" />', fonts(BODY_FONT) + size(22))
styles = replace_style(styles, "Date", '<w:spacing w:before="0" w:after="240" />'
                       '<w:jc w:val="left" />', fonts(BODY_FONT) + size(22))
styles = replace_style(styles, "AbstractTitle", '<w:keepNext /><w:spacing w:before="240" '
                       'w:after="60" />', fonts(HEADING_FONT) + '<w:b />' +
                       f'<w:color w:val="{ACCENT}" />' + size(22))
styles = replace_style(styles, "Abstract",
                       '<w:spacing w:before="0" w:after="240" /><w:jc w:val="both" />'
                       '<w:ind w:left="284" w:right="284" />', size(21))
styles = replace_style(styles, "Heading1",
                       '<w:keepNext /><w:keepLines /><w:spacing w:before="360" w:after="120" />'
                       '<w:outlineLvl w:val="0" />',
                       fonts(HEADING_FONT) + f'<w:b /><w:color w:val="{ACCENT}" />' + size(30))
styles = replace_style(styles, "Heading2",
                       '<w:keepNext /><w:keepLines /><w:spacing w:before="240" w:after="80" />'
                       '<w:outlineLvl w:val="1" />',
                       fonts(HEADING_FONT) + f'<w:b /><w:color w:val="{ACCENT}" />' + size(24))
styles = replace_style(styles, "Heading3",
                       '<w:keepNext /><w:keepLines /><w:spacing w:before="200" w:after="60" />'
                       '<w:outlineLvl w:val="2" />',
                       fonts(HEADING_FONT) + '<w:b /><w:i />' + size(22))
styles = replace_style(styles, "Caption",
                       '<w:spacing w:before="60" w:after="240" w:line="240" w:lineRule="auto" />'
                       '<w:jc w:val="both" />', size(18))
styles = replace_style(styles, "TableCaption",
                       '<w:keepNext /><w:keepLines /><w:spacing w:before="240" w:after="80" '
                       'w:line="240" w:lineRule="auto" /><w:jc w:val="both" />', size(18))
styles = replace_style(styles, "Bibliography",
                       '<w:spacing w:before="0" w:after="60" w:line="240" w:lineRule="auto" />'
                       '<w:ind w:left="397" w:hanging="397" />', size(18))
# Quarto puts each figure image in a Compact paragraph directly above its caption; keepNext
# holds image and caption on one page (it also keeps table rows together, which is fine here).
styles = replace_style(styles, "Compact", '<w:keepNext /><w:spacing w:before="20" w:after="20" '
                       'w:line="240" w:lineRule="auto" />', size(18))
# Quarto also uses Image Caption for table captions (above the table), hence keepNext.
styles = replace_style(styles, "ImageCaption",
                       '<w:keepNext /><w:keepLines /><w:spacing w:before="60" w:after="240" w:line="240" '
                       'w:lineRule="auto" /><w:jc w:val="both" />', size(18))

parts["word/styles.xml"] = styles.encode("utf-8")

# A4, 2.5 cm margins (16 cm text width, the width the figures were made for), page number footer.
footer = ('<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
          '<w:ftr xmlns:w="http://schemas.openxmlformats.org/wordprocessingml/2006/main" '
          'xmlns:r="http://schemas.openxmlformats.org/officeDocument/2006/relationships">'
          '<w:p><w:pPr><w:jc w:val="center" /></w:pPr>'
          '<w:r><w:rPr>' + size(18) + '</w:rPr><w:fldChar w:fldCharType="begin" /></w:r>'
          '<w:r><w:rPr>' + size(18) + '</w:rPr><w:instrText xml:space="preserve"> PAGE '
          '</w:instrText></w:r>'
          '<w:r><w:rPr>' + size(18) + '</w:rPr><w:fldChar w:fldCharType="end" /></w:r>'
          '</w:p></w:ftr>')
parts["word/footer1.xml"] = footer.encode("utf-8")

rels = parts["word/_rels/document.xml.rels"].decode("utf-8")
rels = rels.replace("</Relationships>",
                    '<Relationship Type="http://schemas.openxmlformats.org/officeDocument/2006/'
                    'relationships/footer" Id="rIdFooter1" Target="footer1.xml" />'
                    "</Relationships>")
parts["word/_rels/document.xml.rels"] = rels.encode("utf-8")

types = parts["[Content_Types].xml"].decode("utf-8")
types = types.replace("</Types>",
                      '<Override PartName="/word/footer1.xml" ContentType="application/'
                      'vnd.openxmlformats-officedocument.wordprocessingml.footer+xml" /></Types>')
parts["[Content_Types].xml"] = types.encode("utf-8")

document = parts["word/document.xml"].decode("utf-8")
document = re.sub(
    r'<w:sectPr>.*?</w:sectPr>',
    '<w:sectPr><w:footerReference w:type="default" r:id="rIdFooter1" />'
    '<w:footnotePr><w:numRestart w:val="eachSect" /></w:footnotePr>'
    '<w:pgSz w:w="11906" w:h="16838" />'
    '<w:pgMar w:top="1417" w:right="1417" w:bottom="1417" w:left="1417" w:header="709" '
    'w:footer="709" w:gutter="0" /></w:sectPr>',
    document, flags=re.S)
if 'xmlns:r=' not in document[:2000]:
    raise ValueError("document.xml lacks the relationships namespace")
parts["word/document.xml"] = document.encode("utf-8")

with zipfile.ZipFile(TARGET, "w", zipfile.ZIP_DEFLATED) as archive:
    for name, data in parts.items():
        archive.writestr(name, data)
print(f"Written {TARGET}")
