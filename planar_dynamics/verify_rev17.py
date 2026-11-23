"""Verify archived output bytes and optionally the exact manuscript table/figure contents."""
from hashlib import sha256
from pathlib import Path
import sys
import zipfile
import xml.etree.ElementTree as ET
from common17 import ROOT
from published_tables import table3,table4,table5

FIGURES={1:'fig1_trajectories.png',4:'fig4_crossing_counts.png',5:'fig5_first_crossing.png',
         6:'fig6_uncertainty.png',7:'fig7_manifold_splitting.png',8:'fig8_poincare.png'}

def verify_doc(path):
    ns={'w':'http://schemas.openxmlformats.org/wordprocessingml/2006/main',
        'a':'http://schemas.openxmlformats.org/drawingml/2006/main',
        'r':'http://schemas.openxmlformats.org/officeDocument/2006/relationships'}
    with zipfile.ZipFile(path) as z:
        doc=ET.fromstring(z.read('word/document.xml'))
        tables=doc.findall('.//w:tbl',ns)
        for i,expected in [(2,table3()),(3,table4()),(4,table5())]:
            actual=[[''.join(el.text or '' for el in cell.iter() if el.tag.split('}')[-1]=='t')
                     for cell in row.findall('w:tc',ns)] for row in tables[i].findall('w:tr',ns)[1:]]
            if actual!=expected: raise ValueError(f'Table {i+1} differs')
        rels={e.get('Id'):e.get('Target') for e in ET.fromstring(z.read('word/_rels/document.xml.rels'))}
        for i,name in FIGURES.items():
            blip=doc.findall('.//a:blip',ns)[i-1]
            member='word/'+rels[blip.get('{'+ns['r']+'}embed')]
            if z.read(member)!=(ROOT/'figures'/name).read_bytes(): raise ValueError(f'Figure {i} differs')
    print('Word tables 3–5 and regenerated figure bytes match.')

def run():
    for line in (ROOT/'SHA256SUMS').read_text().splitlines():
        digest,rel=line.split('  ',1)
        if sha256((ROOT/rel).read_bytes()).hexdigest()!=digest:
            raise ValueError(f'Checksum mismatch: {rel}. Compare the recorded software environment.')
    print('All reference file checksums match.')
    for i,rows in ((3,table3()),(4,table4()),(5,table5())):
        print('Table',i)
        for row in rows: print(' | '.join(row))
    if len(sys.argv)>1: verify_doc(sys.argv[1])

if __name__=='__main__':
    run()
