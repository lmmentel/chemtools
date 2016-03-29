AUGTZ_MOLPRO = '''basis={
!
! HYDROGEN       (5s,2p) -> [3s,2p]
! HYDROGEN       (4s,1p) -> [2s,1p]
! HYDROGEN       (1s,1p)
s, H , 13.0100000, 1.9620000, 0.4446000, 0.1220000, 0.0297400
c, 1.3, 0.0196850, 0.1379770, 0.4781480
c, 4.4, 1
c, 5.5, 1
p, H , 0.7270000, 0.1410000
c, 1.1, 1
c, 2.2, 1
! LITHIUM       (10s,5p,2d) -> [4s,3p,2d]
! LITHIUM       (9s,4p,1d) -> [3s,2p,1d]
! LITHIUM       (1s,1p,1d)
s, LI , 1469.0000000, 220.5000000, 50.2600000, 14.2400000, 4.5810000, 1.5800000, 0.5640000, 0.0734500, 0.0280500, 0.0086400
c, 1.8, 0.0007660, 0.0058920, 0.0296710, 0.1091800, 0.2827890, 0.4531230, 0.2747740, 0.0097510
c, 1.8, -0.0001200, -0.0009230, -0.0046890, -0.0176820, -0.0489020, -0.0960090, -0.1363800, 0.5751020
c, 9.9, 1
c, 10.10, 1
p, LI , 1.5340000, 0.2749000, 0.0736200, 0.0240300, 0.0057900
c, 1.3, 0.0227840, 0.1391070, 0.5003750
c, 4.4, 1
c, 5.5, 1
d, LI , 0.1239000, 0.0725000
c, 1.1, 1
c, 2.2, 1
! BERYLLIUM       (10s,5p,2d) -> [4s,3p,2d]
! BERYLLIUM       (9s,4p,1d) -> [3s,2p,1d]
! BERYLLIUM       (1s,1p,1d)
s, BE , 2940.0000000, 441.2000000, 100.5000000, 28.4300000, 9.1690000, 3.1960000, 1.1590000, 0.1811000, 0.0589000, 0.0187700
c, 1.8, 0.0006800, 0.0052360, 0.0266060, 0.0999930, 0.2697020, 0.4514690, 0.2950740, 0.0125870
c, 1.8, -0.0001230, -0.0009660, -0.0048310, -0.0193140, -0.0532800, -0.1207230, -0.1334350, 0.5307670
c, 9.9, 1
c, 10.10, 1
p, BE , 3.6190000, 0.7110000, 0.1951000, 0.0601800, 0.0085000
c, 1.3, 0.0291110, 0.1693650, 0.5134580
c, 4.4, 1
c, 5.5, 1
d, BE , 0.2380000, 0.0740000
c, 1.1, 1
c, 2.2, 1
}'''

AUGTZ_GAMESS = '''$DATA
HYDROGEN
S   3
  1     13.0100000              0.0196850        
  2      1.9620000              0.1379770        
  3      0.4446000              0.4781480        
S   1
  1      0.1220000              1.0000000        
S   1
  1      0.0297400              1.0000000        
P   1
  1      0.7270000              1.0000000        
P   1
  1      0.1410000              1.0000000        

LITHIUM
S   8
  1   1469.0000000              0.0007660        
  2    220.5000000              0.0058920        
  3     50.2600000              0.0296710        
  4     14.2400000              0.1091800        
  5      4.5810000              0.2827890        
  6      1.5800000              0.4531230        
  7      0.5640000              0.2747740        
  8      0.0734500              0.0097510        
S   8
  1   1469.0000000             -0.0001200        
  2    220.5000000             -0.0009230        
  3     50.2600000             -0.0046890        
  4     14.2400000             -0.0176820        
  5      4.5810000             -0.0489020        
  6      1.5800000             -0.0960090        
  7      0.5640000             -0.1363800        
  8      0.0734500              0.5751020        
S   1
  1      0.0280500              1.0000000        
S   1
  1      0.0086400              1.0000000        
P   3
  1      1.5340000              0.0227840        
  2      0.2749000              0.1391070        
  3      0.0736200              0.5003750        
P   1
  1      0.0240300              1.0000000        
P   1
  1      0.0057900              1.0000000        
D   1
  1      0.1239000              1.0000000        
D   1
  1      0.0725000              1.0000000        

BERYLLIUM
S   8
  1   2940.0000000              0.0006800        
  2    441.2000000              0.0052360        
  3    100.5000000              0.0266060        
  4     28.4300000              0.0999930        
  5      9.1690000              0.2697020        
  6      3.1960000              0.4514690        
  7      1.1590000              0.2950740        
  8      0.1811000              0.0125870        
S   8
  1   2940.0000000             -0.0001230        
  2    441.2000000             -0.0009660        
  3    100.5000000             -0.0048310        
  4     28.4300000             -0.0193140        
  5      9.1690000             -0.0532800        
  6      3.1960000             -0.1207230        
  7      1.1590000             -0.1334350        
  8      0.1811000              0.5307670        
S   1
  1      0.0589000              1.0000000        
S   1
  1      0.0187700              1.0000000        
P   3
  1      3.6190000              0.0291110        
  2      0.7110000              0.1693650        
  3      0.1951000              0.5134580        
P   1
  1      0.0601800              1.0000000        
P   1
  1      0.0085000              1.0000000        
D   1
  1      0.2380000              1.0000000        
D   1
  1      0.0740000              1.0000000        
$END'''

AUGTZ_GAUSSIAN = '''****
H     0 
S   3   1.00
     13.0100000              0.0196850        
      1.9620000              0.1379770        
      0.4446000              0.4781480        
S   1   1.00
      0.1220000              1.0000000        
S   1   1.00
      0.0297400              1.0000000        
P   1   1.00
      0.7270000              1.0000000        
P   1   1.00
      0.1410000              1.0000000        
****
Li     0 
S   8   1.00
   1469.0000000              0.0007660        
    220.5000000              0.0058920        
     50.2600000              0.0296710        
     14.2400000              0.1091800        
      4.5810000              0.2827890        
      1.5800000              0.4531230        
      0.5640000              0.2747740        
      0.0734500              0.0097510        
S   8   1.00
   1469.0000000             -0.0001200        
    220.5000000             -0.0009230        
     50.2600000             -0.0046890        
     14.2400000             -0.0176820        
      4.5810000             -0.0489020        
      1.5800000             -0.0960090        
      0.5640000             -0.1363800        
      0.0734500              0.5751020        
S   1   1.00
      0.0280500              1.0000000        
S   1   1.00
      0.0086400              1.0000000        
P   3   1.00
      1.5340000              0.0227840        
      0.2749000              0.1391070        
      0.0736200              0.5003750        
P   1   1.00
      0.0240300              1.0000000        
P   1   1.00
      0.0057900              1.0000000        
D   1   1.00
      0.1239000              1.0000000        
D   1   1.00
      0.0725000              1.0000000        
****
Be     0 
S   8   1.00
   2940.0000000              0.0006800        
    441.2000000              0.0052360        
    100.5000000              0.0266060        
     28.4300000              0.0999930        
      9.1690000              0.2697020        
      3.1960000              0.4514690        
      1.1590000              0.2950740        
      0.1811000              0.0125870        
S   8   1.00
   2940.0000000             -0.0001230        
    441.2000000             -0.0009660        
    100.5000000             -0.0048310        
     28.4300000             -0.0193140        
      9.1690000             -0.0532800        
      3.1960000             -0.1207230        
      1.1590000             -0.1334350        
      0.1811000              0.5307670        
S   1   1.00
      0.0589000              1.0000000        
S   1   1.00
      0.0187700              1.0000000        
P   3   1.00
      3.6190000              0.0291110        
      0.7110000              0.1693650        
      0.1951000              0.5134580        
P   1   1.00
      0.0601800              1.0000000        
P   1   1.00
      0.0085000              1.0000000        
D   1   1.00
      0.2380000              1.0000000        
D   1   1.00
      0.0740000              1.0000000        
****'''

from chemtools.basisset import BasisSet

def test_basis_set_parser_molpro(tmpdir):

    tmpdir.chdir()
    fpath = tmpdir.join('augtz.molpro')
    fpath.write(AUGTZ_MOLPRO)

    bsd = BasisSet.from_file(str(fpath), fmt='molpro')
    assert isinstance(bsd, dict)
    assert 'H' in bsd.keys()
    assert 'Li' in  bsd.keys()
    assert 'Be' in bsd.keys()
    assert isinstance(bsd['H'], BasisSet)
    assert isinstance(bsd['Li'], BasisSet)
    assert isinstance(bsd['Be'], BasisSet)
    assert bsd['H'].nf() == 9
    assert bsd['Li'].nf() == 23
    assert bsd['Be'].nf() == 23

def test_basis_set_parser_gamess(tmpdir):

    tmpdir.chdir()
    fpath = tmpdir.join('augtz.gamess')
    fpath.write(AUGTZ_GAMESS)

    bsd = BasisSet.from_file(str(fpath), fmt='gamessus')
    assert isinstance(bsd, dict)
    assert 'H' in bsd.keys()
    assert 'Li' in  bsd.keys()
    assert 'Be' in bsd.keys()
    assert isinstance(bsd['H'], BasisSet)
    assert isinstance(bsd['Li'], BasisSet)
    assert isinstance(bsd['Be'], BasisSet)
    assert bsd['H'].nf() == 9
    assert bsd['Li'].nf() == 23
    assert bsd['Be'].nf() == 23

def test_basis_set_parser_gaussian(tmpdir):

    tmpdir.chdir()
    fpath = tmpdir.join('augtz.gaussian')
    fpath.write(AUGTZ_GAUSSIAN)

    bsd = BasisSet.from_file(str(fpath), fmt='gaussian')
    assert isinstance(bsd, dict)
    assert 'H' in bsd.keys()
    assert 'Li' in  bsd.keys()
    assert 'Be' in bsd.keys()
    assert isinstance(bsd['H'], BasisSet)
    assert isinstance(bsd['Li'], BasisSet)
    assert isinstance(bsd['Be'], BasisSet)
    assert bsd['H'].nf() == 9
    assert bsd['Li'].nf() == 23
    assert bsd['Be'].nf() == 23

