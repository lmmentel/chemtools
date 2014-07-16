# -*- coding: utf-8 -*-

import unittest
import tempfile

from chemtools.gamessus import GamessLogParser

class TestGLPonHeH2(unittest.TestCase):
    heh2_log = '''----- GAMESS execution script 'rungms' -----
This job is running on host doreen
under operating system Linux at nie, 11 maj 2014, 14:57:29 CEST
Available scratch disk space (Kbyte units) at beginning of the job is
System plików  1K-blocks     użyte dostępne %uż. zamont. na
/dev/sda5      264208320 221156900 29607380  89% /home
GAMESS temporary binary files will be written to /home/lmentel/scratch
GAMESS supplementary output files will be written to /home/lmentel/scratch
Copying input file he-h2_avdz_ormas.inp to your run's scratch directory...
cp he-h2_avdz_ormas.inp /home/lmentel/scratch/he-h2_avdz_ormas.F05
unset echo
/home/lmentel/Programs/gamess-us-may2013/ddikick.x /home/lmentel/Programs/gamess-us-may2013/gamess.01.x he-h2_avdz_ormas -ddi 1 1 doreen -scr /home/lmentel/scratch

 Distributed Data Interface kickoff program.
 Initiating 1 compute processes on 1 nodes to run the following command:
 /home/lmentel/Programs/gamess-us-may2013/gamess.01.x he-h2_avdz_ormas 

          ******************************************************
          *         GAMESS VERSION =  1 MAY 2013 (R1)          *
          *             FROM IOWA STATE UNIVERSITY             *
          * M.W.SCHMIDT, K.K.BALDRIDGE, J.A.BOATZ, S.T.ELBERT, *
          *   M.S.GORDON, J.H.JENSEN, S.KOSEKI, N.MATSUNAGA,   *
          *          K.A.NGUYEN, S.J.SU, T.L.WINDUS,           *
          *       TOGETHER WITH M.DUPUIS, J.A.MONTGOMERY       *
          *         J.COMPUT.CHEM.  14, 1347-1363(1993)        *
          **************** 64 BIT LINUX VERSION ****************

  SINCE 1993, STUDENTS AND POSTDOCS WORKING AT IOWA STATE UNIVERSITY
  AND ALSO IN THEIR VARIOUS JOBS AFTER LEAVING ISU HAVE MADE IMPORTANT
  CONTRIBUTIONS TO THE CODE:
     IVANA ADAMOVIC, CHRISTINE AIKENS, YURI ALEXEEV, POOJA ARORA,
     ANDREY ASADCHEV, ROB BELL, PRADIPTA BANDYOPADHYAY, JONATHAN BENTZ,
     BRETT BODE, GALINA CHABAN, WEI CHEN, CHEOL HO CHOI, PAUL DAY,
     ALBERT DEFUSCO, TIM DUDLEY, DMITRI FEDOROV, GRAHAM FLETCHER,
     MARK FREITAG, KURT GLAESEMANN, DAN KEMP, GRANT MERRILL,
     NORIYUKI MINEZAWA, JONATHAN MULLIN, TAKESHI NAGATA,
     SEAN NEDD, HEATHER NETZLOFF, BOSILJKA NJEGIC, RYAN OLSON, MIKE PAK,
     JIM SHOEMAKER, LYUDMILA SLIPCHENKO, SAROM SOK, JIE SONG,
     TETSUYA TAKETSUGU, SIMON WEBB, SOOHAENG YOO, FEDERICO ZAHARIEV

  ADDITIONAL CODE HAS BEEN PROVIDED BY COLLABORATORS IN OTHER GROUPS:
     IOWA STATE UNIVERSITY:
          JOE IVANIC, LAIMUTIS BYTAUTAS, KLAUS RUEDENBERG
     UNIVERSITY OF TOKYO: KIMIHIKO HIRAO, TAKAHITO NAKAJIMA,
          TAKAO TSUNEDA, MUNEAKI KAMIYA, SUSUMU YANAGISAWA,
          KIYOSHI YAGI, MAHITO CHIBA, SEIKEN TOKURA, NAOAKI KAWAKAMI
     UNIVERSITY OF AARHUS: FRANK JENSEN
     UNIVERSITY OF IOWA: VISVALDAS KAIRYS, HUI LI
     NATIONAL INST. OF STANDARDS AND TECHNOLOGY: WALT STEVENS, DAVID GARMER
     UNIVERSITY OF PISA: BENEDETTA MENNUCCI, JACOPO TOMASI
     UNIVERSITY OF MEMPHIS: HENRY KURTZ, PRAKASHAN KORAMBATH
     UNIVERSITY OF ALBERTA: TOBY ZENG, MARIUSZ KLOBUKOWSKI
     UNIVERSITY OF NEW ENGLAND: MARK SPACKMAN
     MIE UNIVERSITY: HIROAKI UMEDA
     MICHIGAN STATE UNIVERSITY:
          KAROL KOWALSKI, MARTA WLOCH, JEFFREY GOUR, JESSE LUTZ,
          WEI LI, PIOTR PIECUCH
     UNIVERSITY OF SILESIA: MONIKA MUSIAL, STANISLAW KUCHARSKI
     FACULTES UNIVERSITAIRES NOTRE-DAME DE LA PAIX:
          OLIVIER QUINET, BENOIT CHAMPAGNE
     UNIVERSITY OF CALIFORNIA - SANTA BARBARA: BERNARD KIRTMAN
     INSTITUTE FOR MOLECULAR SCIENCE:
          KAZUYA ISHIMURA, MICHIO KATOUDA, AND SHIGERU NAGASE
     UNIVERSITY OF NOTRE DAME: DAN CHIPMAN
     KYUSHU UNIVERSITY:
          HARUYUKI NAKANO,
          FENG LONG GU, JACEK KORCHOWIEC, MARCIN MAKOWSKI, AND YURIKO AOKI,
          HIROTOSHI MORI AND EISAKU MIYOSHI
     PENNSYLVANIA STATE UNIVERSITY:
          TZVETELIN IORDANOV, CHET SWALINA, JONATHAN SKONE,
          SHARON HAMMES-SCHIFFER
     WASEDA UNIVERSITY:
          MASATO KOBAYASHI, TOMOKO AKAMA, TSUGUKI TOUMA,
          TAKESHI YOSHIKAWA, YASUHIRO IKABATA, HIROMI NAKAI
     NANJING UNIVERSITY: SHUHUA LI
     UNIVERSITY OF NEBRASKA:
          PEIFENG SU, DEJUN SI, NANDUN THELLAMUREGE, YALI WANG, HUI LI
     UNIVERSITY OF ZURICH:
          ROBERTO PEVERATI, KIM BALDRIDGE
     N. COPERNICUS UNIVERSITY AND JACKSON STATE UNIVERSITY:
          MARIA BARYSZ

 EXECUTION OF GAMESS BEGUN Sun May 11 14:57:30 2014

            ECHO OF THE FIRST FEW INPUT CARDS -
 INPUT CARD> $CONTRL                                                                        
 INPUT CARD>    scftyp=rhf                                                                  
 INPUT CARD>    cityp=ormas                                                                 
 INPUT CARD>    runtyp=energy                                                               
 INPUT CARD>    maxit=30                                                                    
 INPUT CARD>    mult=1                                                                      
 INPUT CARD>    icut=30                                                                     
 INPUT CARD>    itol=30                                                                     
 INPUT CARD>    ispher=1                                                                    
 INPUT CARD>    units=bohr                                                                  
 INPUT CARD>!    exetyp=check                                                               
 INPUT CARD> $END                                                                           
 INPUT CARD> $SYSTEM timlim=525600 mwords=100 $END                                          
 INPUT CARD> $SCF    conv=1.0d-8    $END                                                    
 INPUT CARD> $TRANS  cuttrf=1.0d-10 $END                                                    
 INPUT CARD> $CIDET                                                                         
 INPUT CARD>    ncore=0                                                                     
 INPUT CARD>    nact=27                                                                     
 INPUT CARD>    nels=4                                                                      
 INPUT CARD>    sz=0                                                                        
 INPUT CARD>    analys=.true.                                                               
 INPUT CARD>    group=c2v                                                                   
 INPUT CARD>    stsym=a1                                                                    
 INPUT CARD>!    nstate=1                                                                   
 INPUT CARD>!    itermx=100                                                                 
 INPUT CARD>!    cvgtol=1.0d-6                                                              
 INPUT CARD>!    nflgdm=1                                                                   
 INPUT CARD> $END                                                                           
 INPUT CARD> $ORMAS                                                                         
 INPUT CARD>    nspace=2                                                                    
 INPUT CARD>    mstart(1)=1,3                                                               
 INPUT CARD>    mine(1)=0,0                                                                 
 INPUT CARD>    maxe(1)=4,4                                                                 
 INPUT CARD>    qcorr=.false.                                                               
 INPUT CARD> $END                                                                           
 INPUT CARD>! R_e(H-H) = 1.448 736 a_0                                                      
 INPUT CARD> $DATA                                                                          
 INPUT CARD>He-H2 FCI                                                                       
 INPUT CARD>cnv 2                                                                           
 INPUT CARD>                                                                                
 INPUT CARD>He     2.0     0.00000     0.000000     -3.000000                               
 INPUT CARD>S   3                                                                           
 INPUT CARD>  1     38.3600000              0.0238090                                       
 INPUT CARD>  2      5.7700000              0.1548910                                       
 INPUT CARD>  3      1.2400000              0.4699870                                       
 INPUT CARD>S   1                                                                           
 INPUT CARD>  1      0.2976000              1.0000000                                       
 INPUT CARD>S   1                                                                           
 INPUT CARD>  1      0.0725500              1.0000000                                       
 INPUT CARD>P   1                                                                           
  100000000 WORDS OF MEMORY AVAILABLE


     RUN TITLE
     ---------
 He-H2 FCI                                                                       

 THE POINT GROUP OF THE MOLECULE IS CNV     
 THE ORDER OF THE PRINCIPAL AXIS IS     2
 *** WARNING! ATOM   1 SHELL    1 TYPE S HAS NORMALIZATION   1.68743229
 *** WARNING! ATOM   2 SHELL    6 TYPE S HAS NORMALIZATION   1.70173870

 ATOM      ATOMIC                      COORDINATES (BOHR)
           CHARGE         X                   Y                   Z
 HE          2.0     0.0000000000        0.0000000000       -3.0000000000
 H           1.0    -0.0000000000       -0.7243680000        3.0000000000
 H           1.0     0.0000000000        0.7243680000        3.0000000000

          INTERNUCLEAR DISTANCES (ANGS.)
          ------------------------------

                1 HE         2 H          3 H     

   1 HE      0.0000000    3.1981185    3.1981185  
   2 H       3.1981185    0.0000000    0.7666381 *
   3 H       3.1981185    0.7666381 *  0.0000000  

  * ... LESS THAN  3.000


     ATOMIC BASIS SET
     ----------------
 THE CONTRACTED PRIMITIVE FUNCTIONS HAVE BEEN UNNORMALIZED
 THE CONTRACTED BASIS FUNCTIONS ARE NOW NORMALIZED TO UNITY

  SHELL TYPE  PRIMITIVE        EXPONENT          CONTRACTION COEFFICIENT(S)

 HE        

      1   S       1            38.3600000    0.040176075383
      1   S       2             5.7700000    0.261368074769
      1   S       3             1.2400000    0.793071239495

      2   S       4             0.2976000    1.000000000000

      3   S       5             0.0725500    1.000000000000

      4   P       6             1.2750000    1.000000000000

      5   P       7             0.2473000    1.000000000000

 H         

     11   S       8            13.0100000    0.033498726390
     11   S       9             1.9620000    0.234800801174
     11   S      10             0.4446000    0.813682957883

     12   S      11             0.1220000    1.000000000000

     13   S      12             0.0297400    1.000000000000

     14   P      13             0.7270000    1.000000000000

     15   P      14             0.1410000    1.000000000000

 TOTAL NUMBER OF BASIS SET SHELLS             =   15
 NUMBER OF CARTESIAN GAUSSIAN BASIS FUNCTIONS =   27
 NOTE: THIS RUN WILL RESTRICT THE MO VARIATION SPACE TO SPHERICAL HARMONICS.
 THE NUMBER OF ORBITALS KEPT IN THE VARIATIONAL SPACE WILL BE PRINTED LATER.
 NUMBER OF ELECTRONS                          =    4
 CHARGE OF MOLECULE                           =    0
 SPIN MULTIPLICITY                            =    1
 NUMBER OF OCCUPIED ORBITALS (ALPHA)          =    2
 NUMBER OF OCCUPIED ORBITALS (BETA )          =    2
 TOTAL NUMBER OF ATOMS                        =    3
 THE NUCLEAR REPULSION ENERGY IS        1.3521176085

     $CONTRL OPTIONS
     ---------------
 SCFTYP=RHF          RUNTYP=ENERGY       EXETYP=RUN     
 MPLEVL=       0     CITYP =ORMAS        CCTYP =NONE         VBTYP =NONE    
 DFTTYP=NONE         TDDFT =NONE    
 MULT  =       1     ICHARG=       0     NZVAR =       0     COORD =UNIQUE  
 PP    =NONE         RELWFN=NONE         LOCAL =NONE         NUMGRD=       F
 ISPHER=       1     NOSYM =       0     MAXIT =      30     UNITS =BOHR    
 PLTORB=       F     MOLPLT=       F     AIMPAC=       F     FRIEND=        
 NPRINT=       7     IREST =       0     GEOM  =INPUT   
 NORMF =       0     NORMP =       0     ITOL  =      30     ICUT  =      30
 INTTYP=BEST         GRDTYP=BEST         QMTTOL= 1.0E-06

     $SYSTEM OPTIONS
     ---------------
  REPLICATED MEMORY=   100000000 WORDS (ON EVERY NODE).
 DISTRIBUTED MEMDDI=           0 MILLION WORDS IN AGGREGATE,
 MEMDDI DISTRIBUTED OVER   1 PROCESSORS IS           0 WORDS/PROCESSOR.
 TOTAL MEMORY REQUESTED ON EACH PROCESSOR=   100000000 WORDS.
 TIMLIM=      525600.00 MINUTES, OR     365.0 DAYS.
 PARALL= F  BALTYP=  DLB     KDIAG=    0  COREFL= F
 MXSEQ2=     300 MXSEQ3=     150

          ----------------
          PROPERTIES INPUT
          ----------------

     MOMENTS            FIELD           POTENTIAL          DENSITY
 IEMOM =       1   IEFLD =       0   IEPOT =       0   IEDEN =       0
 WHERE =COMASS     WHERE =NUCLEI     WHERE =NUCLEI     WHERE =NUCLEI  
 OUTPUT=BOTH       OUTPUT=BOTH       OUTPUT=BOTH       OUTPUT=BOTH    
 IEMINT=       0   IEFINT=       0                     IEDINT=       0
                                                       MORB  =       0
          EXTRAPOLATION IN EFFECT
          SOSCF IN EFFECT
 ORBITAL PRINTING OPTION: NPREO=     1    27     2     1

     -------------------------------
     INTEGRAL TRANSFORMATION OPTIONS
     -------------------------------
     NWORD  =            0
     CUTOFF = 1.0E-10     MPTRAN =       0
     DIRTRF =       F     AOINTS =DUP     

          ----------------------
          INTEGRAL INPUT OPTIONS
          ----------------------
 NOPK  =       1 NORDER=       0 SCHWRZ=       F

     ------------------------------------------
     THE POINT GROUP IS CNV, NAXIS= 2, ORDER= 4
     ------------------------------------------

 -- VARIATIONAL SPACE WILL BE RESTRICTED TO PURE SPHERICAL HARMONICS ONLY --
 AFTER EXCLUDING CONTAMINANT COMBINATIONS FROM THE CARTESIAN GAUSSIAN BASIS
 SET, THE NUMBER OF SPHERICAL HARMONICS KEPT IN THE VARIATION SPACE IS   27

     DIMENSIONS OF THE SYMMETRY SUBSPACES ARE
 A1  =   12     A2  =    2     B1  =    4     B2  =    9

 ..... DONE SETTING UP THE RUN .....
 STEP CPU TIME =     0.02 TOTAL CPU TIME =        0.0 (    0.0 MIN)
 TOTAL WALL CLOCK TIME=        0.6 SECONDS, CPU UTILIZATION IS   3.64%

          ********************
          1 ELECTRON INTEGRALS
          ********************
 ...... END OF ONE-ELECTRON INTEGRALS ......
 STEP CPU TIME =     0.00 TOTAL CPU TIME =        0.0 (    0.0 MIN)
 TOTAL WALL CLOCK TIME=        0.6 SECONDS, CPU UTILIZATION IS   3.64%

          -------------
          GUESS OPTIONS
          -------------
          GUESS =HUCKEL            NORB  =       0          NORDER=       0
          MIX   =       F          PRTMO =       F          PUNMO =       F
          TOLZ  = 1.0E-08          TOLE  = 1.0E-05
          SYMDEN=       F          PURIFY=       F

 INITIAL GUESS ORBITALS GENERATED BY HUCKEL   ROUTINE.
 HUCKEL GUESS REQUIRES      7473 WORDS.

 SYMMETRIES FOR INITIAL GUESS ORBITALS FOLLOW.   BOTH SET(S).
     2 ORBITALS ARE OCCUPIED (    0 CORE ORBITALS).
     1=A1       2=A1       3=B2       4=A1       5=A1       6=A1       7=A1  
     8=A1       9=A1      10=A1      11=A1      12=A1  
 ...... END OF INITIAL ORBITAL SELECTION ......
 STEP CPU TIME =     0.00 TOTAL CPU TIME =        0.0 (    0.0 MIN)
 TOTAL WALL CLOCK TIME=        0.6 SECONDS, CPU UTILIZATION IS   3.28%

                    ----------------------
                    AO INTEGRAL TECHNOLOGY
                    ----------------------
     S,P,L SHELL ROTATED AXIS INTEGRALS, REPROGRAMMED BY
        KAZUYA ISHIMURA (IMS) AND JOSE SIERRA (SYNSTAR).
     S,P,D,L SHELL ROTATED AXIS INTEGRALS PROGRAMMED BY
        KAZUYA ISHIMURA (INSTITUTE FOR MOLECULAR SCIENCE).
     S,P,D,F,G SHELL TO TOTAL QUARTET ANGULAR MOMENTUM SUM 5,
        ERIC PROGRAM BY GRAHAM FLETCHER (ELORET AND NASA ADVANCED
        SUPERCOMPUTING DIVISION, AMES RESEARCH CENTER).
     S,P,D,F,G,L SHELL GENERAL RYS QUADRATURE PROGRAMMED BY
        MICHEL DUPUIS (PACIFIC NORTHWEST NATIONAL LABORATORY).

          --------------------
          2 ELECTRON INTEGRALS
          --------------------

 THE -PK- OPTION IS OFF, THE INTEGRALS ARE NOT IN SUPERMATRIX FORM.
 STORING   15000 INTEGRALS/RECORD ON DISK, USING 12 BYTES/INTEGRAL.
 TWO ELECTRON INTEGRAL EVALUATION REQUIRES   89476 WORDS OF MEMORY.
 II,JST,KST,LST =  1  1  1  1 NREC =         1 INTLOC =    1
 II,JST,KST,LST =  2  1  1  1 NREC =         1 INTLOC =    2
 II,JST,KST,LST =  3  1  1  1 NREC =         1 INTLOC =    7
 II,JST,KST,LST =  4  1  1  1 NREC =         1 INTLOC =   22
 II,JST,KST,LST =  5  1  1  1 NREC =         1 INTLOC =   67
 II,JST,KST,LST =  6  1  1  1 NREC =         1 INTLOC =  214
 II,JST,KST,LST =  7  1  1  1 NREC =         1 INTLOC =  214
 II,JST,KST,LST =  8  1  1  1 NREC =         1 INTLOC =  214
 II,JST,KST,LST =  9  1  1  1 NREC =         1 INTLOC =  214
 II,JST,KST,LST = 10  1  1  1 NREC =         1 INTLOC =  214
 II,JST,KST,LST = 11  1  1  1 NREC =         1 INTLOC =  214
 II,JST,KST,LST = 12  1  1  1 NREC =         1 INTLOC =  632
 II,JST,KST,LST = 13  1  1  1 NREC =         1 INTLOC = 1389
 II,JST,KST,LST = 14  1  1  1 NREC =         1 INTLOC = 2637
 II,JST,KST,LST = 15  1  1  1 NREC =         1 INTLOC =10056
 TOTAL NUMBER OF NONZERO TWO-ELECTRON INTEGRALS =               26488
          2 INTEGRAL RECORDS WERE STORED ON DISK FILE  8.
  ...... END OF TWO-ELECTRON INTEGRALS .....
 STEP CPU TIME =     0.03 TOTAL CPU TIME =        0.0 (    0.0 MIN)
 TOTAL WALL CLOCK TIME=        0.7 SECONDS, CPU UTILIZATION IS   7.25%

          --------------------------
                 RHF SCF CALCULATION
          --------------------------

     NUCLEAR ENERGY =         1.3521176085
     MAXIT =   30     NPUNCH=    2
     EXTRAP=T  DAMP=F  SHIFT=F  RSTRCT=F  DIIS=F  DEM=F  SOSCF=T
     DENSITY MATRIX CONV=  1.00E-08
     SOSCF WILL OPTIMIZE      50 ORBITAL ROTATIONS, SOGTOL=   0.250
     MEMORY REQUIRED FOR RHF ITERS=     35681 WORDS.

 ITER EX DEM     TOTAL ENERGY        E CHANGE  DENSITY CHANGE     ORB. GRAD
   1  0  0       -3.9698188730    -3.9698188730   0.080318525   0.000000000
          ---------------START SECOND ORDER SCF---------------
   2  1  0       -3.9835229927    -0.0137041197   0.014747883   0.012802713
   3  2  0       -3.9841907258    -0.0006677331   0.004123435   0.002353190
   4  3  0       -3.9842265743    -0.0000358485   0.000087014   0.000102178
   5  4  0       -3.9842266054    -0.0000000311   0.000020173   0.000016561
   6  5  0       -3.9842266066    -0.0000000012   0.000000689   0.000000507
   7  6  0       -3.9842266066    -0.0000000000   0.000000202   0.000000103
   8  7  0       -3.9842266066    -0.0000000000   0.000000000   0.000000008

          -----------------
          DENSITY CONVERGED
          -----------------
     TIME TO FORM FOCK OPERATORS=       0.0 SECONDS (       0.0 SEC/ITER)
     TIME TO SOLVE SCF EQUATIONS=       0.0 SECONDS (       0.0 SEC/ITER)

 FINAL RHF ENERGY IS       -3.9842266066 AFTER   8 ITERATIONS

          ------------
          EIGENVECTORS
          ------------

                      1          2          3          4          5
                   -0.9161    -0.5847     0.0609     0.0712     0.1864
                     A1         A1         B2         A1         A1  
    1  HE 1  S    0.594341  -0.009006  -0.000000  -0.030523  -0.141698
    2  HE 1  S    0.496648  -0.009922  -0.000000  -0.040684  -0.629283
    3  HE 1  S    0.023809  -0.001379  -0.000000  -0.174220   1.345883
    4  HE 1  X    0.000000  -0.000000  -0.000000  -0.000000   0.000000
    5  HE 1  Y    0.000000  -0.000000  -0.001945  -0.000000   0.000000
    6  HE 1  Z    0.000037   0.000315  -0.000000   0.004345  -0.000153
    7  HE 1  X    0.000000  -0.000000  -0.000000  -0.000000   0.000000
    8  HE 1  Y    0.000000  -0.000000  -0.003255  -0.000000   0.000000
    9  HE 1  Z   -0.000621   0.000156  -0.000000  -0.023882   0.006909
   10  H  2  S    0.001662   0.401361   0.050078  -0.060858  -0.009821
   11  H  2  S    0.000899   0.173534  -2.168485  -0.375435   0.065272
   12  H  2  S    0.000559   0.007174   4.326998   0.736068  -0.203641
   13  H  2  X    0.000000  -0.000000  -0.000000  -0.000000   0.000000
   14  H  2  Y    0.000164   0.021696  -0.001837  -0.013136   0.001626
   15  H  2  Z    0.000137  -0.000001   0.000054   0.000515  -0.003718
   16  H  2  X    0.000000  -0.000000  -0.000000  -0.000000   0.000000
   17  H  2  Y    0.000768   0.005097  -0.543254   0.071211  -0.055957
   18  H  2  Z   -0.001100   0.000006  -0.001076   0.006755  -0.102910
   19  H  3  S    0.001662   0.401361  -0.050078  -0.060858  -0.009821
   20  H  3  S    0.000899   0.173534   2.168485  -0.375435   0.065272
   21  H  3  S    0.000559   0.007174  -4.326998   0.736068  -0.203641
   22  H  3  X    0.000000  -0.000000  -0.000000  -0.000000   0.000000
   23  H  3  Y   -0.000164  -0.021696  -0.001837   0.013136  -0.001626
   24  H  3  Z    0.000137  -0.000001  -0.000054   0.000515  -0.003718
   25  H  3  X    0.000000  -0.000000  -0.000000  -0.000000   0.000000
   26  H  3  Y   -0.000768  -0.005097  -0.543254  -0.071211   0.055957
   27  H  3  Z   -0.001100   0.000006   0.001076   0.006755  -0.102910

                      6          7          8          9         10
                    0.2298     0.2829     0.3020     0.4213     0.4366
                     B2         B1         A1         A1         B2  
    1  HE 1  S   -0.000000  -0.000000  -0.034259  -0.003623   0.000000
    2  HE 1  S   -0.000000  -0.000000  -0.291087   0.018503   0.000000
    3  HE 1  S   -0.000000  -0.000000   0.694407  -0.179706   0.000000
    4  HE 1  X   -0.000000  -0.003221   0.000000  -0.000000   0.000000
    5  HE 1  Y    0.000731  -0.000000   0.000000  -0.000000   0.011723
    6  HE 1  Z   -0.000000  -0.000000   0.006510   0.002876   0.000000
    7  HE 1  X   -0.000000   0.066186   0.000000  -0.000000   0.000000
    8  HE 1  Y   -0.015908  -0.000000   0.000000  -0.000000  -0.407550
    9  HE 1  Z   -0.000000  -0.000000  -0.245083  -0.010570   0.000000
   10  H  2  S    0.167264  -0.000000   0.030861   0.188676  -0.018450
   11  H  2  S   10.021306  -0.000000  -0.031001  -0.849023   2.570370
   12  H  2  S   -3.148521  -0.000000  -0.101684   0.426417  -0.361006
   13  H  2  X   -0.000000  -0.012420   0.000000  -0.000000   0.000000
   14  H  2  Y   -0.040386  -0.000000  -0.001202  -0.017091  -0.008040
   15  H  2  Z    0.000119  -0.000000  -0.012037  -0.000555   0.042398
   16  H  2  X   -0.000000   0.521289   0.000000  -0.000000   0.000000
   17  H  2  Y    2.104710  -0.000000   0.067388   1.266625   0.696621
   18  H  2  Z   -0.003627  -0.000000   0.493622  -0.063005  -1.685494
   19  H  3  S   -0.167264  -0.000000   0.030861   0.188676   0.018450
   20  H  3  S  -10.021306  -0.000000  -0.031001  -0.849023  -2.570370
   21  H  3  S    3.148521  -0.000000  -0.101684   0.426417   0.361006
   22  H  3  X   -0.000000  -0.012420   0.000000  -0.000000   0.000000
   23  H  3  Y   -0.040386  -0.000000   0.001202   0.017091  -0.008040
   24  H  3  Z   -0.000119  -0.000000  -0.012037  -0.000555  -0.042398
   25  H  3  X   -0.000000   0.521289   0.000000  -0.000000   0.000000
   26  H  3  Y    2.104710  -0.000000  -0.067388  -1.266625   0.696621
   27  H  3  Z    0.003627  -0.000000   0.493622  -0.063005   1.685494

                     11         12         13         14         15
                    0.4572     0.5323     0.5337     0.5700     0.6000
                     A2         B2         B1         A1         B2  
    1  HE 1  S    0.000000  -0.000000  -0.000000  -0.033165   0.000000
    2  HE 1  S    0.000000  -0.000000  -0.000000  -0.086962   0.000000
    3  HE 1  S    0.000000  -0.000000  -0.000000   0.301390   0.000000
    4  HE 1  X    0.000000  -0.000000  -0.022064  -0.000000   0.000000
    5  HE 1  Y    0.000000   0.009831  -0.000000  -0.000000  -0.013499
    6  HE 1  Z    0.000000  -0.000000  -0.000000  -0.016715   0.000000
    7  HE 1  X    0.000000  -0.000000   1.008420  -0.000000   0.000000
    8  HE 1  Y    0.000000  -0.551355  -0.000000  -0.000000   0.759771
    9  HE 1  Z    0.000000  -0.000000  -0.000000   0.899819   0.000000
   10  H  2  S    0.000000  -0.160953  -0.000000  -0.307159  -0.175191
   11  H  2  S    0.000000  22.296970  -0.000000   0.283924  20.815397
   12  H  2  S    0.000000  -1.911125  -0.000000  -0.192633  -1.269433
   13  H  2  X   -0.045823  -0.000000  -0.005917  -0.000000   0.000000
   14  H  2  Y    0.000000  -0.056908  -0.000000  -0.002933  -0.034864
   15  H  2  Z    0.000000  -0.013201  -0.000000   0.017344   0.006693
   16  H  2  X    1.938996  -0.000000  -0.049598  -0.000000   0.000000
   17  H  2  Y    0.000000   5.992050  -0.000000   0.302118   5.566140
   18  H  2  Z    0.000000   0.759439  -0.000000   0.221224  -0.613661
   19  H  3  S    0.000000   0.160953  -0.000000  -0.307159   0.175191
   20  H  3  S    0.000000 -22.296970  -0.000000   0.283924 -20.815397
   21  H  3  S    0.000000   1.911125  -0.000000  -0.192633   1.269433
   22  H  3  X    0.045823  -0.000000  -0.005917  -0.000000   0.000000
   23  H  3  Y    0.000000  -0.056908  -0.000000   0.002933  -0.034864
   24  H  3  Z    0.000000   0.013201  -0.000000   0.017344  -0.006693
   25  H  3  X   -1.938996  -0.000000  -0.049598  -0.000000   0.000000
   26  H  3  Y    0.000000   5.992050  -0.000000  -0.302118   5.566140
   27  H  3  Z    0.000000  -0.759439  -0.000000   0.221224   0.613661

                     16         17         18         19         20
                    0.7632     1.2011     1.6048     1.6470     1.7393
                     A1         B2         B1         A1         A1  
    1  HE 1  S    0.000519  -0.000000   0.000000   0.026668  -1.226086
    2  HE 1  S    0.088750  -0.000000   0.000000   0.028376   1.718415
    3  HE 1  S   -0.187951  -0.000000   0.000000  -0.187392  -0.799222
    4  HE 1  X    0.000000  -0.000000   0.008508   0.000000  -0.000000
    5  HE 1  Y    0.000000   0.004174   0.000000   0.000000  -0.000000
    6  HE 1  Z   -0.003568  -0.000000   0.000000  -0.037713  -0.015068
    7  HE 1  X    0.000000  -0.000000   0.025786   0.000000  -0.000000
    8  HE 1  Y    0.000000   0.088428   0.000000   0.000000  -0.000000
    9  HE 1  Z   -0.464544  -0.000000   0.000000  -0.149214  -0.042403
   10  H  2  S   -0.773676   0.916796   0.000000  -0.027199   0.004211
   11  H  2  S    0.613130  20.174977   0.000000  -0.010967  -0.004482
   12  H  2  S   -0.127049  -0.318909   0.000000   0.056964   0.068303
   13  H  2  X    0.000000  -0.000000   0.675844   0.000000  -0.000000
   14  H  2  Y   -0.020656  -0.487438   0.000000   0.004136   0.001237
   15  H  2  Z   -0.021093  -0.002667   0.000000   0.683700   0.022195
   16  H  2  X    0.000000  -0.000000  -0.293271   0.000000  -0.000000
   17  H  2  Y    0.829062   6.038697   0.000000   0.058913   0.000260
   18  H  2  Z   -0.090738  -0.025895   0.000000  -0.341348  -0.054626
   19  H  3  S   -0.773676  -0.916796   0.000000  -0.027199   0.004211
   20  H  3  S    0.613130 -20.174977   0.000000  -0.010967  -0.004482
   21  H  3  S   -0.127049   0.318909   0.000000   0.056964   0.068303
   22  H  3  X    0.000000  -0.000000   0.675844   0.000000  -0.000000
   23  H  3  Y    0.020656  -0.487438   0.000000  -0.004136  -0.001237
   24  H  3  Z   -0.021093   0.002667   0.000000   0.683700   0.022195
   25  H  3  X    0.000000  -0.000000  -0.293271   0.000000  -0.000000
   26  H  3  Y   -0.829062   6.038697   0.000000  -0.058913  -0.000260
   27  H  3  Z   -0.090738   0.025895   0.000000  -0.341348  -0.054626

                     21         22         23         24         25
                    2.0622     2.2270     2.2295     3.0271     3.0399
                     A1         A2         B2         B1         B2  
    1  HE 1  S    0.001119  -0.000000   0.000000   0.000000  -0.000000
    2  HE 1  S    0.001012  -0.000000   0.000000   0.000000  -0.000000
    3  HE 1  S    0.037254  -0.000000   0.000000   0.000000  -0.000000
    4  HE 1  X    0.000000  -0.000000   0.000000   1.131128  -0.000000
    5  HE 1  Y    0.000000  -0.000000   0.033483   0.000000  -1.132974
    6  HE 1  Z   -0.010978  -0.000000   0.000000   0.000000  -0.000000
    7  HE 1  X    0.000000  -0.000000   0.000000  -0.510042  -0.000000
    8  HE 1  Y    0.000000  -0.000000   0.029780   0.000000   0.527061
    9  HE 1  Z   -0.005671  -0.000000   0.000000   0.000000  -0.000000
   10  H  2  S   -0.395165  -0.000000   0.005671   0.000000   0.073239
   11  H  2  S    0.388528  -0.000000   0.201945   0.000000   1.304203
   12  H  2  S   -0.143052  -0.000000   0.005763   0.000000   0.073998
   13  H  2  X    0.000000   1.044824   0.000000  -0.010971  -0.000000
   14  H  2  Y    0.749385  -0.000000  -0.000245   0.000000   0.035762
   15  H  2  Z   -0.008031  -0.000000   1.044720   0.000000   0.047624
   16  H  2  X    0.000000  -0.695983   0.000000   0.013363  -0.000000
   17  H  2  Y   -0.435805  -0.000000   0.056909   0.000000   0.361711
   18  H  2  Z    0.008294  -0.000000  -0.704362   0.000000  -0.117672
   19  H  3  S   -0.395165  -0.000000  -0.005671   0.000000  -0.073239
   20  H  3  S    0.388528  -0.000000  -0.201945   0.000000  -1.304203
   21  H  3  S   -0.143052  -0.000000  -0.005763   0.000000  -0.073998
   22  H  3  X    0.000000  -1.044824   0.000000  -0.010971  -0.000000
   23  H  3  Y   -0.749385  -0.000000  -0.000245   0.000000   0.035762
   24  H  3  Z   -0.008031  -0.000000  -1.044720   0.000000  -0.047624
   25  H  3  X    0.000000   0.695983   0.000000   0.013363  -0.000000
   26  H  3  Y    0.435805  -0.000000   0.056909   0.000000   0.361711
   27  H  3  Z    0.008294  -0.000000   0.704362   0.000000   0.117672

                     26         27
                    3.0705     3.5679
                     A1         B2  
    1  HE 1  S   -0.028690  -0.000000
    2  HE 1  S    0.078907  -0.000000
    3  HE 1  S   -0.117459  -0.000000
    4  HE 1  X    0.000000  -0.000000
    5  HE 1  Y    0.000000  -0.034255
    6  HE 1  Z    1.138919  -0.000000
    7  HE 1  X    0.000000  -0.000000
    8  HE 1  Y    0.000000   0.037798
    9  HE 1  Z   -0.565041  -0.000000
   10  H  2  S   -0.050539  -2.371819
   11  H  2  S    0.029435   7.476697
   12  H  2  S    0.032646  -0.376190
   13  H  2  X    0.000000  -0.000000
   14  H  2  Y    0.009705  -2.075830
   15  H  2  Z    0.056305   0.003175
   16  H  2  X    0.000000  -0.000000
   17  H  2  Y    0.054232   2.113481
   18  H  2  Z   -0.081618  -0.009568
   19  H  3  S   -0.050539   2.371819
   20  H  3  S    0.029435  -7.476697
   21  H  3  S    0.032646   0.376190
   22  H  3  X    0.000000  -0.000000
   23  H  3  Y   -0.009705  -2.075830
   24  H  3  Z    0.056305  -0.003175
   25  H  3  X    0.000000  -0.000000
   26  H  3  Y   -0.054232   2.113481
   27  H  3  Z   -0.081618   0.009568
 ...... END OF RHF CALCULATION ......
 STEP CPU TIME =     0.01 TOTAL CPU TIME =        0.1 (    0.0 MIN)
 TOTAL WALL CLOCK TIME=        0.8 SECONDS, CPU UTILIZATION IS   7.89%

     ----------------------------------------------------------------
     PROPERTY VALUES FOR THE RHF   SELF-CONSISTENT FIELD WAVEFUNCTION
     ----------------------------------------------------------------

          -----------------
          ENERGY COMPONENTS
          -----------------

         WAVEFUNCTION NORMALIZATION =       1.0000000000

                ONE ELECTRON ENERGY =      -7.6711052217
                TWO ELECTRON ENERGY =       2.3347610066
           NUCLEAR REPULSION ENERGY =       1.3521176085
                                      ------------------
                       TOTAL ENERGY =      -3.9842266066

 ELECTRON-ELECTRON POTENTIAL ENERGY =       2.3347610066
  NUCLEUS-ELECTRON POTENTIAL ENERGY =     -11.5902333599
   NUCLEUS-NUCLEUS POTENTIAL ENERGY =       1.3521176085
                                      ------------------
             TOTAL POTENTIAL ENERGY =      -7.9033547447
               TOTAL KINETIC ENERGY =       3.9191281381
                 VIRIAL RATIO (V/T) =       2.0166104465

  ...... PI ENERGY ANALYSIS ......

 ENERGY ANALYSIS:
            FOCK ENERGY=     -3.0015832037
          BARE H ENERGY=     -7.6711052217
     ELECTRONIC ENERGY =     -5.3363442127
         KINETIC ENERGY=      3.9191281381
          N-N REPULSION=      1.3521176085
           TOTAL ENERGY=     -3.9842266042
        SIGMA PART(1+2)=     -5.3363442127
               (K,V1,2)=      3.9191281381    -11.5902333599      2.3347610090
           PI PART(1+2)=      0.0000000000
               (K,V1,2)=      0.0000000000      0.0000000000      0.0000000000
  SIGMA SKELETON, ERROR=     -3.9842266042      0.0000000000
             MIXED PART= 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00
 ...... END OF PI ENERGY ANALYSIS ......

          ---------------------------------------
          MULLIKEN AND LOWDIN POPULATION ANALYSES
          ---------------------------------------

               ----- POPULATIONS IN EACH AO -----
                             MULLIKEN      LOWDIN
              1  HE 1  S      1.08896     1.05588
              2  HE 1  S      0.88489     0.80662
              3  HE 1  S      0.02540     0.13312
              4  HE 1  X      0.00000     0.00000
              5  HE 1  Y      0.00000     0.00000
              6  HE 1  Z      0.00000     0.00000
              7  HE 1  X      0.00000     0.00000
              8  HE 1  Y      0.00000     0.00000
              9  HE 1  Z      0.00001     0.00006
             10  H  2  S      0.68685     0.56072
             11  H  2  S      0.29160     0.32811
             12  H  2  S      0.00715     0.05159
             13  H  2  X      0.00000     0.00000
             14  H  2  Y      0.01175     0.02774
             15  H  2  Z     -0.00000     0.00000
             16  H  2  X      0.00000     0.00000
             17  H  2  Y      0.00289     0.03362
             18  H  2  Z      0.00013     0.00039
             19  H  3  S      0.68685     0.56072
             20  H  3  S      0.29160     0.32811
             21  H  3  S      0.00715     0.05159
             22  H  3  X      0.00000     0.00000
             23  H  3  Y      0.01175     0.02774
             24  H  3  Z     -0.00000     0.00000
             25  H  3  X      0.00000     0.00000
             26  H  3  Y      0.00289     0.03362
             27  H  3  Z      0.00013     0.00039

          ----- MULLIKEN ATOMIC OVERLAP POPULATIONS -----
          (OFF-DIAGONAL ELEMENTS NEED TO BE MULTIPLIED BY 2)

             1           2           3

    1    1.9992802
    2   -0.0000095   0.5815530
    3   -0.0000095   0.4188259   0.5815530

          TOTAL MULLIKEN AND LOWDIN ATOMIC POPULATIONS
       ATOM         MULL.POP.    CHARGE          LOW.POP.     CHARGE
    1 HE            1.999261    0.000739         1.995675    0.004325
    2 H             1.000369   -0.000369         1.002162   -0.002162
    3 H             1.000369   -0.000369         1.002162   -0.002162

          MULLIKEN SPHERICAL HARMONIC POPULATIONS
       ATOM           S       P       D      F      G      H      I    TOTAL
    1 HE            2.00    0.00    0.00   0.00   0.00   0.00   0.00    2.00
    2 H             0.99    0.01    0.00   0.00   0.00   0.00   0.00    1.00
    3 H             0.99    0.01    0.00   0.00   0.00   0.00   0.00    1.00

          ---------------------
          ELECTROSTATIC MOMENTS
          ---------------------

 POINT   1           X           Y           Z (BOHR)    CHARGE
                -0.000000    0.000000   -0.990462        0.00 (A.U.)
         DX          DY          DZ         /D/  (DEBYE)
     0.000000    0.000000   -0.002899    0.002899
 ...... END OF PROPERTY EVALUATION ......
 STEP CPU TIME =     0.00 TOTAL CPU TIME =        0.1 (    0.0 MIN)
 TOTAL WALL CLOCK TIME=        0.8 SECONDS, CPU UTILIZATION IS   7.69%

 OCCUPATIONALLY RESTRICTED MULTIPLE ACTIVE SPACE CI OPTIONS
                                    NRNFG   NPFLG
 INTEGRAL TRANSFORMATION              1       0
 DIAGONALIZE HAMILTONIAN              1       0
 FORM 1E- DENSITY MATRIX              1       0
 FORM 2E- DENSITY MATRIX              0       0

     ------------------------------------------------------------
       DIRECT DETERMINANT ORMAS-CI INPUT SORTER
       PROGRAM WRITTEN BY JOE IVANIC AND MIKE SCHMIDT
       ORMAS = OCCUPATION RESTRICTED MULTIPLE ACTIVE SPACE
     ------------------------------------------------------------

 THE POINT GROUP                  =   C2V     
 THE STATE SYMMETRY               =   A1      
 NUMBER OF CORE ORBITALS          =    0
 NUMBER OF ACTIVE ORBITALS        =   27
 NUMBER OF ALPHA ELECTRONS        =    2 (   2 ACTIVE)
 NUMBER OF BETA ELECTRONS         =    2 (   2 ACTIVE)
 NUMBER OF OCCUPIED ORBITALS      =   27

 NUMBER OF CI STATES REQUESTED    =    1
 NUMBER OF CI STARTING VECTORS    =    1
 MAX. NO. OF CI EXPANSION VECTORS =   10
 SIZE OF INITIAL CI GUESS MATRIX  =  300

 MAX. NO. OF CI ITERS/STATE       =  100
 CI DIAGONALIZATION CRITERION     = 1.00E-06
 CI PROPERTIES WILL BE FOUND FOR ROOT NUMBER   1
 1E- DENSITY MATRIX OPTIONS (NFLGDM) ARE 1

 CORRELATION ENERGY ANALYSIS      =    T

 FULLY DIRECT OPTION              =    T

 CALC. OF DAVIDSON Q CORRECTION   =    F

 BLOCK SPACE MIXING IN NO DIAG    =    F

 SYMMETRIES FOR THE   0 CORE,  27 ACTIVE ARE

   ACTIVE= A1    A1    B2    A1    A1    B2    B1    A1    A1    B2  
           A2    B2    B1    A1    B2    A1    B2    B1    A1    A1  
           A1    A2    B2    B1    B2    A1    B2  

 THE NUMBER OF SPACES             =   2
 EACH SPACE STARTS AT ORBITAL     =   1   3
 NO OF ORBITALS IN EACH SPACE     =   2  25
 MIN NO OF ELECS IN EACH SPACE    =   0   0
 MAX NO OF ELECS IN EACH SPACE    =   4   4

 SPACE SPECIFICATIONS HAVE PASSED PRELIMINARY CHECKS

 ALL VALID DETERMINANTS WILL BE EXPRESSED AS PAIRS OF ALPHA AND BETA STRINGS.

 MIN NO OF ALPHA ELECS =    0   0
 MAX NO OF ALPHA ELECS =    2   2

 MIN NO OF BETA  ELECS =    0   0
 MAX NO OF BETA  ELECS =    2   2

     --------------------------------------------
     PARTIAL TWO ELECTRON INTEGRAL TRANSFORMATION
     --------------------------------------------

 NUMBER OF CORE MOLECULAR ORBITALS     =    0
 NUMBER OF OCCUPIED MOLECULAR ORBITALS =   27
 TOTAL NUMBER OF MOLECULAR ORBITALS    =   27
 TOTAL NUMBER OF ATOMIC ORBITALS       =   27
 THRESHOLD FOR KEEPING TRANSFORMED 2E- INTEGRALS = 1.000E-10
 AO INTEGRALS WILL BE READ IN FROM DISK...

 PLAN A: REQUIREMENTS FOR FULLY IN-MEMORY TRANSFORMATION:
 # OF WORDS AVAILABLE =            100000000
 # OF WORDS NEEDED    =               339263

 CHOOSING IN MEMORY PARTIAL TRANSFORMATION...
 TOTAL NUMBER OF TRANSFORMED 2E- INTEGRALS KEPT =        20143
 ... END OF INTEGRAL TRANSFORMATION ...
 STEP CPU TIME =     0.02 TOTAL CPU TIME =        0.1 (    0.0 MIN)
 TOTAL WALL CLOCK TIME=        0.8 SECONDS, CPU UTILIZATION IS   9.64%

     --------------------------------------------------
                 DIRECT DETERMINANT ORMAS-CI 
                PROGRAM WRITTEN BY JOE IVANIC
     --------------------------------------------------

 STORAGE OF BINOMIAL COEFFICIENTS REQUIRES          87 WORDS

 TOTAL NUMBER OF ALPHA GROUPS  =            3
 TOTAL NUMBER OF BETA  GROUPS  =            3

 TOTAL NUMBER OF ALPHA STRINGS =          351
 TOTAL NUMBER OF BETA  STRINGS =          351

 STORAGE OF TABLES REQUIRES                      2362 WORDS

 TOTAL NUMBER OF ALPHA-BETA GROUP COMBINATIONS =            9
 OF THESE THE ALLOWED NUMBER OF COMBINATIONS   =            9

 TIME FOR SETTING UP TABLE SET 1 :          0.0

 THE NUMBER OF DETERMINANTS HAVING SPACE SYMMETRY A1 
 IN POINT GROUP C2V  WITH SZ=  0.0 IS          33293

 INTEGRAL STORAGE REQUIRES                  102388 WORDS
 EXTRA ORMAS STORAGE REQUIRES               911744 WORDS
 (EXTRA ORMAS STORAGE INCLUDES THAT FOR MXXPAN = 10)
 TOTAL ORMAS CALCULATION REQUIRES          1016581 WORDS

 INITIAL ORMAS VECTOR GUESS TIME :          0.1
 ORMAS SPIN CALCULATION TIME     :          0.0
 INITIAL ORMAS CI ITERATION TIME :          0.2

 ITERATION      ENERGY           GRADIENT
     0       -3.9997154249     0.46805803
     1       -4.0534714218     0.06634368
     2       -4.0544327169     0.01229337
     3       -4.0544706618     0.00234377
     4       -4.0544722832     0.00048946
     5       -4.0544723516     0.00009713
     6       -4.0544723545     0.00002012
     7       -4.0544723546     0.00000383
     8       -4.0544723546     0.00000085

 CONVERGED STATE    1 ENERGY=       -4.0544723546 IN    8 ITERS

 ALL STATES CONVERGED.

 CI EIGENVECTORS WILL BE LABELED IN GROUP=C2V     
 PRINTING CI COEFFICIENTS LARGER THAN  0.050000

 STATE   1  ENERGY=        -4.0544723546  S=  0.00  SZ=  0.00  SPACE SYM=A1  

   ALPHA     |    BETA     | COEFFICIENT
 1     2     | 1     2     |
-------------|-------------|------------
  1  2       |  1  2       |   0.9864624

 ..... DONE WITH ORMAS-CI COMPUTATION .....
 STEP CPU TIME =     1.92 TOTAL CPU TIME =        2.0 (    0.0 MIN)
 TOTAL WALL CLOCK TIME=        2.8 SECONDS, CPU UTILIZATION IS  72.46%

     ---------------------------
     ONE PARTICLE DENSITY MATRIX
     ---------------------------

 DENSITY MATRIX WILL BE SAVED FOR   1-TH STATE WITH S= 0.00

 CI EIGENSTATE   1 TOTAL ENERGY =       -4.0544723546

          NATURAL ORBITALS IN ATOMIC ORBITAL BASIS
          ----------------------------------------

                      1          2          3          4          5
                    1.9834     1.9627     0.0219     0.0086     0.0059
                     A1         A1         B2         A1         A1  
    1  HE 1  S    0.592711  -0.008415  -0.000000   1.130305   0.141456
    2  HE 1  S    0.495403  -0.009298  -0.000000  -1.129711  -0.164648
    3  HE 1  S    0.028776  -0.001157  -0.000000  -0.039574  -0.017585
    4  HE 1  X    0.000000  -0.000000  -0.000000   0.000000  -0.000000
    5  HE 1  Y    0.000000  -0.000000  -0.000197   0.000000  -0.000000
    6  HE 1  Z    0.000024   0.000282  -0.000000  -0.001760   0.020410
    7  HE 1  X    0.000000  -0.000000  -0.000000   0.000000  -0.000000
    8  HE 1  Y    0.000000  -0.000000   0.000580   0.000000  -0.000000
    9  HE 1  Z   -0.000520   0.000325  -0.000000   0.005559  -0.012487
   10  H  2  S    0.001219   0.407254   0.946887   0.063167  -0.556256
   11  H  2  S    0.000693   0.170222   0.307722  -0.055766   0.631129
   12  H  2  S    0.000335   0.006463  -0.032808  -0.002253  -0.028197
   13  H  2  X    0.000000  -0.000000  -0.000000   0.000000  -0.000000
   14  H  2  Y    0.000127   0.018953  -0.025817   0.024195  -0.260949
   15  H  2  Z    0.000118  -0.000010   0.000004   0.010319   0.005662
   16  H  2  X    0.000000  -0.000000  -0.000000   0.000000  -0.000000
   17  H  2  Y    0.000761   0.000864   0.014054  -0.015021   0.119104
   18  H  2  Z   -0.001020   0.000011  -0.000219   0.011234   0.000001
   19  H  3  S    0.001219   0.407254  -0.946887   0.063167  -0.556256
   20  H  3  S    0.000693   0.170222  -0.307722  -0.055766   0.631129
   21  H  3  S    0.000335   0.006463   0.032808  -0.002253  -0.028197
   22  H  3  X    0.000000  -0.000000  -0.000000   0.000000  -0.000000
   23  H  3  Y   -0.000127  -0.018953  -0.025817  -0.024195   0.260949
   24  H  3  Z    0.000118  -0.000010  -0.000004   0.010319   0.005662
   25  H  3  X    0.000000  -0.000000  -0.000000   0.000000  -0.000000
   26  H  3  Y   -0.000761  -0.000864   0.014054   0.015021  -0.119104
   27  H  3  Z   -0.001020   0.000011   0.000219   0.011234   0.000001

                      6          7          8          9         10
                    0.0043     0.0043     0.0026     0.0026     0.0026
                     A1         B1         A1         B2         B1  
    1  HE 1  S   -0.052719  -0.000000  -0.001940   0.000000  -0.000000
    2  HE 1  S    0.057895  -0.000000  -0.003766   0.000000  -0.000000
    3  HE 1  S    0.041764  -0.000000   0.019070   0.000000  -0.000000
    4  HE 1  X   -0.000000   0.001593  -0.000000   0.000000   0.844989
    5  HE 1  Y   -0.000000  -0.000000  -0.000000   0.843531  -0.000000
    6  HE 1  Z   -0.015321  -0.000000   0.839493   0.000000  -0.000000
    7  HE 1  X   -0.000000  -0.003685  -0.000000   0.000000   0.269904
    8  HE 1  Y   -0.000000  -0.000000  -0.000000   0.273252  -0.000000
    9  HE 1  Z    0.019354  -0.000000   0.282192   0.000000  -0.000000
   10  H  2  S    0.006211  -0.000000   0.024238   0.002131  -0.000000
   11  H  2  S   -0.000661  -0.000000  -0.030090   0.212235  -0.000000
   12  H  2  S   -0.011420  -0.000000  -0.005969   0.018222  -0.000000
   13  H  2  X   -0.000000   0.465899  -0.000000   0.000000  -0.000157
   14  H  2  Y    0.002699  -0.000000   0.011217  -0.002324  -0.000000
   15  H  2  Z    0.461938  -0.000000   0.006417   0.001986  -0.000000
   16  H  2  X   -0.000000   0.167224  -0.000000   0.000000  -0.002067
   17  H  2  Y   -0.010392  -0.000000  -0.013729   0.058463  -0.000000
   18  H  2  Z    0.175263  -0.000000   0.015271  -0.015138  -0.000000
   19  H  3  S    0.006211  -0.000000   0.024238  -0.002131  -0.000000
   20  H  3  S   -0.000661  -0.000000  -0.030090  -0.212235  -0.000000
   21  H  3  S   -0.011420  -0.000000  -0.005969  -0.018222  -0.000000
   22  H  3  X   -0.000000   0.465899  -0.000000   0.000000  -0.000157
   23  H  3  Y   -0.002699  -0.000000  -0.011217  -0.002324  -0.000000
   24  H  3  Z    0.461938  -0.000000   0.006417  -0.001986  -0.000000
   25  H  3  X   -0.000000   0.167224  -0.000000   0.000000  -0.002067
   26  H  3  Y    0.010392  -0.000000   0.013729   0.058463  -0.000000
   27  H  3  Z    0.175263  -0.000000   0.015271   0.015138  -0.000000

                     11         12         13         14         15
                    0.0002     0.0002     0.0002     0.0002     0.0001
                     A1         B2         B2         A2         A1  
    1  HE 1  S   -0.023539  -0.000000   0.000000  -0.000000   0.055245
    2  HE 1  S    0.036403  -0.000000   0.000000  -0.000000  -0.104042
    3  HE 1  S    0.004556  -0.000000   0.000000  -0.000000  -0.007161
    4  HE 1  X    0.000000  -0.000000   0.000000  -0.000000   0.000000
    5  HE 1  Y    0.000000  -0.013728  -0.004590  -0.000000   0.000000
    6  HE 1  Z    0.019612  -0.000000   0.000000  -0.000000  -0.470208
    7  HE 1  X    0.000000  -0.000000   0.000000  -0.000000   0.000000
    8  HE 1  Y    0.000000   0.025849   0.001931  -0.000000   0.000000
    9  HE 1  Z    0.007490  -0.000000   0.000000  -0.000000   0.618009
   10  H  2  S    0.707761   1.826775  -0.109768  -0.000000   0.051741
   11  H  2  S   -0.462865   5.829953  -0.460474  -0.000000  -0.070379
   12  H  2  S    0.059896   0.240664  -0.009716  -0.000000   0.008616
   13  H  2  X    0.000000  -0.000000   0.000000   0.931344   0.000000
   14  H  2  Y   -0.621197   0.879629  -0.055386  -0.000000  -0.031426
   15  H  2  Z   -0.007595   0.048838   0.926515  -0.000000   0.366505
   16  H  2  X    0.000000  -0.000000   0.000000   0.172511   0.000000
   17  H  2  Y   -0.289475   2.157296  -0.155899  -0.000000   0.029883
   18  H  2  Z    0.009491   0.026050   0.186254  -0.000000  -0.408299
   19  H  3  S    0.707761  -1.826775   0.109768  -0.000000   0.051741
   20  H  3  S   -0.462865  -5.829953   0.460474  -0.000000  -0.070379
   21  H  3  S    0.059896  -0.240664   0.009716  -0.000000   0.008616
   22  H  3  X    0.000000  -0.000000   0.000000  -0.931344   0.000000
   23  H  3  Y    0.621197   0.879629  -0.055386  -0.000000   0.031426
   24  H  3  Z   -0.007595  -0.048838  -0.926515  -0.000000   0.366505
   25  H  3  X    0.000000  -0.000000   0.000000  -0.172511   0.000000
   26  H  3  Y    0.289475   2.157296  -0.155899  -0.000000  -0.029883
   27  H  3  Z    0.009491  -0.026050  -0.186254  -0.000000  -0.408299

                     16         17         18         19         20
                    0.0001     0.0001     0.0001     0.0001     0.0000
                     B1         B2         A1         B1         A1  
    1  HE 1  S   -0.000000  -0.000000  -0.063764  -0.000000   0.372052
    2  HE 1  S   -0.000000  -0.000000   0.112584  -0.000000  -1.139846
    3  HE 1  S   -0.000000  -0.000000   0.039594  -0.000000   1.297667
    4  HE 1  X   -0.188887  -0.000000   0.000000  -0.728239   0.000000
    5  HE 1  Y   -0.000000  -0.754577   0.000000  -0.000000   0.000000
    6  HE 1  Z   -0.000000  -0.000000  -0.606827  -0.000000  -0.034441
    7  HE 1  X    0.261002  -0.000000   0.000000   1.068230   0.000000
    8  HE 1  Y   -0.000000   1.104327   0.000000  -0.000000   0.000000
    9  HE 1  Z   -0.000000  -0.000000   0.967846  -0.000000   0.091896
   10  H  2  S   -0.000000  -0.044887   0.074180  -0.000000   0.146794
   11  H  2  S   -0.000000   0.405047  -0.032608  -0.000000   0.010504
   12  H  2  S   -0.000000   0.101734  -0.035427  -0.000000  -0.000562
   13  H  2  X   -0.471762  -0.000000   0.000000   0.132111   0.000000
   14  H  2  Y   -0.000000  -0.043673   0.010458  -0.000000   0.171001
   15  H  2  Z   -0.000000  -0.012611  -0.339156  -0.000000  -0.045088
   16  H  2  X    0.553311  -0.000000   0.000000  -0.162074   0.000000
   17  H  2  Y   -0.000000   0.110206  -0.153513  -0.000000  -0.870481
   18  H  2  Z   -0.000000  -0.018061   0.441183  -0.000000   0.132646
   19  H  3  S   -0.000000   0.044887   0.074180  -0.000000   0.146794
   20  H  3  S   -0.000000  -0.405047  -0.032608  -0.000000   0.010504
   21  H  3  S   -0.000000  -0.101734  -0.035427  -0.000000  -0.000562
   22  H  3  X   -0.471762  -0.000000   0.000000   0.132111   0.000000
   23  H  3  Y   -0.000000  -0.043673  -0.010458  -0.000000  -0.171001
   24  H  3  Z   -0.000000   0.012611  -0.339156  -0.000000  -0.045088
   25  H  3  X    0.553311  -0.000000   0.000000  -0.162074   0.000000
   26  H  3  Y   -0.000000   0.110206   0.153513  -0.000000   0.870481
   27  H  3  Z   -0.000000   0.018061   0.441183  -0.000000   0.132646

                     21         22         23         24         25
                    0.0000     0.0000     0.0000     0.0000     0.0000
                     A1         B2         B2         A2         B2  
    1  HE 1  S    0.288971   0.000000  -0.000000   0.000000   0.000000
    2  HE 1  S   -0.905488   0.000000  -0.000000   0.000000   0.000000
    3  HE 1  S    1.114432   0.000000  -0.000000   0.000000   0.000000
    4  HE 1  X    0.000000   0.000000  -0.000000   0.000000   0.000000
    5  HE 1  Y    0.000000   0.018306  -0.060495   0.000000  -0.033488
    6  HE 1  Z   -0.049550   0.000000  -0.000000   0.000000   0.000000
    7  HE 1  X    0.000000   0.000000  -0.000000   0.000000   0.000000
    8  HE 1  Y    0.000000  -0.030690   0.158490   0.000000   0.112462
    9  HE 1  Z    0.110577   0.000000  -0.000000   0.000000   0.000000
   10  H  2  S   -0.069724  -1.492054  -0.234276   0.000000  -0.166796
   11  H  2  S   -0.157717   4.826172  10.265650   0.000000  32.705207
   12  H  2  S   -0.192179   0.004054  -0.184762   0.000000  -0.720950
   13  H  2  X    0.000000   0.000000  -0.000000  -0.475768   0.000000
   14  H  2  Y   -0.235216  -1.812618  -0.403221   0.000000  -0.564030
   15  H  2  Z   -0.081075  -0.066540   0.463455   0.000000  -0.124678
   16  H  2  X    0.000000   0.000000  -0.000000   2.052885   0.000000
   17  H  2  Y    1.040987   1.801245   2.876079   0.000000   8.934682
   18  H  2  Z    0.166264   0.266588  -1.975758   0.000000   0.542508
   19  H  3  S   -0.069724   1.492054   0.234276   0.000000   0.166796
   20  H  3  S   -0.157717  -4.826172 -10.265650   0.000000 -32.705207
   21  H  3  S   -0.192179  -0.004054   0.184762   0.000000   0.720950
   22  H  3  X    0.000000   0.000000  -0.000000   0.475768   0.000000
   23  H  3  Y    0.235216  -1.812618  -0.403221   0.000000  -0.564030
   24  H  3  Z   -0.081075   0.066540  -0.463455   0.000000   0.124678
   25  H  3  X    0.000000   0.000000  -0.000000  -2.052885   0.000000
   26  H  3  Y   -1.040987   1.801245   2.876079   0.000000   8.934682
   27  H  3  Z    0.166264  -0.266588   1.975758   0.000000  -0.542508

                     26         27
                    0.0000     0.0000
                     A1         B2  
    1  HE 1  S    0.009093   0.000000
    2  HE 1  S   -0.093234   0.000000
    3  HE 1  S    0.489575   0.000000
    4  HE 1  X    0.000000   0.000000
    5  HE 1  Y    0.000000   0.010246
    6  HE 1  Z   -0.022467   0.000000
    7  HE 1  X    0.000000   0.000000
    8  HE 1  Y    0.000000  -0.044143
    9  HE 1  Z    0.100268   0.000000
   10  H  2  S   -0.196020  -0.046445
   11  H  2  S    0.910531  16.553957
   12  H  2  S   -0.902950  -5.801320
   13  H  2  X    0.000000   0.000000
   14  H  2  Y    0.149479  -0.091659
   15  H  2  Z   -0.022243  -0.000861
   16  H  2  X    0.000000   0.000000
   17  H  2  Y   -0.792282   4.138196
   18  H  2  Z    0.062596   0.006103
   19  H  3  S   -0.196020   0.046445
   20  H  3  S    0.910531 -16.553957
   21  H  3  S   -0.902950   5.801320
   22  H  3  X    0.000000   0.000000
   23  H  3  Y   -0.149479  -0.091659
   24  H  3  Z   -0.022243   0.000861
   25  H  3  X    0.000000   0.000000
   26  H  3  Y    0.792282   4.138196
   27  H  3  Z    0.062596  -0.006103
 ..... DONE WITH ONE PARTICLE DENSITY MATRIX .....
 STEP CPU TIME =     0.00 TOTAL CPU TIME =        2.0 (    0.0 MIN)
 TOTAL WALL CLOCK TIME=        2.8 SECONDS, CPU UTILIZATION IS  72.20%

     --------------------------------------------------------
     ORMAS CI PROPERTIES...FOR THE WAVEFUNCTION OF STATE    1
               USING THE EXPECTATION VALUE DENSITY
     --------------------------------------------------------

          -----------------
          ENERGY COMPONENTS
          -----------------

         WAVEFUNCTION NORMALIZATION =       1.0000000000

                ONE ELECTRON ENERGY =      -7.6031513994
                TWO ELECTRON ENERGY =       2.1965614362
           NUCLEAR REPULSION ENERGY =       1.3521176085
                                      ------------------
                       TOTAL ENERGY =      -4.0544723546

 ELECTRON-ELECTRON POTENTIAL ENERGY =       2.1965614362
  NUCLEUS-ELECTRON POTENTIAL ENERGY =     -11.5724441222
   NUCLEUS-NUCLEUS POTENTIAL ENERGY =       1.3521176085
                                      ------------------
             TOTAL POTENTIAL ENERGY =      -8.0237650775
               TOTAL KINETIC ENERGY =       3.9692927229
                 VIRIAL RATIO (V/T) =       2.0214596498

          ---------------------------------------
          MULLIKEN AND LOWDIN POPULATION ANALYSES
          ---------------------------------------

               ----- POPULATIONS IN EACH AO -----
                             MULLIKEN      LOWDIN
              1  HE 1  S      1.07962     1.04660
              2  HE 1  S      0.88088     0.80380
              3  HE 1  S      0.03093     0.13703
              4  HE 1  X      0.00213     0.00203
              5  HE 1  Y      0.00213     0.00204
              6  HE 1  Z      0.00213     0.00204
              7  HE 1  X      0.00051     0.00060
              8  HE 1  Y      0.00052     0.00061
              9  HE 1  Z      0.00055     0.00068
             10  H  2  S      0.69521     0.56285
             11  H  2  S      0.28260     0.32008
             12  H  2  S      0.00615     0.04997
             13  H  2  X      0.00172     0.00161
             14  H  2  Y      0.01116     0.03050
             15  H  2  Z      0.00172     0.00162
             16  H  2  X      0.00053     0.00064
             17  H  2  Y      0.00053     0.03394
             18  H  2  Z      0.00068     0.00107
             19  H  3  S      0.69521     0.56285
             20  H  3  S      0.28260     0.32008
             21  H  3  S      0.00615     0.04997
             22  H  3  X      0.00172     0.00161
             23  H  3  Y      0.01116     0.03050
             24  H  3  Z      0.00172     0.00162
             25  H  3  X      0.00053     0.00064
             26  H  3  Y      0.00053     0.03394
             27  H  3  Z      0.00068     0.00107

          ----- MULLIKEN ATOMIC OVERLAP POPULATIONS -----
          (OFF-DIAGONAL ELEMENTS NEED TO BE MULTIPLIED BY 2)

             1           2           3

    1    1.9995072
    2   -0.0000547   0.6229811
    3   -0.0000547   0.3773747   0.6229811

          TOTAL MULLIKEN AND LOWDIN ATOMIC POPULATIONS
       ATOM         MULL.POP.    CHARGE          LOW.POP.     CHARGE
    1 HE            1.999398    0.000602         1.995435    0.004565
    2 H             1.000301   -0.000301         1.002283   -0.002283
    3 H             1.000301   -0.000301         1.002283   -0.002283

          MULLIKEN SPHERICAL HARMONIC POPULATIONS
       ATOM           S       P       D      F      G      H      I    TOTAL
    1 HE            1.99    0.01    0.00   0.00   0.00   0.00   0.00    2.00
    2 H             0.98    0.02    0.00   0.00   0.00   0.00   0.00    1.00
    3 H             0.98    0.02    0.00   0.00   0.00   0.00   0.00    1.00

          ---------------------
          ELECTROSTATIC MOMENTS
          ---------------------

 POINT   1           X           Y           Z (BOHR)    CHARGE
                -0.000000    0.000000   -0.990462        0.00 (A.U.)
         DX          DY          DZ         /D/  (DEBYE)
     0.000000   -0.000000   -0.002284    0.002284
 ...... END OF PROPERTY EVALUATION ......
 STEP CPU TIME =     0.00 TOTAL CPU TIME =        2.0 (    0.0 MIN)
 TOTAL WALL CLOCK TIME=        2.8 SECONDS, CPU UTILIZATION IS  72.20%
              1019830  WORDS OF DYNAMIC MEMORY USED
 EXECUTION OF GAMESS TERMINATED NORMALLY Sun May 11 14:57:32 2014
 DDI: 263624 bytes (0.3 MB / 0 MWords) used by master data server.

 ----------------------------------------
 CPU timing information for all processes
 ========================================
 0: 1.992 + 0.24 = 2.16
 ----------------------------------------
 ddikick.x: exited gracefully.
unset echo
----- accounting info -----
Files used on the master node doreen were:
-rw-r--r-- 1 lmentel lmentel   26187 maj 11 14:57 /home/lmentel/scratch/he-h2_avdz_ormas.dat
-rw-r--r-- 1 lmentel lmentel    1450 maj 11 14:57 /home/lmentel/scratch/he-h2_avdz_ormas.F05
-rw-r--r-- 1 lmentel lmentel  360032 maj 11 14:57 /home/lmentel/scratch/he-h2_avdz_ormas.F08
-rw-r--r-- 1 lmentel lmentel  363064 maj 11 14:57 /home/lmentel/scratch/he-h2_avdz_ormas.F09
-rw-r--r-- 1 lmentel lmentel 1799600 maj 11 14:57 /home/lmentel/scratch/he-h2_avdz_ormas.F10
-rw-r--r-- 1 lmentel lmentel  266376 maj 11 14:57 /home/lmentel/scratch/he-h2_avdz_ormas.F12
ls: No match.
ls: No match.
ls: No match.
nie, 11 maj 2014, 14:57:35 CEST
0.0u 0.0s 0:06.23 0.4% 0+0k 816+5488io 4pf+0w'''

    def setUp(self):
        temp = tempfile.NamedTemporaryFile(delete=False)
        temp.write(TestGLPonHeH2.heh2_log)
        self.glp = GamessLogParser(temp.name)

    def test_logexists(self):
        self.assertTrue(self.glp.logexists())

    def test_accomplished(self):
        self.assertTrue(self.glp.accomplished())

    def test_get_version(self):
        self.assertEqual(self.glp.get_version(), "1 MAY 2013 (R1)")

    def test_get_charge(self):
        self.assertEqual(self.glp.get_charge(), 0)

    def test_get_electrons(self):
        self.assertEqual(self.glp.get_electrons(), 4)

    def test_get_homo(self):
        self.assertEqual(self.glp.get_homo(), 1)

    def test_get_number_of_atoms(self):
        self.assertEqual(self.glp.get_number_of_atoms(), 3)

    def test_get_number_of_aos(self):
        self.assertEqual(self.glp.get_number_of_aos(), 27)

    def test_get_number_of_mos(self):
        self.assertEqual(self.glp.get_number_of_mos(), 27)

    def test_get_linear_deps(self):
        self.assertIs(self.glp.get_linear_deps(), None)

    def test_get_scf_type(self):
        self.assertEqual(self.glp.get_scf_type(), "RHF")

    def test_get_cc_type(self):
        self.assertEqual(self.glp.get_cc_type(), "NONE")

    def test_get_ci_type(self):
        self.assertEqual(self.glp.get_ci_type(), "ORMAS")

    def test_get_mplevel(self):
        self.assertEqual(self.glp.get_mplevel(), 0)

    def test_get_hf_total_energy(self):
        self.assertAlmostEqual(self.glp.get_hf_total_energy(), -3.9842266066)

    def test_get_ormas_total_energy(self):
        self.assertAlmostEqual(self.glp.get_ormas_total_energy(), -4.0544723546)

    def test_get_energy_components_hf(self):
        d = {"WAVEFUNCTION NORMALIZATION" : 1.0000000000,
             "ONE ELECTRON ENERGY" : -7.6711052217,
             "TWO ELECTRON ENERGY" : 2.3347610066,
             "NUCLEAR REPULSION ENERGY" : 1.3521176085,
             "TOTAL ENERGY" : -3.9842266066,
             "ELECTRON-ELECTRON POTENTIAL ENERGY" : 2.3347610066,
             "NUCLEUS-ELECTRON POTENTIAL ENERGY" : -11.5902333599,
             "NUCLEUS-NUCLEUS POTENTIAL ENERGY" : 1.3521176085,
             "TOTAL POTENTIAL ENERGY" : -7.9033547447,
             "TOTAL KINETIC ENERGY" :3.9191281381,
             "VIRIAL RATIO (V/T)" : 2.0166104465}

        self.assertDictEqual(d, self.glp.get_energy_components('hf'))

    def test_get_energy_components_ormas(self):
        d = {"WAVEFUNCTION NORMALIZATION" : 1.0000000000,
             "ONE ELECTRON ENERGY" : -7.6031513994,
              "TWO ELECTRON ENERGY" : 2.1965614362,
              "NUCLEAR REPULSION ENERGY" : 1.3521176085,
              "TOTAL ENERGY" : -4.0544723546,
              "ELECTRON-ELECTRON POTENTIAL ENERGY" : 2.1965614362,
              "NUCLEUS-ELECTRON POTENTIAL ENERGY" : -11.5724441222,
              "NUCLEUS-NUCLEUS POTENTIAL ENERGY" : 1.3521176085,
              "TOTAL POTENTIAL ENERGY" : -8.0237650775,
              "TOTAL KINETIC ENERGY" : 3.9692927229,
              "VIRIAL RATIO (V/T)" : 2.0214596498}

        self.assertDictEqual(d, self.glp.get_energy_components('ci'))

if __name__ == "__main__":
    unittest.main()
