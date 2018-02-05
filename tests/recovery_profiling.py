# -*- coding: utf-8 -*-

import os
import glob
from pyscenic.rnkdb import FeatherRankingDatabase as RankingDatabase
from pyscenic.utils import load_motif_annotations
from pyscenic.regulome import module2regulome
from pyscenic.genesig import Regulome

DATA_FOLDER = "/Users/bramvandesande/Projects/lcb/tmp"
DATABASE_FOLDER = "/Users/bramvandesande/Projects/lcb/databases/"
FEATHER_GLOB = os.path.join(DATABASE_FOLDER, "mm9-*.feather")
RESOURCES_FOLDER = "/Users/bramvandesande/Projects/lcb/resources"
MOTIF_ANNOTATIONS_FNAME = os.path.join(RESOURCES_FOLDER, "motifs-v9-nr.mgi-m0.001-o0.0.tbl")
MODULE = Regulome(name='Regulome for Ahr',
                  nomenclature='MGI',
                  gene2weights={'Vps33a': 1.1945297414493656, 'Picalm': 0.06954951190775277,
                                'Arhgap12': 0.2191417967550388,
                                'Eif4g1': 0.6522033708234853, '4933434E20Rik': 0.032645005366486374,
                                'Cuta': 0.11830955263689501, 'Dolk': 0.08819765411935789, 'Gtf2h1': 0.09891245083048687,
                                'Srek1': 0.2616072413468973, 'Mthfr': 3.3079611119564483, 'Aff4': 0.3852401021008935,
                                'Rragb': 0.9350505880315496, 'Zzz3': 0.02439228113312613, 'Spag9': 0.12435951278081875,
                                'Ube2q1': 0.915679870965199, 'Nr1d1': 0.2375960570212787, 'Snapc3': 0.0851670128377308,
                                '1700052N19Rik': 0.08376741529946874, 'Cpeb1': 0.09812299270270064,
                                'Gnptg': 0.1334133038249925, 'Usp8': 0.13283546755415526, 'Mfsd1': 0.14023199396320168,
                                'Ubn1': 0.8616985795254117, 'Luc7l2': 0.04299307530991448,
                                'Mcoln1': 0.31392862792431875, 'Phf15': 4.971526455181718, 'Pcdh8': 0.2669756810193168,
                                'Lox': 0.7017261724791286, 'Socs5': 1.3244789367947636, 'Spg21': 0.1050236179972574,
                                'Atp2b1': 0.028907699902541082, 'Erp29': 0.6762035182679148,
                                'Srp54c': 0.3555656329029193, 'Atp11b': 0.0038433337203616826,
                                'Tug1': 3.2026251680547984, 'Abhd11': 0.7446607065618946, 'Aldh3b2': 0.115266748300316,
                                'Mal2': 0.01446846879743427, 'Flrt2': 0.1888844895785373, 'Olfm1': 0.18176116505215054,
                                'Prickle2': 0.034074011686413545, 'Tmem110': 0.1139195676268978,
                                'Psma6': 0.6373087425114722, 'Nptn': 0.02942845096830972, 'U2af2': 0.12370752319947155,
                                'Josd1': 0.3633953929958216, 'Dusp5': 0.3860816513401068, 'Cbfa2t3': 0.2655666147921891,
                                'Ptprs': 0.9282687320971084, 'Trim30d': 0.04656023384483033,
                                'Mtdh': 0.08020983833663681, 'Ralyl': 0.13726478430829145,
                                'Slc36a1': 0.03753240136558692, 'Fbxo10': 0.5558938164503024,
                                'Sla2': 0.04533692716224516, 'Zmat3': 0.0513161638419828,
                                'Morf4l1': 0.055595088722155774, 'Bdnf': 0.06583990152915586,
                                'Jph1': 0.12892438043957755, 'Lin7c': 0.11488447713265965, 'Uvrag': 1.6689135521176337,
                                'Irf2bp2': 17.522352407417436, 'Cbx6': 0.5170388530245495, 'Grm8': 0.1678087627689279,
                                'Npcd': 0.14792667790484118, 'Zfp9': 0.09381393299800976, 'Rfxap': 0.17224099000280213,
                                'Khdrbs1': 0.08601662566993558, 'Esrrg': 0.02169651731153216,
                                'Kat5': 0.7542129250430151, 'Farsa': 0.02671753000385724, 'Foxp1': 0.5111233815033093,
                                '3110002H16Rik': 0.28325916910551197, 'Slc25a31': 0.4155416912611951,
                                'Get4': 0.1246366302020446, 'Atxn3': 4.154750518729518},
                  transcription_factor='Ahr',
                  context=frozenset({'mm9-tss-centered-10kb-10species', 'target weight >= 0.001'}),
                  score=3.3525298840246283)


def load_database():
    db_fnames = glob.glob(FEATHER_GLOB)

    def name(fname):
        return os.path.basename(fname).split(".")[0]

    return RankingDatabase(fname=db_fnames[0], name=name(db_fnames[0]), nomenclature="MGI")


DATABASE = load_database()
MOTIF_ANNOTATIONS = load_motif_annotations(MOTIF_ANNOTATIONS_FNAME)


def recovery():
    module2regulome(DATABASE, MODULE, MOTIF_ANNOTATIONS)


if __name__ == "__main__":
    recovery()
