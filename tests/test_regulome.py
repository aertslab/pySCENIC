# -*- coding: utf-8 -*-

from pyscenic.algo import module2regulome, module2features_auc1st_impl, modules2df

import os, yaml, glob
from functools import partial
from configparser import ConfigParser
from pyscenic.rnkdb import FeatherRankingDatabase as RankingDatabase, SQLiteRankingDatabase
from pyscenic.genesig import GeneSignature, Regulome
from pyscenic.utils import load_motif_annotations


TEST_DATABASE = "hg19-500bp-upstream-10species"
TEST_SIGNATURE = "msigdb_cancer_c6"
RESOURCES_FOLDER="/Users/bramvandesande/Projects/lcb/resources"
MOTIF_ANNOTATIONS_FNAME = os.path.join(RESOURCES_FOLDER, "motifs-v9-nr.mgi-m0.001-o0.0.tbl")
DATA_FOLDER = "/Users/bramvandesande/Projects/lcb/tmp"
DATABASE_FOLDER = "/Users/bramvandesande/Projects/lcb/databases/"
FEATHER_GLOB = os.path.join(DATABASE_FOLDER, "mm9-*.feather")
MODULE1 = Regulome(name='Regulome for Bhlhe41', nomenclature='MGI',
                   gene2weights={'Ddx17': 0, 'Hnrnpa2b1': 1, 'Ubap2l': 2, 'Dagla': 3, '4933434E20Rik': 5,
                                 'Ythdc1': 1092, 'Dctn4': 8, 'Ksr2': 10, 'Arih1': 16, 'Braf': 18, 'Ak4': 19,
                                 'Hnrnph2': 21, 'Josd1': 25, 'Btbd9': 27, 'Brd4': 28, 'Gpr137': 30, 'Gpbp1': 33,
                                 'Picalm': 37, 'Slc25a3': 41, 'Lmbr1': 44, 'Egr4': 1031, 'BC048403': 1040, 'Nrg2': 962,
                                 'Sar1a': 1168, 'Cdk2ap2': 57, 'Ube2q1': 61, 'Cltc': 727, 'Ccdc28a': 501, 'Sez6l2': 69,
                                 'Ell2': 71, 'Impdh2': 72, 'Fam50a': 74, 'Bad': 75, 'Hnrnpa3': 76, 'Ubr4': 633,
                                 'Ube2h': 80, 'Gigyf2': 81, 'Snap25': 83, 'Gas5': 289, 'Ldha': 85, 'Tmem59l': 87,
                                 'Hspa9': 89, 'Magoh': 1198, 'Ppp2ca': 99, 'Dixdc1': 101, 'Mrps18b': 103, 'Uso1': 105,
                                 'Hnrnpd': 161, 'Syncrip': 107, 'Snx2': 108, 'Ddx3x': 116, 'Mesdc1': 117, 'Klhl12': 120,
                                 'Pgrmc2': 121, 'Prrc2c': 122, 'Insig1': 126, 'Ubr5': 128, '3110043O21Rik': 129,
                                 '1700066M21Rik': 132, 'Celf6': 135, 'Nmnat2': 138, 'Igf2bp3': 140, 'Cxxc1': 141,
                                 'Zdbf2': 142, '9530068E07Rik': 144, 'Srsf7': 145, 'Bdp1': 146, 'Gm13363': 149,
                                 'Spred2': 153, 'Tra2a': 157, 'Cdc42se1': 161, 'Clk4': 163, 'Cds1': 164, 'Smarca5': 165,
                                 'Hmgn1': 459, 'Gak': 171, 'Usp48': 172, 'Cnnm4': 173, 'Foxred1': 215, 'Spag9': 177,
                                 'Prkce': 181, 'Creld1': 1083, 'Mdga2': 189, 'Hnrnph1': 190, 'Sgol1': 192, 'Rps9': 195,
                                 'Ilf3': 196, 'Ppig': 377, 'Cnot7': 211, 'Calm2': 213, 'Nrd1': 214, 'Hnrnpu': 215,
                                 'Hspa4': 755, 'Hmg20a': 221, 'Casc3': 224, 'Emd': 228, 'Eif4b': 231, 'Npr3': 235,
                                 'Rps8': 243, 'BC031181': 949, 'Mllt11': 246, 'Gzf1': 249, 'Egr1': 258, 'Stard3nl': 256,
                                 'Srrm2': 258, 'Wdr1': 403, 'Fbxo33': 809, 'Pex10': 272, 'Ghitm': 273, 'Dennd5a': 275,
                                 'Mfsd1': 278, 'Dpm1': 281, 'Appbp2': 281, 'Pitpna': 286, 'Dnajb2': 297, 'Map2k1': 298,
                                 'Enox1': 1016, 'Tmem198': 309, 'Chpf': 310, 'Rad23b': 313, 'Abhd16a': 314,
                                 'Pafah1b1': 366, 'Nrbf2': 320, 'Arhgap12': 323, 'Hat1': 578, 'Mir132': 327, 'Hlf': 328,
                                 'Trpm7': 330, 'Ap2b1': 334, 'Zfand2b': 336, 'Arpc5': 338, 'Snapc1': 339,
                                 'Tmem165': 341, 'Gnao1': 346, 'Pgk1': 348, 'Napa': 355, '2510039O18Rik': 356,
                                 'Commd3': 357, 'Ppp1cb': 361, 'Tmem147': 366, 'Adra1b': 374, 'Ift81': 844,
                                 'Lrrc14': 379, 'Ormdl3': 385, 'Slitrk1': 386, 'Spock2': 391, 'Hcn1': 392, 'Mib2': 393,
                                 'Trp53inp2': 406, 'Clasp2': 407, 'Atp6ap2': 410, 'Mmgt1': 417, 'Efr3a': 419,
                                 'Etv3': 420, 'Hook2': 428, 'Tnrc6b': 429, 'Mbtd1': 440, 'Phf20l1': 940, 'Utp18': 443,
                                 'Dach1': 447, 'Cpsf6': 448, 'Piga': 451, 'Ntn1': 453, 'Celf1': 467, 'Sppl3': 469,
                                 'Vamp3': 475, 'Usp3': 894, 'Ubn1': 481, 'Rprd1a': 485, 'Glyr1': 489, 'G3bp2': 1074,
                                 'Lrrtm1': 496, 'Rad51': 500, 'Ndel1': 502, 'Arhgap21': 504, 'Dcaf11': 505,
                                 'Ywhag': 506, 'Zfp523': 510, 'Gtf2a1l': 521, 'Armcx2': 522, 'Hnrnpa0': 837,
                                 'Helq': 1169, 'Med14': 525, 'Atl2': 527, 'Pls3': 529, 'Dhx36': 533, 'Rrp1': 534,
                                 'Slitrk4': 535, 'Dnaja3': 542, 'Utp23': 547, 'Uqcrfs1': 549, 'Tmem183a': 551,
                                 'Mapk8ip3': 552, 'Ccna2': 555, 'Pdap1': 559, 'Drap1': 562, 'Nrxn1': 563, 'Otub2': 565,
                                 'Tubb2a': 566, 'Vps33a': 568, 'Rbks': 571, 'Gabarap': 574, 'Bre': 577, 'Gtf3c1': 578,
                                 'Ubac1': 580, 'Rbm26': 586, 'Golga5': 591, 'Rab3a': 595, 'Prpf38a': 598, 'Hdac9': 599,
                                 'Dnajb11': 600, 'Golph3': 604, 'Abl2': 608, 'Irgq': 611, 'Sirpa': 613, 'Atf1': 614,
                                 'Eif3k': 912, 'Galnt11': 620, 'Lztfl1': 798, 'Ube2s': 624, 'Ythdc2': 630,
                                 '2610017I09Rik': 634, 'Cdc5l': 636, 'Cxcl16': 639, 'Ttc19': 641, 'Mex3c': 648,
                                 'Csmd3': 649, 'Usp2': 651, 'Nfix': 652, 'Gpcpd1': 653, 'Pqbp1': 654, 'Kdm4b': 657,
                                 'Ddhd1': 659, 'Clpx': 660, 'Ccnl1': 662, 'Kat2b': 664, 'Cacybp': 668, 'Morf4l2': 1023,
                                 'Baz1b': 670, 'Fam160b1': 674, 'Prkar1a': 1057, 'Zranb3': 683, 'Ube2o': 684,
                                 'Utp15': 686, 'Kdm6a': 687, 'Xpo7': 688, 'Bcas2': 702, 'Otud7a': 706, 'Ptpru': 709,
                                 'Cpe': 713, 'Htatsf1': 714, 'Stx6': 718, 'Phlpp1': 720, 'Adam9': 723, 'Gpm6b': 724,
                                 'Orc4': 1038, 'Klf16': 727, 'Tm2d2': 729, 'Safb': 731, 'Nr4a1': 733, 'Tspan31': 734,
                                 'Nbeal1': 735, 'Map3k7': 737, 'Pkp4': 739, 'Mbd5': 1044, 'Cisd2': 744, 'Pdik1l': 745,
                                 'Slco5a1': 750, 'Dusp6': 753, 'Zfp628': 754, 'Ankrd17': 757, 'Nudt9': 761,
                                 'Ssr4': 1160, 'Gm9855': 765, 'Sumo3': 768, 'Ifrd1': 770, 'Psenen': 773, 'Atp1b3': 774,
                                 'St13': 777, 'Pwwp2b': 778, 'Dap3': 785, 'U2af1l4': 787, 'Pcif1': 795, 'Fam60a': 796,
                                 'Rab11fip2': 797, 'Klf13': 799, 'Ppm1a': 800, 'Ndufb6': 803, 'Senp2': 911,
                                 'Necap1': 815, 'Ube2d3': 819, 'Fam174a': 820, 'Dync1h1': 824, 'Bbc3': 827,
                                 'Sstr3': 828, 'Cpeb2': 829, 'Srrm4': 834, 'Psmd3': 836, 'Fgf9': 837, 'Tmf1': 941,
                                 'Trak1': 846, 'Ttc1': 847, 'Elavl2': 853, 'Atp5g1': 854, 'Mapk14': 859, 'Uqcrc2': 861,
                                 'Mid1ip1': 863, 'Narg2': 876, 'Snw1': 878, 'Atg5': 881, 'Cry2': 885, 'Camk2n2': 886,
                                 'Dusp8': 910, 'Lphn1': 913, 'Hexim2': 917, 'Tbc1d23': 918, 'Ociad1': 924, 'Cbx6': 925,
                                 'Adnp2': 926, 'Tob1': 932, 'Dnttip1': 938, 'Mn1': 939, 'Enkur': 940, 'Yme1l1': 941,
                                 'Baz2b': 942, 'Smarcad1': 946, 'Tcf4': 952, 'Rbm12': 955, 'Ddx42': 959, 'Slc35f3': 965,
                                 'Sbds': 966, 'Ube2r2': 969, 'Fam169a': 975, 'Cacna1a': 979, 'Pef1': 990, 'Gpkow': 994,
                                 'Epha4': 997, 'Zfp36l2': 1001, 'Bicd2': 1006, 'Uspl1': 1016, 'Fam76a': 1022,
                                 'Epn1': 1024, 'Fam135a': 1025, 'Clpb': 1027, 'Pdzd8': 1028, 'Cltb': 1032, 'Nono': 1033,
                                 'Fam179b': 1043, 'Sidt2': 1044, 'Sec23b': 1046, 'Fem1a': 1053, 'Ino80e': 1054,
                                 'Ppp2r5a': 1060, 'Kcnab2': 1065, 'Per2': 1067, 'Rbx1': 1068, '1190002N15Rik': 1069,
                                 'Atp6v0d1': 1071, 'Poc1b': 1072, 'Tuba1b': 1075, 'Coq2': 1079, 'Lmo4': 1080,
                                 'Mapk8': 1083, 'Hdac3': 1086, 'Pip5k1a': 1088, 'Tial1': 1089, 'Nfx1': 1090,
                                 'Usp28': 1095, 'Slc25a28': 1097, 'Sik3': 1099, 'Oaz1': 1100, 'Lias': 1101,
                                 'Ahcyl1': 1104, 'Ippk': 6, 'Cog3': 3, 'Lypla2': 20, 'Vps52': 15, 'Ppp4c': 17,
                                 'Gabpa': 63, 'Csnk1g1': 83, 'Trip4': 84, 'Gga1': 35, 'Psmb7': 93, 'Usf1': 44,
                                 'Bag4': 80, 'Ppp1r16a': 66, 'Gtpbp4': 287, 'Mtmr4': 96, 'Insig2': 52, 'Gabarapl2': 208,
                                 'Ccz1': 845, 'Sec31a': 722, 'Zfp384': 188, 'Eif3i': 347, 'Parg': 68, 'Timm23': 69,
                                 'Zranb2': 213, 'Psmd6': 114, 'Edc4': 225, 'Atp6v1e1': 519, 'Kctd5': 94, 'Atp2c1': 314,
                                 'Shoc2': 231, 'Cops3': 197, 'Med6': 268, 'Serp1': 123, 'Dnajc14': 745, 'Urm1': 98,
                                 'Eif3g': 99, 'Srsf5': 149, 'Eif3l': 253, 'Cap1': 131, 'Hccs': 284, 'Edem3': 130,
                                 'Nek9': 136, 'Cbll1': 128, '2700062C07Rik': 129, 'Ipo4': 338, 'Arpc3': 140,
                                 'Snrpd2': 491, 'Dpp8': 207, 'Arpc2': 145, 'Bola2': 360, 'Hcfc1': 296, 'Ndufb5': 633,
                                 'Uba2': 332, 'Trim35': 183, 'Sec11a': 657, '0610010K14Rik': 440, 'Mta2': 189,
                                 'Cnot10': 697, 'Ppid': 191, 'Hspbp1': 193, 'Mlec': 316, 'Cdc123': 270, 'Rsl24d1': 262,
                                 'Atxn2': 248, 'Tada3': 810, 'Chd2': 389, 'Klhl20': 216, 'Pik3c3': 361, 'Snf8': 436,
                                 'Rnf167': 222, 'Crk': 241, 'Incenp': 537, 'Senp7': 232, 'Acin1': 320, 'Thpo': 790,
                                 'Mark3': 502, 'Mrpl43': 403, 'Sh3bp5l': 238, 'Skiv2l': 239, 'Peo1': 405, 'Rac1': 244,
                                 'Arfip1': 245, 'Aftph': 250, 'Snrpb': 443, 'Snrnp27': 268, 'Vps25': 269, 'Dmap1': 272,
                                 'Frs2': 687, 'Fbxo28': 277, 'Nfat5': 1200, 'Asb3': 659, 'Yipf4': 286, 'Aptx': 830,
                                 'Rpl10': 743, 'Tm9sf3': 498, 'Dazap2': 497, 'Pex6': 397, 'Cttnbp2nl': 1015,
                                 'Mrpl10': 890, 'Ap3s2': 314, 'Tomm70a': 321, 'Suds3': 867, 'Immt': 324, 'Utp3': 331,
                                 'Wwp2': 332, 'Ddx39b': 525, 'Ywhab': 339, 'Pde12': 751, 'Sec24c': 640, 'Srsf2': 526,
                                 'Exoc1': 669, 'Kctd10': 427, 'Qars': 523, 'Clasp1': 364, 'Rpl37a': 1048, 'Arfip2': 814,
                                 'Tnks': 379, 'Vps53': 846, 'Mrps18a': 724, 'Fam188a': 380, 'Txnl4b': 381,
                                 'Thap11': 1002, 'Dhx38': 383, 'Ebi3': 1046, 'Itpa': 389, 'Mrpl3': 395, 'l7Rn6': 399,
                                 'Prpf19': 400, 'Uba5': 547, 'Sec61a1': 411, 'Stam2': 1141, 'Flot2': 796, 'Setd8': 423,
                                 '2810403A07Rik': 593, 'Pcnxl3': 428, 'Gba2': 429, 'Eif3m': 431, 'Zfp3': 805,
                                 'A230072C01Rik': 450, 'Rps17': 1004, 'Atm': 582, 'Npat': 585, 'Zfc3h1': 459,
                                 'Gdap2': 823, 'Trim44': 462, 'Ube2m': 470, 'Klhdc3': 622, 'Uggt1': 749, 'Snrpe': 479,
                                 'Cops7a': 483, 'Erh': 485, 'Snip1': 1174, 'Apbb3': 1017, 'Tprgl': 801, 'Det1': 496,
                                 'Cggbp1': 500, 'Prkacb': 656, 'Nars': 894, 'Pld3': 701, '2310022A10Rik': 703,
                                 'Acp2': 762, 'Derl1': 526, 'Uba6': 530, 'Srp72': 532, 'Dcun1d5': 988, 'Setdb1': 559,
                                 'Rnf185': 549, '8430429K09Rik': 552, 'Mrps16': 553, 'Atf6b': 554, 'Eif5a2': 865,
                                 'Ralgapb': 558, 'Aatf': 938, 'Carkd': 567, 'Vta1': 682, 'Ppp2r5e': 573, 'H2-Ke2': 574,
                                 'Wdr48': 794, 'Psmb1': 583, 'Sf3b2': 584, 'Banf1': 821, 'Lemd3': 759, 'Zfp668': 597,
                                 'Trappc4': 603, 'Srp19': 839, 'Idh3a': 613, 'Rars': 813, 'Grcc10': 618, 'Rpl28': 628,
                                 'Stk4': 637, 'Polr2d': 1153, 'Myeov2': 661, 'Ddx24': 650, 'Actr2': 655, 'Rpl19': 658,
                                 'Gm12191': 660, 'Ccdc104': 1065, 'Uba52': 667, 'Ercc3': 669, '1110059E24Rik': 670,
                                 'Cops5': 671, 'Psen1': 678, 'Fbxo32': 679, 'Ubxn6': 685, 'BC003965': 688, 'Ap1ar': 834,
                                 'Atp5h': 696, 'Kctd1': 699, 'Bcl10': 704, 'Fkbp1a': 709, 'Ubr1': 718, 'Smchd1': 719,
                                 'Zkscan4': 1104, 'Bcl2l13': 722, 'Gtl3': 800, 'Actr3': 728, 'Arpc4': 959, 'Rpl37': 734,
                                 'D15Ertd621e': 735, 'Surf2': 736, 'Oraov1': 980, 'Rad23a': 739, 'Fam63b': 752,
                                 'Ndufs7': 919, 'Ttc9': 750, 'Xpc': 751, 'Nuf2': 1178, 'Naa20': 1163, 'Mcm8': 764,
                                 'Tusc2': 766, 'Tfg': 767, 'Cstf2': 770, 'Ugdh': 779, 'Tlk1': 836, 'Slc25a46': 783,
                                 'Zswim3': 785, 'Golga3': 789, 'Ogt': 1107, 'Trim41': 795, 'Snd1': 797, 'B3gat3': 800,
                                 'Psmc6': 802, 'Coro1c': 818, 'Usp24': 806, 'Mga': 811, 'Dzip3': 812, 'Zfp598': 819,
                                 'Ssh2': 1043, 'Timm8b': 826, 'Zkscan6': 1186, 'Eif2s3x': 838, 'Mlycd': 839,
                                 'Rbck1': 845, 'Nfs1': 849, 'Hsph1': 941, 'Poldip3': 851, 'Ovca2': 859, 'Dmtf1': 1156,
                                 'Mapk1ip1l': 866, 'Midn': 1165, 'Glrx5': 868, 'Mogs': 872, 'Ebag9': 878, 'Tcp1': 880,
                                 'Trappc10': 885, 'Txndc9': 887, 'Ahsa1': 890, 'Pcbp2': 892, 'Fnta': 895, 'Sike1': 898,
                                 'Uqcrq': 1053, 'Gpatch2': 900, 'Ndufs2': 1176, 'Zfp60': 911, 'Zkscan1': 917,
                                 'Bag1': 1129, 'Ado': 977, '1110037F02Rik': 933, 'Chmp5': 1138, 'Acvr1': 1111,
                                 'Usp9x': 943, 'U05342': 1185, 'Nup54': 964, 'Laptm5': 949, 'Abi2': 952, 'Pak3': 953,
                                 'Rpl7l1': 954, 'Snx19': 955, 'Zfp780b': 956, 'Ehbp1l1': 957, 'Hbp1': 958,
                                 'Slc25a20': 960, 'Eif4a1': 1173, 'Arhgap1': 964, 'Ltv1': 969, 'Fbxl6': 971,
                                 'Cdk7': 975, 'Zswim4': 976, '1700008J07Rik': 977, 'Wac': 984, 'Zc3h15': 985,
                                 'Atp5c1': 986, 'Cfdp1': 992, '5430417L22Rik': 993, 'Trip12': 998, 'Sp3': 1175,
                                 'Cdc27': 1000, 'Ankrd13d': 1001, 'Rab3gap2': 1004, 'Nop58': 1005, 'Nedd8': 1014,
                                 'Pcbp1': 1019, 'Fubp3': 1022, 'C030016D13Rik': 1024, 'D19Bwg1357e': 1025, 'Ing1': 1030,
                                 'Kdm5a': 1035, 'Ccdc6': 1037, 'Gosr2': 1085, 'Lztr1': 1041, 'Clta': 1043,
                                 'Dcaf12': 1045, 'Junb': 1048, 'Stk19': 1053, 'Usp11': 1059, 'Cnst': 1060,
                                 'Ruvbl1': 1065, 'Snx3': 1066, 'Fam65b': 1071, 'Tnfsf13': 1072, 'Trmt5': 1075,
                                 'Rbm19': 1076, 'Cc2d1a': 1079, 'Ccng2': 1084, 'Sbno1': 1110, 'Krit1': 1090,
                                 'Slc25a14': 1091, 'Faf2': 1092, 'Ppm1g': 1097, 'Rps27a': 1111, 'Hnrnpr': 1115,
                                 'Foxj3': 1117, 'Ptpn1': 1118, 'Ccnc': 1123, 'Purb': 1124, 'Lrrc28': 1126,
                                 'Gm4532': 1127, 'Smc1a': 1143, 'Arf3': 1147, 'Dpy19l1': 1149, 'Cog1': 1150,
                                 'Spred1': 1153, 'Rwdd4a': 1156, 'U2af2': 1164, 'Gpr155': 1166, 'Myg1': 1170,
                                 'Prkag2': 1178, 'Ubp1': 1182, 'Otud7b': 1194, 'Golga7': 1195, 'Naa38': 1199,
                                 'Pank3': 1200, 'Btbd1': 1205, 'Clns1a': 1208, 'D2Wsu81e': 1212, 'Aph1a': 1213,
                                 'Asxl1': 1215, 'Dock4': 1216, 'Zfp277': 1221, 'Exosc6': 1223, 'Ubxn1': 1226,
                                 'Vti1a': 1229, 'Pms2': 1231, '2310057M21Rik': 1237, 'Rnf14': 1239, 'Atp5sl': 1241,
                                 'Gm561': 1242, 'Nosip': 58, 'Mfap1b': 156, 'Sec22b': 196, 'Slc39a6': 201, 'Smu1': 303,
                                 'Zbtb11': 307, 'Ccdc39': 377, 'Gltpd1': 420, 'Trio': 425, 'Wdpcp': 430, 'Mdh1': 432,
                                 'Gpn2': 467, 'Ints1': 501, 'Ddx19b': 507, 'Lrrc41': 508, 'Sub1': 510, 'Zbed6': 512,
                                 'Kdelc1': 515, 'Chtf8': 544, 'Cirh1a': 546, 'Rrp12': 549, 'Tars2': 573, 'Uchl4': 596,
                                 'Skil': 602, 'Slc25a11': 604, 'Dpy19l3': 626, 'Ubr7': 632, 'Ap2a1': 638,
                                 'Tmem109': 646, 'Tmx2': 673, 'Emg1': 694, 'Gprin1': 708, 'Mrpl27': 710, 'Phrf1': 712,
                                 'Ate1': 719, 'Arhgef3': 727, 'Zfp27': 736, 'Pik3r2': 740, 'Eif2s1': 744,
                                 'E130308A19Rik': 785, 'Cd68': 822, 'Ybx1': 840, 'Sgta': 857, 'Smn1': 858,
                                 'Cstf2t': 860, 'Thumpd3': 863, 'Kcnip4': 870, 'Rpl4': 884, 'Pars2': 908, 'Rnf214': 910,
                                 'Ipo13': 926, 'Rbm34': 927, 'Tab3': 940, 'Cdc26': 942, 'Mterfd2': 950, 'Arhgef7': 952,
                                 'Rfc3': 955, 'Bms1': 970, 'Ap2s1': 979, 'Ufc1': 981, 'Rnf181': 1024, 'Serpinb6a': 1033,
                                 'Cdk5rap1': 1056, 'Naa50': 1062, 'Ubl4': 1067, 'Hoxa7': 1072, 'Eya3': 1081,
                                 'Nfkbiz': 1086, '0610030E20Rik': 1091, 'Dync1i2': 1109, 'Rnf10': 1117, 'Etfdh': 1126,
                                 'Ldb2': 1143, 'Mettl3': 1168, 'Cyld': 1170, 'Smad5': 1183, 'Fbrsl1': 1184,
                                 'Ube4b': 1191, 'Coq6': 1194, 'Nt5c': 1196}, transcription_factor='Bhlhe41',
                   context=frozenset({'target weight >= 0.001'}),
                   score=0.17450369340686275)
MODULE2 = Regulome(name='Regulome for Ahr',
                   nomenclature='MGI',
                   gene2weights={'Vps33a': 1.1945297414493656, 'Picalm': 0.06954951190775277,
                                 'Arhgap12': 0.2191417967550388,
                                 'Eif4g1': 0.6522033708234853, '4933434E20Rik': 0.032645005366486374,
                                 'Cuta': 0.11830955263689501, 'Dolk': 0.08819765411935789,
                                 'Gtf2h1': 0.09891245083048687,
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
                                 'Josd1': 0.3633953929958216, 'Dusp5': 0.3860816513401068,
                                 'Cbfa2t3': 0.2655666147921891,
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
                   context=frozenset({'target weight >= 0.001'}),
                   score=3.3525298840246283)


def load_database():
    db_fnames = glob.glob(FEATHER_GLOB)

    def name(fname):
        return os.path.basename(fname).split(".")[0]

    return RankingDatabase(fname=db_fnames[1], name=name(db_fnames[1]), nomenclature="MGI")


DATABASE = load_database()
MOTIF_ANNOTATIONS = load_motif_annotations(MOTIF_ANNOTATIONS_FNAME)


def load_db_info(section):
    config = ConfigParser()
    config.read(os.path.join(os.path.dirname(__file__), 'test_sqlitedb.ini'))
    return config[section]

def load_gs_info(section):
    config = ConfigParser()
    config.read(os.path.join(os.path.dirname(__file__), 'test_genesig.ini'))
    return config[section]

def test_module2regulome():
    gs = GeneSignature.from_gmt(gene_separator="\t", field_separator="\t", **load_gs_info(TEST_SIGNATURE))[3]
    db = SQLiteRankingDatabase(**load_db_info(TEST_DATABASE))
    module = Regulome(gs.name, gs.nomenclature, gs.gene2weights, "TP53")
    motif_annotations = load_motif_annotations(MOTIF_ANNOTATIONS_FNAME)
    reg = module2regulome(db, module, motif_annotations)

def test_modules2df():
    module2features = partial(module2features_auc1st_impl,
                              rank_threshold = 1500, auc_threshold = 0.05, nes_threshold=2.0,
                              avgrcc_sample_frac = None)
    df = modules2df(DATABASE, [MODULE1, MODULE2], MOTIF_ANNOTATIONS, module2features_func=module2features,
                    return_recovery_curves=True)
    print(df)

def test_to_from_yaml():
    gs = GeneSignature.from_gmt(gene_separator="\t", field_separator="\t", **load_gs_info(TEST_SIGNATURE))[3]
    regulome = Regulome(gs.name, gs.nomenclature, gs.gene2weights, "TP53")
    yml = yaml.dump(regulome)
    regulome2 = yaml.load(yml)
    assert regulome.name == regulome2.name
    assert regulome.nomenclature == regulome2.nomenclature
    assert regulome.genes == regulome2.genes
    assert regulome.weights == regulome2.weights
    assert regulome.score == regulome2.score
    assert regulome.context == regulome2.context

