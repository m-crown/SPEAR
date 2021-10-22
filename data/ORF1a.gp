LOCUS       YP_009725295            4405 aa            linear   VRL 18-JUL-2020
DEFINITION  ORF1a polyprotein [Severe acute respiratory syndrome coronavirus
            2].
ACCESSION   YP_009725295
VERSION     YP_009725295.1
DBLINK      BioProject: PRJNA485481
DBSOURCE    REFSEQ: accession NC_045512.2
KEYWORDS    RefSeq.
SOURCE      Severe acute respiratory syndrome coronavirus 2 (SARS-CoV-2)
  ORGANISM  Severe acute respiratory syndrome coronavirus 2
            Viruses; Riboviria; Orthornavirae; Pisuviricota; Pisoniviricetes;
            Nidovirales; Cornidovirineae; Coronaviridae; Orthocoronavirinae;
            Betacoronavirus; Sarbecovirus.
REFERENCE   1  (residues 1 to 4405)
  AUTHORS   Wu,F., Zhao,S., Yu,B., Chen,Y.M., Wang,W., Song,Z.G., Hu,Y.,
            Tao,Z.W., Tian,J.H., Pei,Y.Y., Yuan,M.L., Zhang,Y.L., Dai,F.H.,
            Liu,Y., Wang,Q.M., Zheng,J.J., Xu,L., Holmes,E.C. and Zhang,Y.Z.
  TITLE     A new coronavirus associated with human respiratory disease in
            China
  JOURNAL   Nature 579 (7798), 265-269 (2020)
   PUBMED   32015508
  REMARK    Erratum:[Nature. 2020 Apr;580(7803):E7. PMID: 32296181]
REFERENCE   2  (residues 1 to 4405)
  CONSRTM   NCBI Genome Project
  TITLE     Direct Submission
  JOURNAL   Submitted (17-JAN-2020) National Center for Biotechnology
            Information, NIH, Bethesda, MD 20894, USA
REFERENCE   3  (residues 1 to 4405)
  AUTHORS   Wu,F., Zhao,S., Yu,B., Chen,Y.-M., Wang,W., Hu,Y., Song,Z.-G.,
            Tao,Z.-W., Tian,J.-H., Pei,Y.-Y., Yuan,M.L., Zhang,Y.-L.,
            Dai,F.-H., Liu,Y., Wang,Q.-M., Zheng,J.-J., Xu,L., Holmes,E.C. and
            Zhang,Y.-Z.
  TITLE     Direct Submission
  JOURNAL   Submitted (05-JAN-2020) Shanghai Public Health Clinical Center &
            School of Public Health, Fudan University, Shanghai, China
COMMENT     PROVISIONAL REFSEQ: This record has not yet been subject to final
            NCBI review. The reference sequence is identical to GU280_gp01.
            Annotation was added using homology to SARSr-CoV NC_004718.3. ###
            Formerly called 'Wuhan seafood market pneumonia virus.' If you have
            questions or suggestions, please email us at info@ncbi.nlm.nih.gov
            and include the accession number NC_045512.### Protein structures
            can be found at
            https://www.ncbi.nlm.nih.gov/structure/?term=sars-cov-2.### Find
            all other Severe acute respiratory syndrome coronavirus 2
            (SARS-CoV-2) sequences at
            https://www.ncbi.nlm.nih.gov/genbank/sars-cov-2-seqs/
            
            ##Assembly-Data-START##
            Assembly Method       :: Megahit v. V1.1.3
            Sequencing Technology :: Illumina
            ##Assembly-Data-END##
            COMPLETENESS: full length.
FEATURES             Location/Qualifiers
     source          1..4405
                     /organism="Severe acute respiratory syndrome coronavirus
                     2"
                     /isolate="Wuhan-Hu-1"
                     /host="Homo sapiens"
                     /db_xref="taxon:2697049"
                     /country="China"
                     /collection_date="Dec-2019"
     Protein         1..4405
                     /product="ORF1a polyprotein"
                     /calculated_mol_wt=489861
     mat_peptide     1..180
                     /product="leader protein"
                     /note="nsp1; produced by both pp1a and pp1ab"
                     /protein_id="YP_009742608.1"
                     /calculated_mol_wt=19775
     Region          13..127
                     /region_name="SARS-CoV-like_Nsp1_N"
                     /note="N-terminal domain of non-structural protein 1 from
                     Severe acute respiratory syndrome-related coronavirus and
                     betacoronavirus in the B lineage; cd21796"
                     /db_xref="CDD:409335"
     mat_peptide     181..818
                     /product="nsp2"
                     /note="produced by both pp1a and pp1ab"
                     /protein_id="YP_009742609.1"
                     /calculated_mol_wt=70512
     Region          182..818
                     /region_name="cv_beta_Nsp2_SARS-like"
                     /note="betacoronavirus non-structural protein 2 (Nsp2)
                     similar to SARS-CoV Nsp2, and related proteins from
                     betacoronaviruses in the B lineage; cd21516"
                     /db_xref="CDD:394867"
     mat_peptide     819..2763
                     /product="nsp3"
                     /note="former nsp1; conserved domains are: N-terminal
                     acidic (Ac), predicted phosphoesterase, papain-like
                     proteinase, Y-domain, transmembrane domain 1 (TM1),
                     adenosine diphosphate-ribose 1''-phosphatase (ADRP);
                     produced by both pp1a and pp1ab"
                     /protein_id="YP_009742610.1"
                     /calculated_mol_wt=217254
     Region          880..1050
                     /region_name="DUF3655"
                     /note="Protein of unknown function (DUF3655); pfam12379"
                     /db_xref="CDD:403549"
     Region          1054..1177
                     /region_name="Macro_X_Nsp3-like"
                     /note="X-domain of viral non-structural protein 3 and
                     related macrodomains; cd21557"
                     /db_xref="CDD:394882"
     Site            order(1060,1072,1074,1119,1147..1154)
                     /site_type="other"
                     /note="putative ADP-ribose binding site [chemical
                     binding]"
                     /db_xref="CDD:394882"
     Region          1233..1358
                     /region_name="Macro_SF"
                     /note="macrodomain superfamily; cl00019"
                     /db_xref="CDD:412115"
     Region          1351..1493
                     /region_name="SUD-M"
                     /note="Single-stranded poly(A) binding domain; pfam11633"
                     /db_xref="CDD:314498"
     Region          1496..1561
                     /region_name="SUD_C_SARS-CoV_Nsp3"
                     /note="C-terminal SARS-Unique Domain (SUD) of
                     non-structural protein 3 (Nsp3) from Severe Acute
                     Respiratory Syndrome coronavirus and related
                     betacoronaviruses in the B lineage; cd21525"
                     /db_xref="CDD:394841"
     Region          1566..1868
                     /region_name="betaCoV_PLPro"
                     /note="betacoronavirus papain-like protease; cd21732"
                     /db_xref="CDD:409649"
     Site            order(1672..1675,1725..1727,1729..1730,1733,1762,
                     1785..1786,1788,1811,1827,1834..1836,1864)
                     /site_type="other"
                     /note="ubiquitin binding site [polypeptide binding]"
                     /db_xref="CDD:409649"
     Site            order(1725..1727,1810..1811,1827,1830,1836,1864)
                     /site_type="other"
                     /note="polypeptide substrate binding site [polypeptide
                     binding]"
                     /db_xref="CDD:409649"
     Region          1913..2019
                     /region_name="SARS-CoV-like_Nsp3_NAB"
                     /note="nucleic acid binding domain of non-structural
                     protein 3 from Severe acute respiratory syndrome-related
                     coronavirus and betacoronavirus in the B lineage; cd21822"
                     /db_xref="CDD:409348"
     Region          2044..2159
                     /region_name="betaCoV_Nsp3_betaSM"
                     /note="betacoronavirus-specific marker of betacoronavirus
                     non-structural protein 3; cl41743"
                     /db_xref="CDD:425374"
     Region          2232..2762
                     /region_name="TM_Y_SARS-CoV-like_Nsp3_C"
                     /note="C-terminus of non-structural protein 3, including
                     transmembrane and Y domains, from Severe acute respiratory
                     syndrome-related coronavirus and betacoronavirus in the B
                     lineage; cd21717"
                     /db_xref="CDD:409665"
     Region          2232..2253
                     /region_name="TM1"
                     /note="TM1 [structural motif]"
                     /db_xref="CDD:409665"
     Region          2337..2359
                     /region_name="TM2"
                     /note="TM2 [structural motif]"
                     /db_xref="CDD:409665"
     mat_peptide     2764..3263
                     /product="nsp4"
                     /note="nsp4B_TM; contains transmembrane domain 2 (TM2);
                     produced by both pp1a and pp1ab"
                     /protein_id="YP_009742611.1"
                     /calculated_mol_wt=56184
     Region          2777..3157
                     /region_name="cv_Nsp4_TM"
                     /note="coronavirus non-structural protein 4 (Nsp4)
                     transmembrane domain; cd21473"
                     /db_xref="CDD:394836"
     Region          2777..2799
                     /region_name="putative TM helix 1"
                     /note="putative TM helix 1 [structural motif]"
                     /db_xref="CDD:394836"
     Region          3043..3064
                     /region_name="putative TM helix 2"
                     /note="putative TM helix 2 [structural motif]"
                     /db_xref="CDD:394836"
     Region          3079..3100
                     /region_name="putative TM helix 3"
                     /note="putative TM helix 3 [structural motif]"
                     /db_xref="CDD:394836"
     Region          3127..3149
                     /region_name="putative TM helix 4"
                     /note="putative TM helix 4 [structural motif]"
                     /db_xref="CDD:394836"
     Region          3169..3261
                     /region_name="Corona_NSP4_C"
                     /note="Coronavirus nonstructural protein 4 C-terminus;
                     pfam16348"
                     /db_xref="CDD:406690"
     mat_peptide     3264..3569
                     /product="3C-like proteinase"
                     /note="nsp5A_3CLpro and nsp5B_3CLpro; main proteinase
                     (Mpro); mediates cleavages downstream of nsp4. 3D
                     structure of the SARSr-CoV homolog has been determined
                     (Yang et al., 2003); produced by both pp1a and pp1ab"
                     /protein_id="YP_009742612.1"
                     /calculated_mol_wt=33797
     Region          3267..3563
                     /region_name="betaCoV_Nsp5_Mpro"
                     /note="betacoronavirus non-structural protein 5, also
                     called Main protease (Mpro); cd21666"
                     /db_xref="CDD:394887"
     Site            order(3267..3274,3277,3381,3385..3391,3400..3404,3429,
                     3435,3549,3553,3561..3562)
                     /site_type="other"
                     /note="homodimer interface [polypeptide binding]"
                     /db_xref="CDD:394887"
     Site            order(3288..3290,3304,3403..3408,3426..3429,3435,
                     3450..3455)
                     /site_type="other"
                     /note="polypeptide substrate binding site [polypeptide
                     binding]"
                     /db_xref="CDD:394887"
     mat_peptide     3570..3859
                     /product="nsp6"
                     /note="nsp6_TM; putative transmembrane domain; produced by
                     both pp1a and pp1ab"
                     /protein_id="YP_009742613.1"
                     /calculated_mol_wt=33034
     Region          3570..3859
                     /region_name="betaCoV-Nsp6"
                     /note="betacoronavirus non-structural protein 6; cd21560"
                     /db_xref="CDD:394846"
     mat_peptide     3860..3942
                     /product="nsp7"
                     /note="produced by both pp1a and pp1ab"
                     /protein_id="YP_009742614.1"
                     /calculated_mol_wt=9240
     Region          3860..3942
                     /region_name="betaCoV_Nsp7"
                     /note="betacoronavirus non-structural protein 7; cd21827"
                     /db_xref="CDD:409253"
     Site            order(3861,3864..3867,3870..3872,3874..3875,3878,3887,
                     3890,3896,3908..3913,3915..3920,3927..3931)
                     /site_type="other"
                     /note="oligomer interface [polypeptide binding]"
                     /db_xref="CDD:409253"
     mat_peptide     3943..4140
                     /product="nsp8"
                     /note="produced by both pp1a and pp1ab"
                     /protein_id="YP_009742615.1"
                     /calculated_mol_wt=21881
     Region          3943..4139
                     /region_name="nsp8"
                     /note="nsp8 replicase; pfam08717"
                     /db_xref="CDD:400866"
     mat_peptide     4141..4253
                     /product="nsp9"
                     /note="ssRNA-binding protein; produced by both pp1a and
                     pp1ab"
                     /protein_id="YP_009742616.1"
                     /calculated_mol_wt=12378
     Region          4141..4253
                     /region_name="betaCoV_Nsp9"
                     /note="betacoronavirus non-structural protein 9; cd21898"
                     /db_xref="CDD:409331"
     Site            4141..4146
                     /site_type="other"
                     /note="N-finger"
                     /db_xref="CDD:409331"
     Site            order(4143..4144,4147..4148,4213..4214,4236..4237,
                     4239..4241,4243..4245,4247..4248)
                     /site_type="other"
                     /note="homodimer interface [polypeptide binding]"
                     /db_xref="CDD:409331"
     mat_peptide     4254..4392
                     /product="nsp10"
                     /note="nsp10_CysHis; formerly known as growth-factor-like
                     protein (GFL); produced by both pp1a and pp1ab"
                     /protein_id="YP_009742617.1"
                     /calculated_mol_wt=14790
     Region          4254..4384
                     /region_name="alpha_betaCoV_Nsp10"
                     /note="alphacoronavirus and betacoronavirus non-structural
                     protein 14; cd21901"
                     /db_xref="CDD:409326"
     Site            order(4254..4261,4265,4267..4269,4271..4273,4278..4279,
                     4282..4283,4286,4293..4298,4311..4312,4322,4324..4325,
                     4329,4331..4336,4341..4343,4346..4349)
                     /site_type="other"
                     /note="Nsp14 interface [polypeptide binding]"
                     /db_xref="CDD:409326"
     Site            order(4267..4269,4271..4273,4278,4293,4295..4298,
                     4311..4313,4331..4334,4337,4348..4349,4368)
                     /site_type="other"
                     /note="oligomer interface [polypeptide binding]"
                     /db_xref="CDD:409326"
     Site            order(4278,4295..4298,4311..4313,4337,4348..4349,4368)
                     /site_type="other"
                     /note="homotrimer interface [polypeptide binding]"
                     /db_xref="CDD:409326"
     Site            order(4293..4300,4310..4312,4322..4325,4330..4331,4333,
                     4346..4349)
                     /site_type="other"
                     /note="Nsp16 interface [polypeptide binding]"
                     /db_xref="CDD:409326"
     mat_peptide     4393..4405
                     /product="nsp11"
                     /note="produced by pp1a only"
                     /protein_id="YP_009725312.1"
                     /calculated_mol_wt=1326
     CDS             1..4405
                     /gene="ORF1ab"
                     /locus_tag="GU280_gp01"
                     /coded_by="NC_045512.2:266..13483"
                     /note="pp1a"
                     /db_xref="GeneID:43740578"
ORIGIN      
        1 meslvpgfne kthvqlslpv lqvrdvlvrg fgdsveevls earqhlkdgt cglvevekgv
       61 lpqleqpyvf ikrsdartap hghvmvelva elegiqygrs getlgvlvph vgeipvayrk
      121 vllrkngnkg agghsygadl ksfdlgdelg tdpyedfqen wntkhssgvt relmrelngg
      181 aytryvdnnf cgpdgyplec ikdllaragk asctlseqld fidtkrgvyc creheheiaw
      241 yterseksye lqtpfeikla kkfdtfngec pnfvfplnsi iktiqprvek kkldgfmgri
      301 rsvypvaspn ecnqmclstl mkcdhcgets wqtgdfvkat cefcgtenlt kegattcgyl
      361 pqnavvkiyc pachnsevgp ehslaeyhne sglktilrkg grtiafggcv fsyvgchnkc
      421 aywvprasan igcnhtgvvg egseglndnl leilqkekvn inivgdfkln eeiaiilasf
      481 sastsafvet vkgldykafk qivescgnfk vtkgkakkga wnigeqksil splyafasea
      541 arvvrsifsr tletaqnsvr vlqkaaitil dgisqyslrl idammftsdl atnnlvvmay
      601 itggvvqlts qwltnifgtv yeklkpvldw leekfkegve flrdgweivk fistcaceiv
      661 ggqivtcake ikesvqtffk lvnkflalca dsiiiggakl kalnlgetfv thskglyrkc
      721 vksreetgll mplkapkeii flegetlpte vlteevvlkt gdlqpleqpt seaveaplvg
      781 tpvcinglml leikdtekyc alapnmmvtn ntftlkggap tkvtfgddtv ievqgyksvn
      841 itfelderid kvlnekcsay tvelgtevne facvvadavi ktlqpvsell tplgidldew
      901 smatyylfde sgefklashm ycsfyppded eeegdceeee fepstqyeyg teddyqgkpl
      961 efgatsaalq peeeqeedwl dddsqqtvgq qdgsednqtt tiqtivevqp qlemeltpvv
     1021 qtievnsfsg ylkltdnvyi knadiveeak kvkptvvvna anvylkhggg vagalnkatn
     1081 namqvesddy iatngplkvg gscvlsghnl akhclhvvgp nvnkgediql lksayenfnq
     1141 hevllaplls agifgadpih slrvcvdtvr tnvylavfdk nlydklvssf lemksekqve
     1201 qkiaeipkee vkpfiteskp sveqrkqddk kikacveevt ttleetkflt enlllyidin
     1261 gnlhpdsatl vsdiditflk kdapyivgdv vqegvltavv iptkkaggtt emlakalrkv
     1321 ptdnyittyp gqglngytve eaktvlkkck safyilpsii snekqeilgt vswnlremla
     1381 haeetrklmp vcvetkaivs tiqrkykgik iqegvvdyga rfyfytsktt vaslintlnd
     1441 lnetlvtmpl gyvthglnle eaarymrslk vpatvsvssp davtayngyl tsssktpeeh
     1501 fietislags ykdwsysgqs tqlgieflkr gdksvyytsn pttfhldgev itfdnlktll
     1561 slrevrtikv fttvdninlh tqvvdmsmty gqqfgptyld gadvtkikph nshegktfyv
     1621 lpnddtlrve afeyyhttdp sflgrymsal nhtkkwkypq vngltsikwa dnncylatal
     1681 ltlqqielkf nppalqdayy rarageaanf calilaycnk tvgelgdvre tmsylfqhan
     1741 ldsckrvlnv vcktcgqqqt tlkgveavmy mgtlsyeqfk kgvqipctcg kqatkylvqq
     1801 espfvmmsap paqyelkhgt ftcaseytgn yqcghykhit sketlycidg alltksseyk
     1861 gpitdvfyke nsytttikpv tykldgvvct eidpkldnyy kkdnsyfteq pidlvpnqpy
     1921 pnasfdnfkf vcdnikfadd lnqltgykkp asrelkvtff pdlngdvvai dykhytpsfk
     1981 kgakllhkpi vwhvnnatnk atykpntwci rclwstkpve tsnsfdvlks edaqgmdnla
     2041 cedlkpvsee vvenptiqkd vlecnvktte vvgdiilkpa nnslkiteev ghtdlmaayv
     2101 dnssltikkp nelsrvlglk tlathglaav nsvpwdtian yakpflnkvv stttnivtrc
     2161 lnrvctnymp yfftlllqlc tftrstnsri kasmpttiak ntvksvgkfc leasfnylks
     2221 pnfsklinii iwflllsvcl gsliystaal gvlmsnlgmp syctgyregy lnstnvtiat
     2281 yctgsipcsv clsgldsldt ypsletiqit issfkwdlta fglvaewfla yilftrffyv
     2341 lglaaimqlf fsyfavhfis nswlmwliin lvqmapisam vrmyiffasf yyvwksyvhv
     2401 vdgcnsstcm mcykrnratr vecttivngv rrsfyvyang gkgfcklhnw ncvncdtfca
     2461 gstfisdeva rdlslqfkrp inptdqssyi vdsvtvkngs ihlyfdkagq ktyerhslsh
     2521 fvnldnlran ntkgslpinv ivfdgkskce essaksasvy ysqlmcqpil lldqalvsdv
     2581 gdsaevavkm fdayvntfss tfnvpmeklk tlvataeael aknvsldnvl stfisaarqg
     2641 fvdsdvetkd vveclklshq sdievtgdsc nnymltynkv enmtprdlga cidcsarhin
     2701 aqvakshnia liwnvkdfms lseqlrkqir saakknnlpf kltcattrqv vnvvttkial
     2761 kggkivnnwl kqlikvtlvf lfvaaifyli tpvhvmskht dfsseiigyk aidggvtrdi
     2821 astdtcfank hadfdtwfsq rggsytndka cpliaavitr evgfvvpglp gtilrttngd
     2881 flhflprvfs avgnicytps klieytdfat sacvlaaect ifkdasgkpv pycydtnvle
     2941 gsvayeslrp dtryvlmdgs iiqfpntyle gsvrvvttfd seycrhgtce rseagvcvst
     3001 sgrwvlnndy yrslpgvfcg vdavnlltnm ftpliqpiga ldisasivag givaivvtcl
     3061 ayyfmrfrra fgeyshvvaf ntllflmsft vlcltpvysf lpgvysviyl yltfyltndv
     3121 sflahiqwmv mftplvpfwi tiayiicist khfywffsny lkrrvvfngv sfstfeeaal
     3181 ctfllnkemy lklrsdvllp ltqynrylal ynkykyfsga mdttsyreaa cchlakalnd
     3241 fsnsgsdvly qppqtsitsa vlqsgfrkma fpsgkvegcm vqvtcgtttl nglwlddvvy
     3301 cprhvictse dmlnpnyedl lirksnhnfl vqagnvqlrv ighsmqncvl klkvdtanpk
     3361 tpkykfvriq pgqtfsvlac yngspsgvyq camrpnftik gsflngscgs vgfnidydcv
     3421 sfcymhhmel ptgvhagtdl egnfygpfvd rqtaqaagtd ttitvnvlaw lyaavingdr
     3481 wflnrftttl ndfnlvamky nyepltqdhv dilgplsaqt giavldmcas lkellqngmn
     3541 grtilgsall edeftpfdvv rqcsgvtfqs avkrtikgth hwllltilts llvlvqstqw
     3601 slffflyena flpfamgiia msafammfvk hkhaflclfl lpslatvayf nmvympaswv
     3661 mrimtwldmv dtslsgfklk dcvmyasavv llilmtartv yddgarrvwt lmnvltlvyk
     3721 vyygnaldqa ismwaliisv tsnysgvvtt vmflargivf mcveycpiff itgntlqcim
     3781 lvycflgyfc tcyfglfcll nryfrltlgv ydylvstqef rymnsqgllp pknsidafkl
     3841 nikllgvggk pcikvatvqs kmsdvkctsv vllsvlqqlr vesssklwaq cvqlhndill
     3901 akdtteafek mvsllsvlls mqgavdinkl ceemldnrat lqaiasefss lpsyaafata
     3961 qeayeqavan gdsevvlkkl kkslnvakse fdrdaamqrk lekmadqamt qmykqarsed
     4021 krakvtsamq tmlftmlrkl dndalnniin nardgcvpln iiplttaakl mvvipdynty
     4081 kntcdgttft yasalweiqq vvdadskivq lseismdnsp nlawplivta lransavklq
     4141 nnelspvalr qmscaagttq tactddnala yynttkggrf vlallsdlqd lkwarfpksd
     4201 gtgtiytele ppcrfvtdtp kgpkvkylyf ikglnnlnrg mvlgslaatv rlqagnatev
     4261 panstvlsfc afavdaakay kdylasggqp itncvkmlct htgtgqaitv tpeanmdqes
     4321 fggascclyc rchidhpnpk gfcdlkgkyv qipttcandp vgftlkntvc tvcgmwkgyg
     4381 cscdqlrepm lqsadaqsfl ngfav
//

