############################################
## -------------------------------------- ##
## ------------- Description ------------ ##
## -------------------------------------- ##
############################################

#Output Notes:
#All fastq data files that correspond to samples with index ending in 23 and 25 as indicated in metadata file are control reads and will be put in a directory ./Raw_Control_Reads
 #Sample file names will be changed to be "SampleID"_forward.fq (e.g. for sample HV_001_25 reads it would be ./Raw_Control_Reads/HV_001_25_forward.fq and ./Raw_Control_Reads/HV_001_25_reverse.fq)
#All other fastq data files that are not control reads will be put in a different directory ./Raw_Reads
 #Sample file names will be changed to be "SampleID"_forward.fq (e.g. for sample HV_001_01 reads it would be ./Raw_Reads/HV_001_01_forward.fq and ./Raw_Reads/HV_001_01_reverse.fq)
#The next step in this pipeline is The next step in pipeline is 2. Virome Assembly (2_Virome_Assembly.sh)

#General Notes:
#This pipeline is designed to be run using the Holland Computing Center at the University of Nebraska. Some tool commands may differ depending on installation of the tool. Please refer to the listed Githubs for each tool used as mentioned in script for further information if issues arise 
#This script is to help organize all files so that they are in an easier format to process through the pipeline. 
#Warning: Files need to be downloaded from the SRA website (Bioproject PRJNA754140) which is >4TB worth of data and will take some time to download.  

#############################################
## --------------------------------------- ##
## ----------- Make Directories ---------- ##
## --------------------------------------- ##
#############################################

mkdir Raw_Reads
mkdir Raw_Control_Reads

#####################################################
## ----------------------------------------------- ##
## - Obtain Fastq Data Files Using SRA Accession - ##
## ----------------------------------------------- ##
#####################################################

#The bioproject that contains all files is PRJNA754140. Samples HV_001 - HV_030 were used for this manuscript.
#There are many ways to download the data from SRA. You can either manually download from the SRA website or use fastq-dump from the sratoolkit.
#manuals on fastq-dump and strtoolkit can be found here: https://hpc.nih.gov/apps/sratoolkit.html

fastq-dump SRR15442278
fastq-dump SRR15442277
fastq-dump SRR15442533
fastq-dump SRR15442813
fastq-dump SRR15442456
fastq-dump SRR15442217
fastq-dump SRR15442505
fastq-dump SRR15442721
fastq-dump SRR15442300
fastq-dump SRR15442289
fastq-dump SRR15442276
fastq-dump SRR15442265
fastq-dump SRR15442254
fastq-dump SRR15442243
fastq-dump SRR15442482
fastq-dump SRR15442652
fastq-dump SRR15442609
fastq-dump SRR15442598
fastq-dump SRR15442587
fastq-dump SRR15442544
fastq-dump SRR15442532
fastq-dump SRR15442521
fastq-dump SRR15442965
fastq-dump SRR15442954
fastq-dump SRR15442943
fastq-dump SRR15442932
fastq-dump SRR15442921
fastq-dump SRR15442910
fastq-dump SRR15442835
fastq-dump SRR15442824
fastq-dump SRR15442812
fastq-dump SRR15442801
fastq-dump SRR15442790
fastq-dump SRR15442779
fastq-dump SRR15442768
fastq-dump SRR15442757
fastq-dump SRR15442746
fastq-dump SRR15442703
fastq-dump SRR15442692
fastq-dump SRR15442681
fastq-dump SRR15442455
fastq-dump SRR15442444
fastq-dump SRR15442401
fastq-dump SRR15442390
fastq-dump SRR15442379
fastq-dump SRR15442336
fastq-dump SRR15442325
fastq-dump SRR15442314
fastq-dump SRR15442239
fastq-dump SRR15442228
fastq-dump SRR15442216
fastq-dump SRR15442205
fastq-dump SRR15442194
fastq-dump SRR15442183
fastq-dump SRR15442635
fastq-dump SRR15442624
fastq-dump SRR15442613
fastq-dump SRR15442570
fastq-dump SRR15442559
fastq-dump SRR15442548
fastq-dump SRR15442504
fastq-dump SRR15442493
fastq-dump SRR15442481
fastq-dump SRR15442894
fastq-dump SRR15442883
fastq-dump SRR15442872
fastq-dump SRR15442861
fastq-dump SRR15442850
fastq-dump SRR15442743
fastq-dump SRR15442732
fastq-dump SRR15442720
fastq-dump SRR15442677
fastq-dump SRR15442666
fastq-dump SRR15442473
fastq-dump SRR15442430
fastq-dump SRR15442419
fastq-dump SRR15442408
fastq-dump SRR15442365
fastq-dump SRR15442354
fastq-dump SRR15442343
fastq-dump SRR15442299
fastq-dump SRR15442298
fastq-dump SRR15442297
fastq-dump SRR15442296
fastq-dump SRR15442295
fastq-dump SRR15442294
fastq-dump SRR15442293
fastq-dump SRR15442292
fastq-dump SRR15442291
fastq-dump SRR15442290
fastq-dump SRR15442288
fastq-dump SRR15442287
fastq-dump SRR15442286
fastq-dump SRR15442285
fastq-dump SRR15442284
fastq-dump SRR15442283
fastq-dump SRR15442282
fastq-dump SRR15442281
fastq-dump SRR15442280
fastq-dump SRR15442279
fastq-dump SRR15442275
fastq-dump SRR15442274
fastq-dump SRR15442273
fastq-dump SRR15442272
fastq-dump SRR15442271
fastq-dump SRR15442270
fastq-dump SRR15442269
fastq-dump SRR15442268
fastq-dump SRR15442267
fastq-dump SRR15442266
fastq-dump SRR15442264
fastq-dump SRR15442263
fastq-dump SRR15442262
fastq-dump SRR15442261
fastq-dump SRR15442260
fastq-dump SRR15442259
fastq-dump SRR15442258
fastq-dump SRR15442257
fastq-dump SRR15442256
fastq-dump SRR15442255
fastq-dump SRR15442253
fastq-dump SRR15442252
fastq-dump SRR15442251
fastq-dump SRR15442250
fastq-dump SRR15442249
fastq-dump SRR15442248
fastq-dump SRR15442247
fastq-dump SRR15442246
fastq-dump SRR15442245
fastq-dump SRR15442244
fastq-dump SRR15442178
fastq-dump SRR15442177
fastq-dump SRR15442176
fastq-dump SRR15442175
fastq-dump SRR15442174
fastq-dump SRR15442173
fastq-dump SRR15442172
fastq-dump SRR15442171
fastq-dump SRR15442170
fastq-dump SRR15442169
fastq-dump SRR15442662
fastq-dump SRR15442661
fastq-dump SRR15442660
fastq-dump SRR15442659
fastq-dump SRR15442658
fastq-dump SRR15442657
fastq-dump SRR15442656
fastq-dump SRR15442655
fastq-dump SRR15442654
fastq-dump SRR15442653
fastq-dump SRR15442651
fastq-dump SRR15442650
fastq-dump SRR15442649
fastq-dump SRR15442648
fastq-dump SRR15442647
fastq-dump SRR15442646
fastq-dump SRR15442645
fastq-dump SRR15442644
fastq-dump SRR15442643
fastq-dump SRR15442642
fastq-dump SRR15442608
fastq-dump SRR15442607
fastq-dump SRR15442606
fastq-dump SRR15442605
fastq-dump SRR15442604
fastq-dump SRR15442603
fastq-dump SRR15442602
fastq-dump SRR15442601
fastq-dump SRR15442600
fastq-dump SRR15442599
fastq-dump SRR15442597
fastq-dump SRR15442596
fastq-dump SRR15442595
fastq-dump SRR15442594
fastq-dump SRR15442593
fastq-dump SRR15442592
fastq-dump SRR15442591
fastq-dump SRR15442590
fastq-dump SRR15442589
fastq-dump SRR15442588
fastq-dump SRR15442586
fastq-dump SRR15442585
fastq-dump SRR15442584
fastq-dump SRR15442583
fastq-dump SRR15442582
fastq-dump SRR15442581
fastq-dump SRR15442580
fastq-dump SRR15442579
fastq-dump SRR15442578
fastq-dump SRR15442545
fastq-dump SRR15442543
fastq-dump SRR15442542
fastq-dump SRR15442541
fastq-dump SRR15442540
fastq-dump SRR15442539
fastq-dump SRR15442538
fastq-dump SRR15442537
fastq-dump SRR15442536
fastq-dump SRR15442535
fastq-dump SRR15442534
fastq-dump SRR15442531
fastq-dump SRR15442530
fastq-dump SRR15442529
fastq-dump SRR15442528
fastq-dump SRR15442527
fastq-dump SRR15442526
fastq-dump SRR15442525
fastq-dump SRR15442524
fastq-dump SRR15442523
fastq-dump SRR15442522
fastq-dump SRR15442520
fastq-dump SRR15442519
fastq-dump SRR15442518
fastq-dump SRR15442517
fastq-dump SRR15442516
fastq-dump SRR15442515
fastq-dump SRR15442514
fastq-dump SRR15442968
fastq-dump SRR15442967
fastq-dump SRR15442966
fastq-dump SRR15442964
fastq-dump SRR15442963
fastq-dump SRR15442962
fastq-dump SRR15442961
fastq-dump SRR15442960
fastq-dump SRR15442959
fastq-dump SRR15442958
fastq-dump SRR15442957
fastq-dump SRR15442956
fastq-dump SRR15442955
fastq-dump SRR15442953
fastq-dump SRR15442952
fastq-dump SRR15442951
fastq-dump SRR15442950
fastq-dump SRR15442949
fastq-dump SRR15442948
fastq-dump SRR15442947
fastq-dump SRR15442946
fastq-dump SRR15442945
fastq-dump SRR15442944
fastq-dump SRR15442942
fastq-dump SRR15442941
fastq-dump SRR15442940
fastq-dump SRR15442939
fastq-dump SRR15442938
fastq-dump SRR15442937
fastq-dump SRR15442936
fastq-dump SRR15442935
fastq-dump SRR15442934
fastq-dump SRR15442933
fastq-dump SRR15442931
fastq-dump SRR15442930
fastq-dump SRR15442929
fastq-dump SRR15442928
fastq-dump SRR15442927
fastq-dump SRR15442926
fastq-dump SRR15442925
fastq-dump SRR15442924
fastq-dump SRR15442923
fastq-dump SRR15442922
fastq-dump SRR15442920
fastq-dump SRR15442919
fastq-dump SRR15442918
fastq-dump SRR15442917
fastq-dump SRR15442916
fastq-dump SRR15442915
fastq-dump SRR15442914
fastq-dump SRR15442913
fastq-dump SRR15442912
fastq-dump SRR15442911
fastq-dump SRR15442909
fastq-dump SRR15442908
fastq-dump SRR15442907
fastq-dump SRR15442906
fastq-dump SRR15442905
fastq-dump SRR15442840
fastq-dump SRR15442839
fastq-dump SRR15442838
fastq-dump SRR15442837
fastq-dump SRR15442836
fastq-dump SRR15442834
fastq-dump SRR15442833
fastq-dump SRR15442832
fastq-dump SRR15442831
fastq-dump SRR15442830
fastq-dump SRR15442829
fastq-dump SRR15442828
fastq-dump SRR15442827
fastq-dump SRR15442826
fastq-dump SRR15442825
fastq-dump SRR15442823
fastq-dump SRR15442822
fastq-dump SRR15442821
fastq-dump SRR15442820
fastq-dump SRR15442819
fastq-dump SRR15442818
fastq-dump SRR15442817
fastq-dump SRR15442816
fastq-dump SRR15442815
fastq-dump SRR15442814
fastq-dump SRR15442811
fastq-dump SRR15442810
fastq-dump SRR15442809
fastq-dump SRR15442808
fastq-dump SRR15442807
fastq-dump SRR15442806
fastq-dump SRR15442805
fastq-dump SRR15442804
fastq-dump SRR15442803
fastq-dump SRR15442802
fastq-dump SRR15442800
fastq-dump SRR15442799
fastq-dump SRR15442798
fastq-dump SRR15442797
fastq-dump SRR15442796
fastq-dump SRR15442795
fastq-dump SRR15442794
fastq-dump SRR15442793
fastq-dump SRR15442792
fastq-dump SRR15442791
fastq-dump SRR15442789
fastq-dump SRR15442788
fastq-dump SRR15442787
fastq-dump SRR15442786
fastq-dump SRR15442785
fastq-dump SRR15442784
fastq-dump SRR15442783
fastq-dump SRR15442782
fastq-dump SRR15442781
fastq-dump SRR15442780
fastq-dump SRR15442778
fastq-dump SRR15442777
fastq-dump SRR15442776
fastq-dump SRR15442775
fastq-dump SRR15442774
fastq-dump SRR15442773
fastq-dump SRR15442772
fastq-dump SRR15442771
fastq-dump SRR15442770
fastq-dump SRR15442769
fastq-dump SRR15442767
fastq-dump SRR15442766
fastq-dump SRR15442765
fastq-dump SRR15442764
fastq-dump SRR15442763
fastq-dump SRR15442762
fastq-dump SRR15442761
fastq-dump SRR15442760
fastq-dump SRR15442759
fastq-dump SRR15442758
fastq-dump SRR15442756
fastq-dump SRR15442755
fastq-dump SRR15442754
fastq-dump SRR15442753
fastq-dump SRR15442752
fastq-dump SRR15442751
fastq-dump SRR15442750
fastq-dump SRR15442749
fastq-dump SRR15442748
fastq-dump SRR15442747
fastq-dump SRR15442745
fastq-dump SRR15442712
fastq-dump SRR15442711
fastq-dump SRR15442710
fastq-dump SRR15442709
fastq-dump SRR15442708
fastq-dump SRR15442707
fastq-dump SRR15442706
fastq-dump SRR15442705
fastq-dump SRR15442704
fastq-dump SRR15442702
fastq-dump SRR15442701
fastq-dump SRR15442700
fastq-dump SRR15442699
fastq-dump SRR15442698
fastq-dump SRR15442697
fastq-dump SRR15442696
fastq-dump SRR15442695
fastq-dump SRR15442694
fastq-dump SRR15442693
fastq-dump SRR15442691
fastq-dump SRR15442690
fastq-dump SRR15442689
fastq-dump SRR15442688
fastq-dump SRR15442687
fastq-dump SRR15442686
fastq-dump SRR15442685
fastq-dump SRR15442684
fastq-dump SRR15442683
fastq-dump SRR15442682
fastq-dump SRR15442466
fastq-dump SRR15442465
fastq-dump SRR15442464
fastq-dump SRR15442463
fastq-dump SRR15442462
fastq-dump SRR15442461
fastq-dump SRR15442460
fastq-dump SRR15442459
fastq-dump SRR15442458
fastq-dump SRR15442457
fastq-dump SRR15442454
fastq-dump SRR15442453
fastq-dump SRR15442452
fastq-dump SRR15442451
fastq-dump SRR15442450
fastq-dump SRR15442449
fastq-dump SRR15442448
fastq-dump SRR15442447
fastq-dump SRR15442446
fastq-dump SRR15442445
fastq-dump SRR15442443
fastq-dump SRR15442442
fastq-dump SRR15442441
fastq-dump SRR15442440
fastq-dump SRR15442439
fastq-dump SRR15442438
fastq-dump SRR15442437
fastq-dump SRR15442436
fastq-dump SRR15442435
fastq-dump SRR15442402
fastq-dump SRR15442400
fastq-dump SRR15442399
fastq-dump SRR15442398
fastq-dump SRR15442397
fastq-dump SRR15442396
fastq-dump SRR15442395
fastq-dump SRR15442394
fastq-dump SRR15442393
fastq-dump SRR15442392
fastq-dump SRR15442391
fastq-dump SRR15442389
fastq-dump SRR15442388
fastq-dump SRR15442387
fastq-dump SRR15442386
fastq-dump SRR15442385
fastq-dump SRR15442384
fastq-dump SRR15442383
fastq-dump SRR15442382
fastq-dump SRR15442381
fastq-dump SRR15442380
fastq-dump SRR15442378
fastq-dump SRR15442377
fastq-dump SRR15442376
fastq-dump SRR15442375
fastq-dump SRR15442374
fastq-dump SRR15442373
fastq-dump SRR15442372
fastq-dump SRR15442371
fastq-dump SRR15442338
fastq-dump SRR15442337
fastq-dump SRR15442335
fastq-dump SRR15442334
fastq-dump SRR15442333
fastq-dump SRR15442332
fastq-dump SRR15442331
fastq-dump SRR15442330
fastq-dump SRR15442329
fastq-dump SRR15442328
fastq-dump SRR15442327
fastq-dump SRR15442326
fastq-dump SRR15442324
fastq-dump SRR15442323
fastq-dump SRR15442322
fastq-dump SRR15442321
fastq-dump SRR15442320
fastq-dump SRR15442319
fastq-dump SRR15442318
fastq-dump SRR15442317
fastq-dump SRR15442316
fastq-dump SRR15442315
fastq-dump SRR15442313
fastq-dump SRR15442312
fastq-dump SRR15442311
fastq-dump SRR15442310
fastq-dump SRR15442309
fastq-dump SRR15442308
fastq-dump SRR15442307
fastq-dump SRR15442242
fastq-dump SRR15442241
fastq-dump SRR15442240
fastq-dump SRR15442238
fastq-dump SRR15442237
fastq-dump SRR15442236
fastq-dump SRR15442235
fastq-dump SRR15442234
fastq-dump SRR15442233
fastq-dump SRR15442232
fastq-dump SRR15442231
fastq-dump SRR15442230
fastq-dump SRR15442229
fastq-dump SRR15442227
fastq-dump SRR15442226
fastq-dump SRR15442225
fastq-dump SRR15442224
fastq-dump SRR15442223
fastq-dump SRR15442222
fastq-dump SRR15442221
fastq-dump SRR15442220
fastq-dump SRR15442219
fastq-dump SRR15442218
fastq-dump SRR15442215
fastq-dump SRR15442214
fastq-dump SRR15442213
fastq-dump SRR15442212
fastq-dump SRR15442211
fastq-dump SRR15442210
fastq-dump SRR15442209
fastq-dump SRR15442208
fastq-dump SRR15442207
fastq-dump SRR15442206
fastq-dump SRR15442204
fastq-dump SRR15442203
fastq-dump SRR15442202
fastq-dump SRR15442201
fastq-dump SRR15442200
fastq-dump SRR15442199
fastq-dump SRR15442198
fastq-dump SRR15442197
fastq-dump SRR15442196
fastq-dump SRR15442195
fastq-dump SRR15442193
fastq-dump SRR15442192
fastq-dump SRR15442191
fastq-dump SRR15442190
fastq-dump SRR15442189
fastq-dump SRR15442188
fastq-dump SRR15442187
fastq-dump SRR15442186
fastq-dump SRR15442185
fastq-dump SRR15442184
fastq-dump SRR15442182
fastq-dump SRR15442181
fastq-dump SRR15442180
fastq-dump SRR15442179
fastq-dump SRR15442641
fastq-dump SRR15442640
fastq-dump SRR15442639
fastq-dump SRR15442638
fastq-dump SRR15442637
fastq-dump SRR15442636
fastq-dump SRR15442634
fastq-dump SRR15442633
fastq-dump SRR15442632
fastq-dump SRR15442631
fastq-dump SRR15442630
fastq-dump SRR15442629
fastq-dump SRR15442628
fastq-dump SRR15442627
fastq-dump SRR15442626
fastq-dump SRR15442625
fastq-dump SRR15442623
fastq-dump SRR15442622
fastq-dump SRR15442621
fastq-dump SRR15442620
fastq-dump SRR15442619
fastq-dump SRR15442618
fastq-dump SRR15442617
fastq-dump SRR15442616
fastq-dump SRR15442615
fastq-dump SRR15442614
fastq-dump SRR15442612
fastq-dump SRR15442611
fastq-dump SRR15442610
fastq-dump SRR15442577
fastq-dump SRR15442576
fastq-dump SRR15442575
fastq-dump SRR15442574
fastq-dump SRR15442573
fastq-dump SRR15442572
fastq-dump SRR15442571
fastq-dump SRR15442569
fastq-dump SRR15442568
fastq-dump SRR15442567
fastq-dump SRR15442566
fastq-dump SRR15442565
fastq-dump SRR15442564
fastq-dump SRR15442563
fastq-dump SRR15442562
fastq-dump SRR15442561
fastq-dump SRR15442560
fastq-dump SRR15442558
fastq-dump SRR15442557
fastq-dump SRR15442556
fastq-dump SRR15442555
fastq-dump SRR15442554
fastq-dump SRR15442553
fastq-dump SRR15442552
fastq-dump SRR15442551
fastq-dump SRR15442550
fastq-dump SRR15442549
fastq-dump SRR15442547
fastq-dump SRR15442546
fastq-dump SRR15442513
fastq-dump SRR15442512
fastq-dump SRR15442511
fastq-dump SRR15442510
fastq-dump SRR15442509
fastq-dump SRR15442508
fastq-dump SRR15442507
fastq-dump SRR15442506
fastq-dump SRR15442503
fastq-dump SRR15442502
fastq-dump SRR15442501
fastq-dump SRR15442500
fastq-dump SRR15442499
fastq-dump SRR15442498
fastq-dump SRR15442497
fastq-dump SRR15442496
fastq-dump SRR15442495
fastq-dump SRR15442494
fastq-dump SRR15442492
fastq-dump SRR15442491
fastq-dump SRR15442490
fastq-dump SRR15442489
fastq-dump SRR15442488
fastq-dump SRR15442487
fastq-dump SRR15442486
fastq-dump SRR15442485
fastq-dump SRR15442484
fastq-dump SRR15442483
fastq-dump SRR15442904
fastq-dump SRR15442903
fastq-dump SRR15442902
fastq-dump SRR15442901
fastq-dump SRR15442900
fastq-dump SRR15442899
fastq-dump SRR15442898
fastq-dump SRR15442897
fastq-dump SRR15442896
fastq-dump SRR15442895
fastq-dump SRR15442893
fastq-dump SRR15442892
fastq-dump SRR15442891
fastq-dump SRR15442890
fastq-dump SRR15442889
fastq-dump SRR15442888
fastq-dump SRR15442887
fastq-dump SRR15442886
fastq-dump SRR15442885
fastq-dump SRR15442884
fastq-dump SRR15442882
fastq-dump SRR15442881
fastq-dump SRR15442880
fastq-dump SRR15442879
fastq-dump SRR15442878
fastq-dump SRR15442877
fastq-dump SRR15442876
fastq-dump SRR15442875
fastq-dump SRR15442874
fastq-dump SRR15442873
fastq-dump SRR15442871
fastq-dump SRR15442870
fastq-dump SRR15442869
fastq-dump SRR15442868
fastq-dump SRR15442867
fastq-dump SRR15442866
fastq-dump SRR15442865
fastq-dump SRR15442864
fastq-dump SRR15442863
fastq-dump SRR15442862
fastq-dump SRR15442860
fastq-dump SRR15442859
fastq-dump SRR15442858
fastq-dump SRR15442857
fastq-dump SRR15442856
fastq-dump SRR15442855
fastq-dump SRR15442854
fastq-dump SRR15442853
fastq-dump SRR15442852
fastq-dump SRR15442851
fastq-dump SRR15442849
fastq-dump SRR15442848
fastq-dump SRR15442847
fastq-dump SRR15442846
fastq-dump SRR15442845
fastq-dump SRR15442844
fastq-dump SRR15442843
fastq-dump SRR15442842
fastq-dump SRR15442841
fastq-dump SRR15442744
fastq-dump SRR15442742
fastq-dump SRR15442741
fastq-dump SRR15442740
fastq-dump SRR15442739
fastq-dump SRR15442738
fastq-dump SRR15442737
fastq-dump SRR15442736
fastq-dump SRR15442735
fastq-dump SRR15442734
fastq-dump SRR15442733
fastq-dump SRR15442731
fastq-dump SRR15442730
fastq-dump SRR15442729
fastq-dump SRR15442728
fastq-dump SRR15442727
fastq-dump SRR15442726
fastq-dump SRR15442725
fastq-dump SRR15442724
fastq-dump SRR15442723
fastq-dump SRR15442722
fastq-dump SRR15442719
fastq-dump SRR15442718
fastq-dump SRR15442717
fastq-dump SRR15442716
fastq-dump SRR15442715
fastq-dump SRR15442714
fastq-dump SRR15442713
fastq-dump SRR15442680

############################################
## -------------------------------------- ##
## ------------ Forward Reads ----------- ##
## -------------------------------------- ##
############################################

mv EG_Sflab_241_1_USPD16092227-1_HV5CHBBXX_L2_1.fq.gz Raw_Reads/HV_001_01_forward.fq.gz
mv EG_Sflab_241_2_USPD16092227-2_HV5CHBBXX_L2_1.fq.gz Raw_Reads/HV_001_02_forward.fq.gz
mv EG_Sflab_241_3_USPD16092227-3_HV5CHBBXX_L2_1.fq.gz Raw_Reads/HV_001_03_forward.fq.gz
mv EG_Sflab_241_4_USPD16092227-4_HV5CHBBXX_L2_1.fq.gz Raw_Reads/HV_001_04_forward.fq.gz
mv EG_Sflab_241_5_USPD16092227-5_HV5CHBBXX_L2_1.fq.gz Raw_Reads/HV_001_05_forward.fq.gz
mv EG_Sflab_241_6_USPD16092227-6_HV5CHBBXX_L2_1.fq.gz Raw_Reads/HV_001_06_forward.fq.gz
mv EG_Sflab_241_7_USPD16092227-7_HV5CHBBXX_L2_1.fq.gz Raw_Reads/HV_001_07_forward.fq.gz
mv EG_Sflab_241_8_USPD16092227-8_HV5CHBBXX_L2_1.fq.gz Raw_Reads/HV_001_08_forward.fq.gz
mv EG_Sflab_241_9_USPD16092227-9_HV5CHBBXX_L2_1.fq.gz Raw_Reads/HV_001_09_forward.fq.gz
mv EG_Sflab_241_10_USPD16092227-10_HV5CHBBXX_L2_1.fq.gz Raw_Reads/HV_001_10_forward.fq.gz
mv EG_Sflab_241_11_USPD16092227-11_HV5CHBBXX_L2_1.fq.gz Raw_Reads/HV_001_11_forward.fq.gz
mv EG_Sflab_241_12_USPD16092227-12_HV5CHBBXX_L2_1.fq.gz Raw_Reads/HV_001_12_forward.fq.gz
mv EG_Sflab_241_13_USPD16092227-13_HV5CHBBXX_L2_1.fq.gz Raw_Reads/HV_001_13_forward.fq.gz
mv EG_Sflab_241_14_USPD16092227-14_HV5CHBBXX_L2_1.fq.gz Raw_Reads/HV_001_14_forward.fq.gz
mv EG_Sflab_241_15_USPD16092227-15_HV5CHBBXX_L2_1.fq.gz Raw_Reads/HV_001_15_forward.fq.gz
mv EG_Sflab_241_16_USPD16092227-16_HV5CHBBXX_L2_1.fq.gz Raw_Reads/HV_001_16_forward.fq.gz
mv EG_Sflab_241_18_USPD16092227-18_HV5CHBBXX_L2_1.fq.gz Raw_Reads/HV_001_18_forward.fq.gz
mv EG_Sflab_241_19_USPD16092227-19_HV5CHBBXX_L2_1.fq.gz Raw_Reads/HV_001_19_forward.fq.gz
mv EG_Sflab_241_20_USPD16092227-20_HV5CHBBXX_L2_1.fq.gz Raw_Reads/HV_001_20_forward.fq.gz
mv EG_Sflab_241_21_USPD16092227-21_HV5CHBBXX_L2_1.fq.gz Raw_Reads/HV_001_21_forward.fq.gz
mv EG_Sflab_241_22_USPD16092227-22_HV5CHBBXX_L2_1.fq.gz Raw_Reads/HV_001_22_forward.fq.gz
mv EG_Sflab_241_23_USPD16092227-23_HV5CHBBXX_L2_1.fq.gz Raw_Reads/HV_001_23_forward.fq.gz
mv EG_Sflab_241_25_USPD16092227-25_HV5CHBBXX_L2_1.fq.gz Raw_Reads/HV_001_25_forward.fq.gz
mv EG_Sflab_241_27_USPD16092227-27_HV5CHBBXX_L2_1.fq.gz Raw_Reads/HV_001_27_forward.fq.gz
mv HV002_1_USPD16094887-1_HWJJHBBXX_L2_1.fq.gz Raw_Reads/HV_002_01_forward.fq.gz
mv HV002_2_USPD16094887-2_HWJJHBBXX_L2_1.fq.gz Raw_Reads/HV_002_02_forward.fq.gz
mv HV002_3_USPD16094887-3_HWJJHBBXX_L2_1.fq.gz Raw_Reads/HV_002_03_forward.fq.gz
mv HV002_4_USPD16094887-4_HWJJHBBXX_L2_1.fq.gz Raw_Reads/HV_002_04_forward.fq.gz
mv HV002_5_USPD16094887-5_HWJJHBBXX_L2_1.fq.gz Raw_Reads/HV_002_05_forward.fq.gz
mv HV002_6_USPD16094887-6_HWJJHBBXX_L2_1.fq.gz Raw_Reads/HV_002_06_forward.fq.gz
mv HV002_7_USPD16094887-7_HWJJHBBXX_L2_1.fq.gz Raw_Reads/HV_002_07_forward.fq.gz
mv HV002_8_USPD16094887-8_HWJJHBBXX_L2_1.fq.gz Raw_Reads/HV_002_08_forward.fq.gz
mv HV002_9_USPD16094887-9_HWJJHBBXX_L2_1.fq.gz Raw_Reads/HV_002_09_forward.fq.gz
mv HV002_10_USPD16094887-10_HWJJHBBXX_L2_1.fq.gz Raw_Reads/HV_002_10_forward.fq.gz
mv HV002_11_USPD16094887-11_HWJJHBBXX_L2_1.fq.gz Raw_Reads/HV_002_11_forward.fq.gz
mv HV002_12_USPD16094887-12_HWJJHBBXX_L2_1.fq.gz Raw_Reads/HV_002_12_forward.fq.gz
mv HV002_13_USPD16094887-13_HWJJHBBXX_L2_1.fq.gz Raw_Reads/HV_002_13_forward.fq.gz
mv HV002_14_USPD16094887-14_HWJJHBBXX_L2_1.fq.gz Raw_Reads/HV_002_14_forward.fq.gz
mv HV002_15_USPD16094887-15_HWJJHBBXX_L2_1.fq.gz Raw_Reads/HV_002_15_forward.fq.gz
mv HV002_18_USPD16094887-18_HWJJHBBXX_L2_1.fq.gz Raw_Reads/HV_002_18_forward.fq.gz
mv HV002_19_USPD16094887-19_HWJJHBBXX_L2_1.fq.gz Raw_Reads/HV_002_19_forward.fq.gz
mv HV002_20_USPD16094887-20_HWJJHBBXX_L2_1.fq.gz Raw_Reads/HV_002_20_forward.fq.gz
mv HV002_21_USPD16094887-21_HWJJHBBXX_L2_1.fq.gz Raw_Reads/HV_002_21_forward.fq.gz
mv HV002_22_USPD16094887-22_HWJJHBBXX_L2_1.fq.gz Raw_Reads/HV_002_22_forward.fq.gz
mv HV002_23_USPD16094887-23_HWJJHBBXX_L2_1.fq.gz Raw_Reads/HV_002_23_forward.fq.gz
mv HV003_1_USPD16094888-1_HWJJHBBXX_L1_1.fq.gz Raw_Reads/HV_003_01_forward.fq.gz
mv HV003_2_USPD16094888-2_HWJJHBBXX_L1_1.fq.gz Raw_Reads/HV_003_02_forward.fq.gz
mv HV003_3_USPD16094888-3_HWJJHBBXX_L1_1.fq.gz Raw_Reads/HV_003_03_forward.fq.gz
mv HV003_4_USPD16094888-4_HWJJHBBXX_L1_1.fq.gz Raw_Reads/HV_003_04_forward.fq.gz
mv HV003_5_USPD16094888-5_HWJJHBBXX_L1_1.fq.gz Raw_Reads/HV_003_05_forward.fq.gz
mv HV003_6_USPD16094888-6_HWJJHBBXX_L1_1.fq.gz Raw_Reads/HV_003_06_forward.fq.gz
mv HV003_7_USPD16094888-7_HWJJHBBXX_L1_1.fq.gz Raw_Reads/HV_003_07_forward.fq.gz
mv HV003_8_USPD16094888-8_HWJJHBBXX_L1_1.fq.gz Raw_Reads/HV_003_08_forward.fq.gz
mv HV003_9_USPD16094888-9_HWJJHBBXX_L1_1.fq.gz Raw_Reads/HV_003_09_forward.fq.gz
mv HV003_10_USPD16094888-10_HWJJHBBXX_L1_1.fq.gz Raw_Reads/HV_003_10_forward.fq.gz
mv HV003_11_USPD16094888-11_HWJJHBBXX_L1_1.fq.gz Raw_Reads/HV_003_11_forward.fq.gz
mv HV003_12_USPD16094888-12_HWJJHBBXX_L1_1.fq.gz Raw_Reads/HV_003_12_forward.fq.gz
mv HV003_13_USPD16094888-13_HWJJHBBXX_L1_1.fq.gz Raw_Reads/HV_003_13_forward.fq.gz
mv HV003_14_USPD16094888-14_HWJJHBBXX_L1_1.fq.gz Raw_Reads/HV_003_14_forward.fq.gz
mv HV003_15_USPD16094888-15_HWJJHBBXX_L1_1.fq.gz Raw_Reads/HV_003_15_forward.fq.gz
mv HV003_16_USPD16094888-16_HWJJHBBXX_L1_1.fq.gz Raw_Reads/HV_003_16_forward.fq.gz
mv HV003_18_USPD16094888-18_HWJJHBBXX_L1_1.fq.gz Raw_Reads/HV_003_18_forward.fq.gz
mv HV003_19_USPD16094888-19_HWJJHBBXX_L1_1.fq.gz Raw_Reads/HV_003_19_forward.fq.gz
mv HV003_20_USPD16094888-20_HWJJHBBXX_L1_1.fq.gz Raw_Reads/HV_003_20_forward.fq.gz
mv HV003_21_USPD16094888-21_HWJJHBBXX_L1_1.fq.gz Raw_Reads/HV_003_21_forward.fq.gz
mv HV003_22_USPD16094888-22_HWJJHBBXX_L1_1.fq.gz Raw_Reads/HV_003_22_forward.fq.gz
mv HV004_1_USPD16094889-1_HWKHHBBXX_L8_1.fq.gz Raw_Reads/HV_004_01_forward.fq.gz
mv HV004_2_USPD16094889-2_HWKHHBBXX_L8_1.fq.gz Raw_Reads/HV_004_02_forward.fq.gz
mv HV004_3_USPD16094889-3_HWKHHBBXX_L8_1.fq.gz Raw_Reads/HV_004_03_forward.fq.gz
mv HV004_4_USPD16094889-4_HWKHHBBXX_L8_1.fq.gz Raw_Reads/HV_004_04_forward.fq.gz
mv HV004_5_USPD16094889-5_HWKHHBBXX_L8_1.fq.gz Raw_Reads/HV_004_05_forward.fq.gz
mv HV004_6_USPD16094889-6_HWKHHBBXX_L8_1.fq.gz Raw_Reads/HV_004_06_forward.fq.gz
mv HV004_7_USPD16094889-7_HWKHHBBXX_L8_1.fq.gz Raw_Reads/HV_004_07_forward.fq.gz
mv HV004_8_USPD16094889-8_HWKHHBBXX_L8_1.fq.gz Raw_Reads/HV_004_08_forward.fq.gz
mv HV004_9_USPD16094889-9_HWKHHBBXX_L8_1.fq.gz Raw_Reads/HV_004_09_forward.fq.gz
mv HV004_10_USPD16094889-10_HWKHHBBXX_L8_1.fq.gz Raw_Reads/HV_004_10_forward.fq.gz
mv HV004_11_USPD16094889-11_HWKHHBBXX_L8_1.fq.gz Raw_Reads/HV_004_11_forward.fq.gz
mv HV004_12_USPD16094889-12_HWKHHBBXX_L8_1.fq.gz Raw_Reads/HV_004_12_forward.fq.gz
mv HV004_13_USPD16094889-13_HWKHHBBXX_L8_1.fq.gz Raw_Reads/HV_004_13_forward.fq.gz
mv HV004_14_USPD16094889-14_HWKHHBBXX_L8_1.fq.gz Raw_Reads/HV_004_14_forward.fq.gz
mv HV004_15_USPD16094889-15_HWKHHBBXX_L8_1.fq.gz Raw_Reads/HV_004_15_forward.fq.gz
mv HV004_16_USPD16094889-16_HWKHHBBXX_L8_1.fq.gz Raw_Reads/HV_004_16_forward.fq.gz
mv HV004_18_USPD16094889-18_HWKHHBBXX_L8_1.fq.gz Raw_Reads/HV_004_18_forward.fq.gz
mv HV004_19_USPD16094889-19_HWKHHBBXX_L8_1.fq.gz Raw_Reads/HV_004_19_forward.fq.gz
mv HV004_20_USPD16094889-20_HWKHHBBXX_L8_1.fq.gz Raw_Reads/HV_004_20_forward.fq.gz
mv HV004_21_USPD16094889-21_HWKHHBBXX_L8_1.fq.gz Raw_Reads/HV_004_21_forward.fq.gz
mv HV004_22_USPD16094889-22_HWKHHBBXX_L8_1.fq.gz Raw_Reads/HV_004_22_forward.fq.gz
mv HV004_23_USPD16094889-23_HWKHHBBXX_L8_1.fq.gz Raw_Reads/HV_004_23_forward.fq.gz
mv HV004_25_USPD16094889-25_HWKHHBBXX_L8_1.fq.gz Raw_Reads/HV_004_25_forward.fq.gz
mv HV004_27_USPD16094889-27_HWKHHBBXX_L8_1.fq.gz Raw_Reads/HV_004_27_forward.fq.gz
mv HV005_1_USPD16094891-1_HWKHHBBXX_L7_1.fq.gz Raw_Reads/HV_005_01_forward.fq.gz
mv HV005_3_USPD16094891-3_HWKHHBBXX_L7_1.fq.gz Raw_Reads/HV_005_03_forward.fq.gz
mv HV005_4_USPD16094891-4_HWKHHBBXX_L7_1.fq.gz Raw_Reads/HV_005_04_forward.fq.gz
mv HV005_5_USPD16094891-5_HWKHHBBXX_L7_1.fq.gz Raw_Reads/HV_005_05_forward.fq.gz
mv HV005_6_USPD16094891-6_HWKHHBBXX_L7_1.fq.gz Raw_Reads/HV_005_06_forward.fq.gz
mv HV005_7_USPD16094891-7_HWKHHBBXX_L7_1.fq.gz Raw_Reads/HV_005_07_forward.fq.gz
mv HV005_8_USPD16094891-8_HWKHHBBXX_L7_1.fq.gz Raw_Reads/HV_005_08_forward.fq.gz
mv HV005_9_USPD16094891-9_HWKHHBBXX_L7_1.fq.gz Raw_Reads/HV_005_09_forward.fq.gz
mv HV005_10_USPD16094891-10_HWKHHBBXX_L7_1.fq.gz Raw_Reads/HV_005_10_forward.fq.gz
mv HV005_11_USPD16094891-11_HWKHHBBXX_L7_1.fq.gz Raw_Reads/HV_005_11_forward.fq.gz
mv HV005_12_USPD16094891-12_HWKHHBBXX_L7_1.fq.gz Raw_Reads/HV_005_12_forward.fq.gz
mv HV005_13_USPD16094891-13_HWKHHBBXX_L7_1.fq.gz Raw_Reads/HV_005_13_forward.fq.gz
mv HV005_14_USPD16094891-14_HWKHHBBXX_L7_1.fq.gz Raw_Reads/HV_005_14_forward.fq.gz
mv HV005_15_USPD16094891-15_HWKHHBBXX_L7_1.fq.gz Raw_Reads/HV_005_15_forward.fq.gz
mv HV005_16_USPD16094891-16_HWKHHBBXX_L7_1.fq.gz Raw_Reads/HV_005_16_forward.fq.gz
mv HV005_18_USPD16094891-18_HWKHHBBXX_L7_1.fq.gz Raw_Reads/HV_005_18_forward.fq.gz
mv HV005_19_USPD16094891-19_HWKHHBBXX_L7_1.fq.gz Raw_Reads/HV_005_19_forward.fq.gz
mv HV005_20_USPD16094891-20_HWKHHBBXX_L7_1.fq.gz Raw_Reads/HV_005_20_forward.fq.gz
mv HV005_22_USPD16094891-22_HWKHHBBXX_L7_1.fq.gz Raw_Reads/HV_005_22_forward.fq.gz
mv HV005_23_USPD16094891-23_HWKHHBBXX_L7_1.fq.gz Raw_Reads/HV_005_23_forward.fq.gz
mv HV005_25_USPD16094891-25_HWKHHBBXX_L7_1.fq.gz Raw_Reads/HV_005_25_forward.fq.gz
mv HV005_27_USPD16094891-27_HWKHHBBXX_L7_1.fq.gz Raw_Reads/HV_005_27_forward.fq.gz
mv HV_006_1_USPD16099387-1_H2FJ5BBXX_L6_1.fq.gz Raw_Reads/HV_006_01_forward.fq.gz
mv HV_006_2_USPD16099387-2_H2FJ5BBXX_L6_1.fq.gz Raw_Reads/HV_006_02_forward.fq.gz
mv HV_006_3_USPD16099387-3_H2FJ5BBXX_L6_1.fq.gz Raw_Reads/HV_006_03_forward.fq.gz
mv HV_006_4_USPD16099387-4_H2FJ5BBXX_L6_1.fq.gz Raw_Reads/HV_006_04_forward.fq.gz
mv HV_006_5_USPD16099387-5_H2FJ5BBXX_L6_1.fq.gz Raw_Reads/HV_006_05_forward.fq.gz
mv HV_006_6_USPD16099387-6_H2FJ5BBXX_L6_1.fq.gz Raw_Reads/HV_006_06_forward.fq.gz
mv HV_006_7_USPD16099387-7_H2FJ5BBXX_L6_1.fq.gz Raw_Reads/HV_006_07_forward.fq.gz
mv HV_006_8_USPD16099387-8_H2FJ5BBXX_L6_1.fq.gz Raw_Reads/HV_006_08_forward.fq.gz
mv HV_006_9_USPD16099387-9_H2FJ5BBXX_L6_1.fq.gz Raw_Reads/HV_006_09_forward.fq.gz
mv HV_006_10_USPD16099387-10_H2FJ5BBXX_L6_1.fq.gz Raw_Reads/HV_006_10_forward.fq.gz
mv HV_006_11_USPD16099387-11_H2FJ5BBXX_L6_1.fq.gz Raw_Reads/HV_006_11_forward.fq.gz
mv HV_006_12_USPD16099387-12_H2FJ5BBXX_L6_1.fq.gz Raw_Reads/HV_006_12_forward.fq.gz
mv HV_006_13_USPD16099387-13_H2FJ5BBXX_L6_1.fq.gz Raw_Reads/HV_006_13_forward.fq.gz
mv HV_006_14_USPD16099387-14_H2FJ5BBXX_L6_1.fq.gz Raw_Reads/HV_006_14_forward.fq.gz
mv HV_006_15_USPD16099387-15_H2FJ5BBXX_L6_1.fq.gz Raw_Reads/HV_006_15_forward.fq.gz
mv HV_006_16_USPD16099387-16_H2FJ5BBXX_L6_1.fq.gz Raw_Reads/HV_006_16_forward.fq.gz
mv HV_006_18_USPD16099387-18_H2FJ5BBXX_L6_1.fq.gz Raw_Reads/HV_006_18_forward.fq.gz
mv HV_006_19_USPD16099387-19_H2FJ5BBXX_L6_1.fq.gz Raw_Reads/HV_006_19_forward.fq.gz
mv HV_006_20_USPD16099387-20_H2FJ5BBXX_L6_1.fq.gz Raw_Reads/HV_006_20_forward.fq.gz
mv HV_006_21_USPD16099387-21_H2FJ5BBXX_L6_1.fq.gz Raw_Reads/HV_006_21_forward.fq.gz
mv HV_006_22_USPD16099387-22_H2FJ5BBXX_L6_1.fq.gz Raw_Reads/HV_006_22_forward.fq.gz
mv HV_006_23_USPD16099387-23_H2FJ5BBXX_L6_1.fq.gz Raw_Reads/HV_006_23_forward.fq.gz
mv HV_006_25_USPD16099387-25_H2FJ5BBXX_L6_1.fq.gz Raw_Reads/HV_006_25_forward.fq.gz
mv HV_006_27_USPD16099387-27_H2FJ5BBXX_L6_1.fq.gz Raw_Reads/HV_006_27_forward.fq.gz
mv HV_007_1_USPD16099388-1_H2FJ5BBXX_L7_1.fq.gz Raw_Reads/HV_007_01_forward.fq.gz
mv HV_007_2_USPD16099388-2_H2FJ5BBXX_L7_1.fq.gz Raw_Reads/HV_007_02_forward.fq.gz
mv HV_007_3_USPD16099388-3_H2FJ5BBXX_L7_1.fq.gz Raw_Reads/HV_007_03_forward.fq.gz
mv HV_007_4_USPD16099388-4_H2FJ5BBXX_L7_1.fq.gz Raw_Reads/HV_007_04_forward.fq.gz
mv HV_007_5_USPD16099388-5_H2FJ5BBXX_L7_1.fq.gz Raw_Reads/HV_007_05_forward.fq.gz
mv HV_007_6_USPD16099388-6_H2FJ5BBXX_L7_1.fq.gz Raw_Reads/HV_007_06_forward.fq.gz
mv HV_007_7_USPD16099388-7_H2FJ5BBXX_L7_1.fq.gz Raw_Reads/HV_007_07_forward.fq.gz
mv HV_007_8_USPD16099388-8_H2FJ5BBXX_L7_1.fq.gz Raw_Reads/HV_007_08_forward.fq.gz
mv HV_007_9_USPD16099388-9_H2FJ5BBXX_L7_1.fq.gz Raw_Reads/HV_007_09_forward.fq.gz
mv HV_007_10_USPD16099388-10_H2FJ5BBXX_L7_1.fq.gz Raw_Reads/HV_007_10_forward.fq.gz
mv HV_007_11_USPD16099388-11_H2FJ5BBXX_L7_1.fq.gz Raw_Reads/HV_007_11_forward.fq.gz
mv HV_007_12_USPD16099388-12_H2FJ5BBXX_L7_1.fq.gz Raw_Reads/HV_007_12_forward.fq.gz
mv HV_007_13_USPD16099388-13_H2FJ5BBXX_L7_1.fq.gz Raw_Reads/HV_007_13_forward.fq.gz
mv HV_007_14_USPD16099388-14_H2FJ5BBXX_L7_1.fq.gz Raw_Reads/HV_007_14_forward.fq.gz
mv HV_007_15_USPD16099388-15_H2FJ5BBXX_L7_1.fq.gz Raw_Reads/HV_007_15_forward.fq.gz
mv HV_007_16_USPD16099388-16_H2FJ5BBXX_L7_1.fq.gz Raw_Reads/HV_007_16_forward.fq.gz
mv HV_007_18_USPD16099388-18_H2FJ5BBXX_L7_1.fq.gz Raw_Reads/HV_007_18_forward.fq.gz
mv HV_007_19_USPD16099388-19_H2FJ5BBXX_L7_1.fq.gz Raw_Reads/HV_007_19_forward.fq.gz
mv HV_007_20_USPD16099388-20_H2FJ5BBXX_L7_1.fq.gz Raw_Reads/HV_007_20_forward.fq.gz
mv HV_007_21_USPD16099388-21_H2FJ5BBXX_L7_1.fq.gz Raw_Reads/HV_007_21_forward.fq.gz
mv HV_007_22_USPD16099388-22_H2FJ5BBXX_L7_1.fq.gz Raw_Reads/HV_007_22_forward.fq.gz
mv HV_007_23_USPD16099388-23_H2FJ5BBXX_L7_1.fq.gz Raw_Reads/HV_007_23_forward.fq.gz
mv HV_007_25_USPD16099388-25_H2FJ5BBXX_L7_1.fq.gz Raw_Reads/HV_007_25_forward.fq.gz
mv HV_007_27_USPD16099388-27_H2FJ5BBXX_L7_1.fq.gz Raw_Reads/HV_007_27_forward.fq.gz
mv HV_008_1_USPD16099385-1_H2FJ5BBXX_L4_1.fq.gz Raw_Reads/HV_008_01_forward.fq.gz
mv HV_008_2_USPD16099385-2_H2FJ5BBXX_L4_1.fq.gz Raw_Reads/HV_008_02_forward.fq.gz
mv HV_008_3_USPD16099385-3_H2FJ5BBXX_L4_1.fq.gz Raw_Reads/HV_008_03_forward.fq.gz
mv HV_008_4_USPD16099385-4_H2FJ5BBXX_L4_1.fq.gz Raw_Reads/HV_008_04_forward.fq.gz
mv HV_008_5_USPD16099385-5_H2FJ5BBXX_L4_1.fq.gz Raw_Reads/HV_008_05_forward.fq.gz
mv HV_008_6_USPD16099385-6_H2FJ5BBXX_L4_1.fq.gz Raw_Reads/HV_008_06_forward.fq.gz
mv HV_008_7_USPD16099385-7_H2FJ5BBXX_L4_1.fq.gz Raw_Reads/HV_008_07_forward.fq.gz
mv HV_008_8_USPD16099385-8_H2FJ5BBXX_L4_1.fq.gz Raw_Reads/HV_008_08_forward.fq.gz
mv HV_008_9_USPD16099385-9_H2FJ5BBXX_L4_1.fq.gz Raw_Reads/HV_008_09_forward.fq.gz
mv HV_008_10_USPD16099385-10_H2FJ5BBXX_L4_1.fq.gz Raw_Reads/HV_008_10_forward.fq.gz
mv HV_008_11_USPD16099385-11_H2FJ5BBXX_L4_1.fq.gz Raw_Reads/HV_008_11_forward.fq.gz
mv HV_008_12_USPD16099385-12_H2FJ5BBXX_L4_1.fq.gz Raw_Reads/HV_008_12_forward.fq.gz
mv HV_008_13_USPD16099385-13_H2FJ5BBXX_L4_1.fq.gz Raw_Reads/HV_008_13_forward.fq.gz
mv HV_008_14_USPD16099385-14_H2FJ5BBXX_L4_1.fq.gz Raw_Reads/HV_008_14_forward.fq.gz
mv HV_008_15_USPD16099385-15_H2FJ5BBXX_L4_1.fq.gz Raw_Reads/HV_008_15_forward.fq.gz
mv HV_008_16_USPD16099385-16_H2FJ5BBXX_L4_1.fq.gz Raw_Reads/HV_008_16_forward.fq.gz
mv HV_008_18_USPD16099385-18_H2FJ5BBXX_L4_1.fq.gz Raw_Reads/HV_008_18_forward.fq.gz
mv HV_008_19_USPD16099385-19_H2FJ5BBXX_L4_1.fq.gz Raw_Reads/HV_008_19_forward.fq.gz
mv HV_008_20_USPD16099385-20_H2FJ5BBXX_L4_1.fq.gz Raw_Reads/HV_008_20_forward.fq.gz
mv HV_008_21_USPD16099385-21_H2FJ5BBXX_L4_1.fq.gz Raw_Reads/HV_008_21_forward.fq.gz
mv HV_008_22_USPD16099385-22_H2FJ5BBXX_L4_1.fq.gz Raw_Reads/HV_008_22_forward.fq.gz
mv HV_008_23_USPD16099385-23_H2FJ5BBXX_L4_1.fq.gz Raw_Reads/HV_008_23_forward.fq.gz
mv HV_008_25_USPD16099385-25_H2FJ5BBXX_L4_1.fq.gz Raw_Reads/HV_008_25_forward.fq.gz
mv HV_008_27_USPD16099385-27_H2FJ5BBXX_L4_1.fq.gz Raw_Reads/HV_008_27_forward.fq.gz
mv HV_009_1_USPD16099386-1_H2FJ5BBXX_L5_1.fq.gz Raw_Reads/HV_009_01_forward.fq.gz
mv HV_009_3_USPD16099386-3_H2FJ5BBXX_L5_1.fq.gz Raw_Reads/HV_009_03_forward.fq.gz
mv HV_009_4_USPD16099386-4_H2FJ5BBXX_L5_1.fq.gz Raw_Reads/HV_009_04_forward.fq.gz
mv HV_009_5_USPD16099386-5_H2FJ5BBXX_L5_1.fq.gz Raw_Reads/HV_009_05_forward.fq.gz
mv HV_009_6_USPD16099386-6_H2FJ5BBXX_L5_1.fq.gz Raw_Reads/HV_009_06_forward.fq.gz
mv HV_009_7_USPD16099386-7_H2FJ5BBXX_L5_1.fq.gz Raw_Reads/HV_009_07_forward.fq.gz
mv HV_009_8_USPD16099386-8_H2FJ5BBXX_L5_1.fq.gz Raw_Reads/HV_009_08_forward.fq.gz
mv HV_009_9_USPD16099386-9_H2FJ5BBXX_L5_1.fq.gz Raw_Reads/HV_009_09_forward.fq.gz
mv HV_009_10_USPD16099386-10_H2FJ5BBXX_L5_1.fq.gz Raw_Reads/HV_009_10_forward.fq.gz
mv HV_009_11_USPD16099386-11_H2FJ5BBXX_L5_1.fq.gz Raw_Reads/HV_009_11_forward.fq.gz
mv HV_009_12_USPD16099386-12_H2FJ5BBXX_L5_1.fq.gz Raw_Reads/HV_009_12_forward.fq.gz
mv HV_009_13_USPD16099386-13_H2FJ5BBXX_L5_1.fq.gz Raw_Reads/HV_009_13_forward.fq.gz
mv HV_009_14_USPD16099386-14_H2FJ5BBXX_L5_1.fq.gz Raw_Reads/HV_009_14_forward.fq.gz
mv HV_009_15_USPD16099386-15_H2FJ5BBXX_L5_1.fq.gz Raw_Reads/HV_009_15_forward.fq.gz
mv HV_009_16_USPD16099386-16_H2FJ5BBXX_L5_1.fq.gz Raw_Reads/HV_009_16_forward.fq.gz
mv HV_009_18_USPD16099386-18_H2FJ5BBXX_L5_1.fq.gz Raw_Reads/HV_009_18_forward.fq.gz
mv HV_009_19_USPD16099386-19_H2FJ5BBXX_L5_1.fq.gz Raw_Reads/HV_009_19_forward.fq.gz
mv HV_009_20_USPD16099386-20_H2FJ5BBXX_L5_1.fq.gz Raw_Reads/HV_009_20_forward.fq.gz
mv HV_009_21_USPD16099386-21_H2FJ5BBXX_L5_1.fq.gz Raw_Reads/HV_009_21_forward.fq.gz
mv HV_009_22_USPD16099386-22_H2FJ5BBXX_L5_1.fq.gz Raw_Reads/HV_009_22_forward.fq.gz
mv HV_009_23_USPD16099386-23_H2FJ5BBXX_L5_1.fq.gz Raw_Reads/HV_009_23_forward.fq.gz
mv HV_009_25_USPD16099386-25_H2FJ5BBXX_L5_1.fq.gz Raw_Reads/HV_009_25_forward.fq.gz
mv HV_009_27_USPD16099386-27_H2FJ5BBXX_L5_1.fq.gz Raw_Reads/HV_009_27_forward.fq.gz
mv HV_010_1_USPD16099384-1_H2FJ5BBXX_L3_1.fq.gz Raw_Reads/HV_010_01_forward.fq.gz
mv HV_010_2_USPD16099384-2_H2FJ5BBXX_L3_1.fq.gz Raw_Reads/HV_010_02_forward.fq.gz
mv HV_010_3_USPD16099384-3_H2FJ5BBXX_L3_1.fq.gz Raw_Reads/HV_010_03_forward.fq.gz
mv HV_010_4_USPD16099384-4_H2FJ5BBXX_L3_1.fq.gz Raw_Reads/HV_010_04_forward.fq.gz
mv HV_010_5_USPD16099384-5_H2FJ5BBXX_L3_1.fq.gz Raw_Reads/HV_010_05_forward.fq.gz
mv HV_010_6_USPD16099384-6_H2FJ5BBXX_L3_1.fq.gz Raw_Reads/HV_010_06_forward.fq.gz
mv HV_010_7_USPD16099384-7_H2FJ5BBXX_L3_1.fq.gz Raw_Reads/HV_010_07_forward.fq.gz
mv HV_010_8_USPD16099384-8_H2FJ5BBXX_L3_1.fq.gz Raw_Reads/HV_010_08_forward.fq.gz
mv HV_010_9_USPD16099384-9_H2FJ5BBXX_L3_1.fq.gz Raw_Reads/HV_010_09_forward.fq.gz
mv HV_010_10_USPD16099384-10_H2FJ5BBXX_L3_1.fq.gz Raw_Reads/HV_010_10_forward.fq.gz
mv HV_010_11_USPD16099384-11_H2FJ5BBXX_L3_1.fq.gz Raw_Reads/HV_010_11_forward.fq.gz
mv HV_010_12_USPD16099384-12_H2FJ5BBXX_L3_1.fq.gz Raw_Reads/HV_010_12_forward.fq.gz
mv HV_010_13_USPD16099384-13_H2FJ5BBXX_L3_1.fq.gz Raw_Reads/HV_010_13_forward.fq.gz
mv HV_010_14_USPD16099384-14_H2FJ5BBXX_L3_1.fq.gz Raw_Reads/HV_010_14_forward.fq.gz
mv HV_010_15_USPD16099384-15_H2FJ5BBXX_L3_1.fq.gz Raw_Reads/HV_010_15_forward.fq.gz
mv HV_010_16_USPD16099384-16_H2FJ5BBXX_L3_1.fq.gz Raw_Reads/HV_010_16_forward.fq.gz
mv HV_010_18_USPD16099384-18_H2FJ5BBXX_L3_1.fq.gz Raw_Reads/HV_010_18_forward.fq.gz
mv HV_010_19_USPD16099384-19_H2FJ5BBXX_L3_1.fq.gz Raw_Reads/HV_010_19_forward.fq.gz
mv HV_010_20_USPD16099384-20_H2FJ5BBXX_L3_1.fq.gz Raw_Reads/HV_010_20_forward.fq.gz
mv HV_010_21_USPD16099384-21_H2FJ5BBXX_L3_1.fq.gz Raw_Reads/HV_010_21_forward.fq.gz
mv HV_010_22_USPD16099384-22_H2FJ5BBXX_L3_1.fq.gz Raw_Reads/HV_010_22_forward.fq.gz
mv HV_010_25_USPD16099384-25_H2FJ5BBXX_L3_1.fq.gz Raw_Reads/HV_010_25_forward.fq.gz
mv HV_010_27_USPD16099384-27_H2FJ5BBXX_L3_1.fq.gz Raw_Reads/HV_010_27_forward.fq.gz
mv HV_011_1_USPD16099389-1_H2FJ5BBXX_L8_1.fq.gz Raw_Reads/HV_011_01_forward.fq.gz
mv HV_011_2_USPD16099389-2_H2FJ5BBXX_L8_1.fq.gz Raw_Reads/HV_011_02_forward.fq.gz
mv HV_011_3_USPD16099389-3_H2FJ5BBXX_L8_1.fq.gz Raw_Reads/HV_011_03_forward.fq.gz
mv HV_011_4_USPD16099389-4_H2FJ5BBXX_L8_1.fq.gz Raw_Reads/HV_011_04_forward.fq.gz
mv HV_011_5_USPD16099389-5_H2FJ5BBXX_L8_1.fq.gz Raw_Reads/HV_011_05_forward.fq.gz
mv HV_011_6_USPD16099389-6_H2FJ5BBXX_L8_1.fq.gz Raw_Reads/HV_011_06_forward.fq.gz
mv HV_011_7_USPD16099389-7_H2FJ5BBXX_L8_1.fq.gz Raw_Reads/HV_011_07_forward.fq.gz
mv HV_011_8_USPD16099389-8_H2FJ5BBXX_L8_1.fq.gz Raw_Reads/HV_011_08_forward.fq.gz
mv HV_011_9_USPD16099389-9_H2FJ5BBXX_L8_1.fq.gz Raw_Reads/HV_011_09_forward.fq.gz
mv HV_011_10_USPD16099389-10_H2FJ5BBXX_L8_1.fq.gz Raw_Reads/HV_011_10_forward.fq.gz
mv HV_011_11_USPD16099389-11_H2FJ5BBXX_L8_1.fq.gz Raw_Reads/HV_011_11_forward.fq.gz
mv HV_011_12_USPD16099389-12_H2FJ5BBXX_L8_1.fq.gz Raw_Reads/HV_011_12_forward.fq.gz
mv HV_011_13_USPD16099389-13_H2FJ5BBXX_L8_1.fq.gz Raw_Reads/HV_011_13_forward.fq.gz
mv HV_011_14_USPD16099389-14_H2FJ5BBXX_L8_1.fq.gz Raw_Reads/HV_011_14_forward.fq.gz
mv HV_011_15_USPD16099389-15_H2FJ5BBXX_L8_1.fq.gz Raw_Reads/HV_011_15_forward.fq.gz
mv HV_011_16_USPD16099389-16_H2FJ5BBXX_L8_1.fq.gz Raw_Reads/HV_011_16_forward.fq.gz
mv HV_011_18_USPD16099389-18_H2FJ5BBXX_L8_1.fq.gz Raw_Reads/HV_011_18_forward.fq.gz
mv HV_011_19_USPD16099389-19_H2FJ5BBXX_L8_1.fq.gz Raw_Reads/HV_011_19_forward.fq.gz
mv HV_011_20_USPD16099389-20_H2FJ5BBXX_L8_1.fq.gz Raw_Reads/HV_011_20_forward.fq.gz
mv HV_011_21_USPD16099389-21_H2FJ5BBXX_L8_1.fq.gz Raw_Reads/HV_011_21_forward.fq.gz
mv HV_011_22_USPD16099389-22_H2FJ5BBXX_L8_1.fq.gz Raw_Reads/HV_011_22_forward.fq.gz
mv HV_011_23_USPD16099389-23_H2FJ5BBXX_L8_1.fq.gz Raw_Reads/HV_011_23_forward.fq.gz
mv HV_011_25_USPD16099389-25_H2FJ5BBXX_L8_1.fq.gz Raw_Reads/HV_011_25_forward.fq.gz
mv HV_011_27_USPD16099389-27_H2FJ5BBXX_L8_1.fq.gz Raw_Reads/HV_011_27_forward.fq.gz
mv HV012_1_USPD16102686-1_H333VBBXX_L2_1.fq.gz Raw_Reads/HV_012_01_forward.fq.gz
mv HV012_2_USPD16102686-2_H333VBBXX_L2_1.fq.gz Raw_Reads/HV_012_02_forward.fq.gz
mv HV012_3_USPD16102686-3_H333VBBXX_L2_1.fq.gz Raw_Reads/HV_012_03_forward.fq.gz
mv HV012_4_USPD16102686-4_H333VBBXX_L2_1.fq.gz Raw_Reads/HV_012_04_forward.fq.gz
mv HV012_5_USPD16102686-5_H333VBBXX_L2_1.fq.gz Raw_Reads/HV_012_05_forward.fq.gz
mv HV012_6_USPD16102686-6_H333VBBXX_L2_1.fq.gz Raw_Reads/HV_012_06_forward.fq.gz
mv HV012_7_USPD16102686-7_H333VBBXX_L2_1.fq.gz Raw_Reads/HV_012_07_forward.fq.gz
mv HV012_8_USPD16102686-8_H333VBBXX_L2_1.fq.gz Raw_Reads/HV_012_08_forward.fq.gz
mv HV012_9_USPD16102686-9_H333VBBXX_L2_1.fq.gz Raw_Reads/HV_012_09_forward.fq.gz
mv HV012_10_USPD16102686-10_H333VBBXX_L2_1.fq.gz Raw_Reads/HV_012_10_forward.fq.gz
mv HV012_11_USPD16102686-11_H333VBBXX_L2_1.fq.gz Raw_Reads/HV_012_11_forward.fq.gz
mv HV012_12_USPD16102686-12_H333VBBXX_L2_1.fq.gz Raw_Reads/HV_012_12_forward.fq.gz
mv HV012_13_USPD16102686-13_H333VBBXX_L2_1.fq.gz Raw_Reads/HV_012_13_forward.fq.gz
mv HV012_14_USPD16102686-14_H333VBBXX_L2_1.fq.gz Raw_Reads/HV_012_14_forward.fq.gz
mv HV012_15_USPD16102686-15_H333VBBXX_L2_1.fq.gz Raw_Reads/HV_012_15_forward.fq.gz
mv HV012_16_USPD16102686-16_H333VBBXX_L2_1.fq.gz Raw_Reads/HV_012_16_forward.fq.gz
mv HV012_18_USPD16102686-18_H333VBBXX_L2_1.fq.gz Raw_Reads/HV_012_18_forward.fq.gz
mv HV012_19_USPD16102686-19_H333VBBXX_L2_1.fq.gz Raw_Reads/HV_012_19_forward.fq.gz
mv HV012_20_USPD16102686-20_H333VBBXX_L2_1.fq.gz Raw_Reads/HV_012_20_forward.fq.gz
mv HV012_21_USPD16102686-21_H333VBBXX_L2_1.fq.gz Raw_Reads/HV_012_21_forward.fq.gz
mv HV012_22_USPD16102686-22_H333VBBXX_L2_1.fq.gz Raw_Reads/HV_012_22_forward.fq.gz
mv HV012_23_USPD16102686-23_H333VBBXX_L2_1.fq.gz Raw_Reads/HV_012_23_forward.fq.gz
mv HV012_25_USPD16102686-25_H333VBBXX_L2_1.fq.gz Raw_Reads/HV_012_25_forward.fq.gz
mv HV012_27_USPD16102686-27_H333VBBXX_L2_1.fq.gz Raw_Reads/HV_012_27_forward.fq.gz
mv HV013_1_USPD16102680-1_H3FVNBBXX_L1_1.fq.gz Raw_Reads/HV_013_01_forward.fq.gz
mv HV013_2_USPD16102680-2_H3FVNBBXX_L1_1.fq.gz Raw_Reads/HV_013_02_forward.fq.gz
mv HV013_3_USPD16102680-3_H3FVNBBXX_L1_1.fq.gz Raw_Reads/HV_013_03_forward.fq.gz
mv HV013_4_USPD16102680-4_H3FVNBBXX_L1_1.fq.gz Raw_Reads/HV_013_04_forward.fq.gz
mv HV013_5_USPD16102680-5_H3FVNBBXX_L1_1.fq.gz Raw_Reads/HV_013_05_forward.fq.gz
mv HV013_7_USPD16102680-7_H3FVNBBXX_L1_1.fq.gz Raw_Reads/HV_013_07_forward.fq.gz
mv HV013_8_USPD16102680-8_H3FVNBBXX_L1_1.fq.gz Raw_Reads/HV_013_08_forward.fq.gz
mv HV013_9_USPD16102680-9_H3FVNBBXX_L1_1.fq.gz Raw_Reads/HV_013_09_forward.fq.gz
mv HV013_10_USPD16102680-10_H3FVNBBXX_L1_1.fq.gz Raw_Reads/HV_013_10_forward.fq.gz
mv HV013_11_USPD16102680-11_H3FVNBBXX_L1_1.fq.gz Raw_Reads/HV_013_11_forward.fq.gz
mv HV013_12_USPD16102680-12_H3FVNBBXX_L1_1.fq.gz Raw_Reads/HV_013_12_forward.fq.gz
mv HV013_13_USPD16102680-13_H3FVNBBXX_L1_1.fq.gz Raw_Reads/HV_013_13_forward.fq.gz
mv HV013_14_USPD16102680-14_H3FVNBBXX_L1_1.fq.gz Raw_Reads/HV_013_14_forward.fq.gz
mv HV013_15_USPD16102680-15_H3FVNBBXX_L1_1.fq.gz Raw_Reads/HV_013_15_forward.fq.gz
mv HV013_16_USPD16102680-16_H3FVNBBXX_L1_1.fq.gz Raw_Reads/HV_013_16_forward.fq.gz
mv HV013_18_USPD16102680-18_H3FVNBBXX_L1_1.fq.gz Raw_Reads/HV_013_18_forward.fq.gz
mv HV013_19_USPD16102680-19_H3FVNBBXX_L1_1.fq.gz Raw_Reads/HV_013_19_forward.fq.gz
mv HV013_20_USPD16102680-20_H3FVNBBXX_L1_1.fq.gz Raw_Reads/HV_013_20_forward.fq.gz
mv HV013_21_USPD16102680-21_H3FVNBBXX_L1_1.fq.gz Raw_Reads/HV_013_21_forward.fq.gz
mv HV013_22_USPD16102680-22_H3FVNBBXX_L1_1.fq.gz Raw_Reads/HV_013_22_forward.fq.gz
mv HV013_23_USPD16102680-23_H3FVNBBXX_L1_1.fq.gz Raw_Reads/HV_013_23_forward.fq.gz
mv HV013_25_USPD16102680-25_H3FVNBBXX_L1_1.fq.gz Raw_Reads/HV_013_25_forward.fq.gz
mv HV013_27_USPD16102680-27_H3FVNBBXX_L1_1.fq.gz Raw_Reads/HV_013_27_forward.fq.gz
mv HV014_1_USPD16102681-1_H3FVNBBXX_L2_1.fq.gz Raw_Reads/HV_014_01_forward.fq.gz
mv HV014_2_USPD16102681-2_H3FVNBBXX_L2_1.fq.gz Raw_Reads/HV_014_02_forward.fq.gz
mv HV014_3_USPD16102681-3_H3FVNBBXX_L2_1.fq.gz Raw_Reads/HV_014_03_forward.fq.gz
mv HV014_4_USPD16102681-4_H3FVNBBXX_L2_1.fq.gz Raw_Reads/HV_014_04_forward.fq.gz
mv HV014_5_USPD16102681-5_H3FVNBBXX_L2_1.fq.gz Raw_Reads/HV_014_05_forward.fq.gz
mv HV014_6_USPD16102681-6_H3FVNBBXX_L2_1.fq.gz Raw_Reads/HV_014_06_forward.fq.gz
mv HV014_7_USPD16102681-7_H3FVNBBXX_L2_1.fq.gz Raw_Reads/HV_014_07_forward.fq.gz
mv HV014_8_USPD16102681-8_H3FVNBBXX_L2_1.fq.gz Raw_Reads/HV_014_08_forward.fq.gz
mv HV014_9_USPD16102681-9_H3FVNBBXX_L2_1.fq.gz Raw_Reads/HV_014_09_forward.fq.gz
mv HV014_10_USPD16102681-10_H3FVNBBXX_L2_1.fq.gz Raw_Reads/HV_014_10_forward.fq.gz
mv HV014_11_USPD16102681-11_H3FVNBBXX_L2_1.fq.gz Raw_Reads/HV_014_11_forward.fq.gz
mv HV014_12_USPD16102681-12_H3FVNBBXX_L2_1.fq.gz Raw_Reads/HV_014_12_forward.fq.gz
mv HV014_13_USPD16102681-13_H3FVNBBXX_L2_1.fq.gz Raw_Reads/HV_014_13_forward.fq.gz
mv HV014_14_USPD16102681-14_H3FVNBBXX_L2_1.fq.gz Raw_Reads/HV_014_14_forward.fq.gz
mv HV014_15_USPD16102681-15_H3FVNBBXX_L2_1.fq.gz Raw_Reads/HV_014_15_forward.fq.gz
mv HV014_16_USPD16102681-16_H3FVNBBXX_L2_1.fq.gz Raw_Reads/HV_014_16_forward.fq.gz
mv HV014_18_USPD16102681-18_H3FVNBBXX_L2_1.fq.gz Raw_Reads/HV_014_18_forward.fq.gz
mv HV014_19_USPD16102681-19_H3FVNBBXX_L2_1.fq.gz Raw_Reads/HV_014_19_forward.fq.gz
mv HV014_20_USPD16102681-20_H3FVNBBXX_L2_1.fq.gz Raw_Reads/HV_014_20_forward.fq.gz
mv HV014_21_USPD16102681-21_H3FVNBBXX_L2_1.fq.gz Raw_Reads/HV_014_21_forward.fq.gz
mv HV014_22_USPD16102681-22_H3FVNBBXX_L2_1.fq.gz Raw_Reads/HV_014_22_forward.fq.gz
mv HV014_23_USPD16102681-23_H3FVNBBXX_L2_1.fq.gz Raw_Reads/HV_014_23_forward.fq.gz
mv HV014_25_USPD16102681-25_H3FVNBBXX_L2_1.fq.gz Raw_Reads/HV_014_25_forward.fq.gz
mv HV014_27_USPD16102681-27_H3FVNBBXX_L2_1.fq.gz Raw_Reads/HV_014_27_forward.fq.gz
mv HV015_1_USPD16102683-1_H3FVNBBXX_L4_1.fq.gz Raw_Reads/HV_015_01_forward.fq.gz
mv HV015_2_USPD16102683-2_H3FVNBBXX_L4_1.fq.gz Raw_Reads/HV_015_02_forward.fq.gz
mv HV015_3_USPD16102683-3_H3FVNBBXX_L4_1.fq.gz Raw_Reads/HV_015_03_forward.fq.gz
mv HV015_4_USPD16102683-4_H3FVNBBXX_L4_1.fq.gz Raw_Reads/HV_015_04_forward.fq.gz
mv HV015_5_USPD16102683-5_H3FVNBBXX_L4_1.fq.gz Raw_Reads/HV_015_05_forward.fq.gz
mv HV015_6_USPD16102683-6_H3FVNBBXX_L4_1.fq.gz Raw_Reads/HV_015_06_forward.fq.gz
mv HV015_7_USPD16102683-7_H3FVNBBXX_L4_1.fq.gz Raw_Reads/HV_015_07_forward.fq.gz
mv HV015_8_USPD16102683-8_H3FVNBBXX_L4_1.fq.gz Raw_Reads/HV_015_08_forward.fq.gz
mv HV015_9_USPD16102683-9_H3FVNBBXX_L4_1.fq.gz Raw_Reads/HV_015_09_forward.fq.gz
mv HV015_10_USPD16102683-10_H3FVNBBXX_L4_1.fq.gz Raw_Reads/HV_015_10_forward.fq.gz
mv HV015_11_USPD16102683-11_H3FVNBBXX_L4_1.fq.gz Raw_Reads/HV_015_11_forward.fq.gz
mv HV015_12_USPD16102683-12_H3FVNBBXX_L4_1.fq.gz Raw_Reads/HV_015_12_forward.fq.gz
mv HV015_13_USPD16102683-13_H3FVNBBXX_L4_1.fq.gz Raw_Reads/HV_015_13_forward.fq.gz
mv HV015_14_USPD16102683-14_H3FVNBBXX_L4_1.fq.gz Raw_Reads/HV_015_14_forward.fq.gz
mv HV015_15_USPD16102683-15_H3FVNBBXX_L4_1.fq.gz Raw_Reads/HV_015_15_forward.fq.gz
mv HV015_16_USPD16102683-16_H3FVNBBXX_L4_1.fq.gz Raw_Reads/HV_015_16_forward.fq.gz
mv HV015_18_USPD16102683-18_H3FVNBBXX_L4_1.fq.gz Raw_Reads/HV_015_18_forward.fq.gz
mv HV015_19_USPD16102683-19_H3FVNBBXX_L4_1.fq.gz Raw_Reads/HV_015_19_forward.fq.gz
mv HV015_20_USPD16102683-20_H3FVNBBXX_L4_1.fq.gz Raw_Reads/HV_015_20_forward.fq.gz
mv HV015_21_USPD16102683-21_H3FVNBBXX_L4_1.fq.gz Raw_Reads/HV_015_21_forward.fq.gz
mv HV015_23_USPD16102683-23_H3FVNBBXX_L4_1.fq.gz Raw_Reads/HV_015_23_forward.fq.gz
mv HV015_25_USPD16102683-25_H3FVNBBXX_L4_1.fq.gz Raw_Reads/HV_015_25_forward.fq.gz
mv HV015_27_USPD16102683-27_H3FVNBBXX_L4_1.fq.gz Raw_Reads/HV_015_27_forward.fq.gz
mv HV016_1_USPD16102682-1_H3FVNBBXX_L3_1.fq.gz Raw_Reads/HV_016_01_forward.fq.gz
mv HV016_2_USPD16102682-2_H3FVNBBXX_L3_1.fq.gz Raw_Reads/HV_016_02_forward.fq.gz
mv HV016_3_USPD16102682-3_H3FVNBBXX_L3_1.fq.gz Raw_Reads/HV_016_03_forward.fq.gz
mv HV016_4_USPD16102682-4_H3FVNBBXX_L3_1.fq.gz Raw_Reads/HV_016_04_forward.fq.gz
mv HV016_5_USPD16102682-5_H3FVNBBXX_L3_1.fq.gz Raw_Reads/HV_016_05_forward.fq.gz
mv HV016_6_USPD16102682-6_H3FVNBBXX_L3_1.fq.gz Raw_Reads/HV_016_06_forward.fq.gz
mv HV016_7_USPD16102682-7_H3FVNBBXX_L3_1.fq.gz Raw_Reads/HV_016_07_forward.fq.gz
mv HV016_8_USPD16102682-8_H3FVNBBXX_L3_1.fq.gz Raw_Reads/HV_016_08_forward.fq.gz
mv HV016_9_USPD16102682-9_H3FVNBBXX_L3_1.fq.gz Raw_Reads/HV_016_09_forward.fq.gz
mv HV016_10_USPD16102682-10_H3FVNBBXX_L3_1.fq.gz Raw_Reads/HV_016_10_forward.fq.gz
mv HV016_11_USPD16102682-11_H3FVNBBXX_L3_1.fq.gz Raw_Reads/HV_016_11_forward.fq.gz
mv HV016_12_USPD16102682-12_H3FVNBBXX_L3_1.fq.gz Raw_Reads/HV_016_12_forward.fq.gz
mv HV016_13_USPD16102682-13_H3FVNBBXX_L3_1.fq.gz Raw_Reads/HV_016_13_forward.fq.gz
mv HV016_14_USPD16102682-14_H3FVNBBXX_L3_1.fq.gz Raw_Reads/HV_016_14_forward.fq.gz
mv HV016_15_USPD16102682-15_H3FVNBBXX_L3_1.fq.gz Raw_Reads/HV_016_15_forward.fq.gz
mv HV016_16_USPD16102682-16_H3FVNBBXX_L3_1.fq.gz Raw_Reads/HV_016_16_forward.fq.gz
mv HV016_18_USPD16102682-18_H3FVNBBXX_L3_1.fq.gz Raw_Reads/HV_016_18_forward.fq.gz
mv HV016_19_USPD16102682-19_H3FVNBBXX_L3_1.fq.gz Raw_Reads/HV_016_19_forward.fq.gz
mv HV016_20_USPD16102682-20_H3FVNBBXX_L3_1.fq.gz Raw_Reads/HV_016_20_forward.fq.gz
mv HV016_21_USPD16102682-21_H3FVNBBXX_L3_1.fq.gz Raw_Reads/HV_016_21_forward.fq.gz
mv HV016_22_USPD16102682-22_H3FVNBBXX_L3_1.fq.gz Raw_Reads/HV_016_22_forward.fq.gz
mv HV016_23_USPD16102682-23_H3FVNBBXX_L3_1.fq.gz Raw_Reads/HV_016_23_forward.fq.gz
mv HV016_25_USPD16102682-25_H3FVNBBXX_L3_1.fq.gz Raw_Reads/HV_016_25_forward.fq.gz
mv HV016_27_USPD16102682-27_H3FVNBBXX_L3_1.fq.gz Raw_Reads/HV_016_27_forward.fq.gz
mv HV017_1_USPD16102684-1_H3FVNBBXX_L5_1.fq.gz Raw_Reads/HV_017_01_forward.fq.gz
mv HV017_2_USPD16102684-2_H3FVNBBXX_L5_1.fq.gz Raw_Reads/HV_017_02_forward.fq.gz
mv HV017_3_USPD16102684-3_H3FVNBBXX_L5_1.fq.gz Raw_Reads/HV_017_03_forward.fq.gz
mv HV017_4_USPD16102684-4_H3FVNBBXX_L5_1.fq.gz Raw_Reads/HV_017_04_forward.fq.gz
mv HV017_5_USPD16102684-5_H3FVNBBXX_L5_1.fq.gz Raw_Reads/HV_017_05_forward.fq.gz
mv HV017_6_USPD16102684-6_H3FVNBBXX_L5_1.fq.gz Raw_Reads/HV_017_06_forward.fq.gz
mv HV017_7_USPD16102684-7_H3FVNBBXX_L5_1.fq.gz Raw_Reads/HV_017_07_forward.fq.gz
mv HV017_8_USPD16102684-8_H3FVNBBXX_L5_1.fq.gz Raw_Reads/HV_017_08_forward.fq.gz
mv HV017_9_USPD16102684-9_H3FVNBBXX_L5_1.fq.gz Raw_Reads/HV_017_09_forward.fq.gz
mv HV017_10_USPD16102684-10_H3FVNBBXX_L5_1.fq.gz Raw_Reads/HV_017_10_forward.fq.gz
mv HV017_11_USPD16102684-11_H3FVNBBXX_L5_1.fq.gz Raw_Reads/HV_017_11_forward.fq.gz
mv HV017_12_USPD16102684-12_H3FVNBBXX_L5_1.fq.gz Raw_Reads/HV_017_12_forward.fq.gz
mv HV017_13_USPD16102684-13_H3FVNBBXX_L5_1.fq.gz Raw_Reads/HV_017_13_forward.fq.gz
mv HV017_14_USPD16102684-14_H3FVNBBXX_L5_1.fq.gz Raw_Reads/HV_017_14_forward.fq.gz
mv HV017_15_USPD16102684-15_H3FVNBBXX_L5_1.fq.gz Raw_Reads/HV_017_15_forward.fq.gz
mv HV017_16_USPD16102684-16_H3FVNBBXX_L5_1.fq.gz Raw_Reads/HV_017_16_forward.fq.gz
mv HV017_18_USPD16102684-18_H3FVNBBXX_L5_1.fq.gz Raw_Reads/HV_017_18_forward.fq.gz
mv HV017_19_USPD16102684-19_H3FVNBBXX_L5_1.fq.gz Raw_Reads/HV_017_19_forward.fq.gz
mv HV017_20_USPD16102684-20_H3FVNBBXX_L5_1.fq.gz Raw_Reads/HV_017_20_forward.fq.gz
mv HV017_21_USPD16102684-21_H3FVNBBXX_L5_1.fq.gz Raw_Reads/HV_017_21_forward.fq.gz
mv HV017_22_USPD16102684-22_H3FVNBBXX_L5_1.fq.gz Raw_Reads/HV_017_22_forward.fq.gz
mv HV017_23_USPD16102684-23_H3FVNBBXX_L5_1.fq.gz Raw_Reads/HV_017_23_forward.fq.gz
mv HV017_25_USPD16102684-25_H3FVNBBXX_L5_1.fq.gz Raw_Reads/HV_017_25_forward.fq.gz
mv HV017_27_USPD16102684-27_H3FVNBBXX_L5_1.fq.gz Raw_Reads/HV_017_27_forward.fq.gz
mv HV018_01_1.fq.gz Raw_Reads/HV_018_01_forward.fq.gz
mv HV018_02_1.fq.gz Raw_Reads/HV_018_02_forward.fq.gz
mv HV018_03_1.fq.gz Raw_Reads/HV_018_03_forward.fq.gz
mv HV018_04_1.fq.gz Raw_Reads/HV_018_04_forward.fq.gz
mv HV018_05_1.fq.gz Raw_Reads/HV_018_05_forward.fq.gz
mv HV018_06_1.fq.gz Raw_Reads/HV_018_06_forward.fq.gz
mv HV018_07_1.fq.gz Raw_Reads/HV_018_07_forward.fq.gz
mv HV018_08_1.fq.gz Raw_Reads/HV_018_08_forward.fq.gz
mv HV018_09_1.fq.gz Raw_Reads/HV_018_09_forward.fq.gz
mv HV018_10_1.fq.gz Raw_Reads/HV_018_10_forward.fq.gz
mv HV018_11_1.fq.gz Raw_Reads/HV_018_11_forward.fq.gz
mv HV018_12_1.fq.gz Raw_Reads/HV_018_12_forward.fq.gz
mv HV018_13_1.fq.gz Raw_Reads/HV_018_13_forward.fq.gz
mv HV018_14_1.fq.gz Raw_Reads/HV_018_14_forward.fq.gz
mv HV018_15_1.fq.gz Raw_Reads/HV_018_15_forward.fq.gz
mv HV018_16_1.fq.gz Raw_Reads/HV_018_16_forward.fq.gz
mv HV018_18_1.fq.gz Raw_Reads/HV_018_18_forward.fq.gz
mv HV018_19_1.fq.gz Raw_Reads/HV_018_19_forward.fq.gz
mv HV018_20_1.fq.gz Raw_Reads/HV_018_20_forward.fq.gz
mv HV018_21_1.fq.gz Raw_Reads/HV_018_21_forward.fq.gz
mv HV018_22_1.fq.gz Raw_Reads/HV_018_22_forward.fq.gz
mv HV018_23_1.fq.gz Raw_Reads/HV_018_23_forward.fq.gz
mv HV018_25_1.fq.gz Raw_Reads/HV_018_25_forward.fq.gz
mv HV018_27_1.fq.gz Raw_Reads/HV_018_27_forward.fq.gz
mv HV019_01_1.fq.gz Raw_Reads/HV_019_01_forward.fq.gz
mv HV019_02_1.fq.gz Raw_Reads/HV_019_02_forward.fq.gz
mv HV019_03_1.fq.gz Raw_Reads/HV_019_03_forward.fq.gz
mv HV019_04_1.fq.gz Raw_Reads/HV_019_04_forward.fq.gz
mv HV019_05_1.fq.gz Raw_Reads/HV_019_05_forward.fq.gz
mv HV019_06_1.fq.gz Raw_Reads/HV_019_06_forward.fq.gz
mv HV019_07_1.fq.gz Raw_Reads/HV_019_07_forward.fq.gz
mv HV019_08_1.fq.gz Raw_Reads/HV_019_08_forward.fq.gz
mv HV019_09_1.fq.gz Raw_Reads/HV_019_09_forward.fq.gz
mv HV019_10_1.fq.gz Raw_Reads/HV_019_10_forward.fq.gz
mv HV019_11_1.fq.gz Raw_Reads/HV_019_11_forward.fq.gz
mv HV019_12_1.fq.gz Raw_Reads/HV_019_12_forward.fq.gz
mv HV019_13_1.fq.gz Raw_Reads/HV_019_13_forward.fq.gz
mv HV019_14_1.fq.gz Raw_Reads/HV_019_14_forward.fq.gz
mv HV019_15_1.fq.gz Raw_Reads/HV_019_15_forward.fq.gz
mv HV019_16_1.fq.gz Raw_Reads/HV_019_16_forward.fq.gz
mv HV019_18_1.fq.gz Raw_Reads/HV_019_18_forward.fq.gz
mv HV019_19_1.fq.gz Raw_Reads/HV_019_19_forward.fq.gz
mv HV019_20_1.fq.gz Raw_Reads/HV_019_20_forward.fq.gz
mv HV019_21_1.fq.gz Raw_Reads/HV_019_21_forward.fq.gz
mv HV019_22_1.fq.gz Raw_Reads/HV_019_22_forward.fq.gz
mv HV019_23_1.fq.gz Raw_Reads/HV_019_23_forward.fq.gz
mv HV019_25_1.fq.gz Raw_Reads/HV_019_25_forward.fq.gz
mv HV019_27_1.fq.gz Raw_Reads/HV_019_27_forward.fq.gz
mv HV020_01_1.fq.gz Raw_Reads/HV_020_01_forward.fq.gz
mv HV020_02_1.fq.gz Raw_Reads/HV_020_02_forward.fq.gz
mv HV020_03_1.fq.gz Raw_Reads/HV_020_03_forward.fq.gz
mv HV020_04_1.fq.gz Raw_Reads/HV_020_04_forward.fq.gz
mv HV020_05_1.fq.gz Raw_Reads/HV_020_05_forward.fq.gz
mv HV020_06_1.fq.gz Raw_Reads/HV_020_06_forward.fq.gz
mv HV020_07_1.fq.gz Raw_Reads/HV_020_07_forward.fq.gz
mv HV020_08_1.fq.gz Raw_Reads/HV_020_08_forward.fq.gz
mv HV020_09_1.fq.gz Raw_Reads/HV_020_09_forward.fq.gz
mv HV020_10_1.fq.gz Raw_Reads/HV_020_10_forward.fq.gz
mv HV020_11_1.fq.gz Raw_Reads/HV_020_11_forward.fq.gz
mv HV020_12_1.fq.gz Raw_Reads/HV_020_12_forward.fq.gz
mv HV020_13_1.fq.gz Raw_Reads/HV_020_13_forward.fq.gz
mv HV020_14_1.fq.gz Raw_Reads/HV_020_14_forward.fq.gz
mv HV020_15_1.fq.gz Raw_Reads/HV_020_15_forward.fq.gz
mv HV020_16_1.fq.gz Raw_Reads/HV_020_16_forward.fq.gz
mv HV020_18_1.fq.gz Raw_Reads/HV_020_18_forward.fq.gz
mv HV020_19_1.fq.gz Raw_Reads/HV_020_19_forward.fq.gz
mv HV020_20_1.fq.gz Raw_Reads/HV_020_20_forward.fq.gz
mv HV020_21_1.fq.gz Raw_Reads/HV_020_21_forward.fq.gz
mv HV020_22_1.fq.gz Raw_Reads/HV_020_22_forward.fq.gz
mv HV020_23_1.fq.gz Raw_Reads/HV_020_23_forward.fq.gz
mv HV020_25_1.fq.gz Raw_Reads/HV_020_25_forward.fq.gz
mv HV020_27_1.fq.gz Raw_Reads/HV_020_27_forward.fq.gz
mv HV021_01_1.fq.gz Raw_Reads/HV_021_01_forward.fq.gz
mv HV021_02_1.fq.gz Raw_Reads/HV_021_02_forward.fq.gz
mv HV021_03_1.fq.gz Raw_Reads/HV_021_03_forward.fq.gz
mv HV021_04_1.fq.gz Raw_Reads/HV_021_04_forward.fq.gz
mv HV021_05_1.fq.gz Raw_Reads/HV_021_05_forward.fq.gz
mv HV021_06_1.fq.gz Raw_Reads/HV_021_06_forward.fq.gz
mv HV021_07_1.fq.gz Raw_Reads/HV_021_07_forward.fq.gz
mv HV021_08_1.fq.gz Raw_Reads/HV_021_08_forward.fq.gz
mv HV021_09_1.fq.gz Raw_Reads/HV_021_09_forward.fq.gz
mv HV021_10_1.fq.gz Raw_Reads/HV_021_10_forward.fq.gz
mv HV021_11_1.fq.gz Raw_Reads/HV_021_11_forward.fq.gz
mv HV021_12_1.fq.gz Raw_Reads/HV_021_12_forward.fq.gz
mv HV021_13_1.fq.gz Raw_Reads/HV_021_13_forward.fq.gz
mv HV021_14_1.fq.gz Raw_Reads/HV_021_14_forward.fq.gz
mv HV021_15_1.fq.gz Raw_Reads/HV_021_15_forward.fq.gz
mv HV021_16_1.fq.gz Raw_Reads/HV_021_16_forward.fq.gz
mv HV021_18_1.fq.gz Raw_Reads/HV_021_18_forward.fq.gz
mv HV021_19_1.fq.gz Raw_Reads/HV_021_19_forward.fq.gz
mv HV021_20_1.fq.gz Raw_Reads/HV_021_20_forward.fq.gz
mv HV021_21_1.fq.gz Raw_Reads/HV_021_21_forward.fq.gz
mv HV021_22_1.fq.gz Raw_Reads/HV_021_22_forward.fq.gz
mv HV021_23_1.fq.gz Raw_Reads/HV_021_23_forward.fq.gz
mv HV021_25_1.fq.gz Raw_Reads/HV_021_25_forward.fq.gz
mv HV021_27_1.fq.gz Raw_Reads/HV_021_27_forward.fq.gz
mv HV_022_CKDL190143346-1a-1_H75V2BBXX_L2_1.fq.gz Raw_Reads/HV_022_01_forward.fq.gz
mv HV_022_CKDL190143346-1a-2_H75V2BBXX_L2_1.fq.gz Raw_Reads/HV_022_02_forward.fq.gz
mv HV_022_CKDL190143346-1a-3_H75V2BBXX_L2_1.fq.gz Raw_Reads/HV_022_03_forward.fq.gz
mv HV_022_CKDL190143346-1a-4_H75V2BBXX_L2_1.fq.gz Raw_Reads/HV_022_04_forward.fq.gz
mv HV_022_CKDL190143346-1a-5_H75V2BBXX_L2_1.fq.gz Raw_Reads/HV_022_05_forward.fq.gz
mv HV_022_CKDL190143346-1a-6_H75V2BBXX_L2_1.fq.gz Raw_Reads/HV_022_06_forward.fq.gz
mv HV_022_CKDL190143346-1a-7_H75V2BBXX_L2_1.fq.gz Raw_Reads/HV_022_07_forward.fq.gz
mv HV_022_CKDL190143346-1a-8_H75V2BBXX_L2_1.fq.gz Raw_Reads/HV_022_08_forward.fq.gz
mv HV_022_CKDL190143346-1a-9_H75V2BBXX_L2_1.fq.gz Raw_Reads/HV_022_09_forward.fq.gz
mv HV_022_CKDL190143346-1a-10_H75V2BBXX_L2_1.fq.gz Raw_Reads/HV_022_10_forward.fq.gz
mv HV_022_CKDL190143346-1a-11_H75V2BBXX_L2_1.fq.gz Raw_Reads/HV_022_11_forward.fq.gz
mv HV_022_CKDL190143346-1a-12_H75V2BBXX_L2_1.fq.gz Raw_Reads/HV_022_12_forward.fq.gz
mv HV_022_CKDL190143346-1a-13_H75V2BBXX_L2_1.fq.gz Raw_Reads/HV_022_13_forward.fq.gz
mv HV_022_CKDL190143346-1a-14_H75V2BBXX_L2_1.fq.gz Raw_Reads/HV_022_14_forward.fq.gz
mv HV_022_CKDL190143346-1a-15_H75V2BBXX_L2_1.fq.gz Raw_Reads/HV_022_15_forward.fq.gz
mv HV_022_CKDL190143346-1a-16_H75V2BBXX_L2_1.fq.gz Raw_Reads/HV_022_16_forward.fq.gz
mv HV_022_CKDL190143346-1a-18_H75V2BBXX_L2_1.fq.gz Raw_Reads/HV_022_18_forward.fq.gz
mv HV_022_CKDL190143346-1a-19_H75V2BBXX_L2_1.fq.gz Raw_Reads/HV_022_19_forward.fq.gz
mv HV_022_CKDL190143346-1a-20_H75V2BBXX_L2_1.fq.gz Raw_Reads/HV_022_20_forward.fq.gz
mv HV_022_CKDL190143346-1a-21_H75V2BBXX_L2_1.fq.gz Raw_Reads/HV_022_21_forward.fq.gz
mv HV_022_CKDL190143346-1a-22_H75V2BBXX_L2_1.fq.gz Raw_Reads/HV_022_22_forward.fq.gz
mv HV_022_CKDL190143346-1a-23_H75V2BBXX_L2_1.fq.gz Raw_Reads/HV_022_23_forward.fq.gz
mv HV_022_CKDL190143346-1a-25_H75V2BBXX_L2_1.fq.gz Raw_Reads/HV_022_25_forward.fq.gz
mv HV_022_CKDL190143346-1a-27_H75V2BBXX_L2_1.fq.gz Raw_Reads/HV_022_27_forward.fq.gz
mv HV_023_CKDL190143347-1a-1_H75V2BBXX_L3_1.fq.gz Raw_Reads/HV_023_01_forward.fq.gz
mv HV_023_CKDL190143347-1a-2_H75V2BBXX_L3_1.fq.gz Raw_Reads/HV_023_02_forward.fq.gz
mv HV_023_CKDL190143347-1a-3_H75V2BBXX_L3_1.fq.gz Raw_Reads/HV_023_03_forward.fq.gz
mv HV_023_CKDL190143347-1a-4_H75V2BBXX_L3_1.fq.gz Raw_Reads/HV_023_04_forward.fq.gz
mv HV_023_CKDL190143347-1a-5_H75V2BBXX_L3_1.fq.gz Raw_Reads/HV_023_05_forward.fq.gz
mv HV_023_CKDL190143347-1a-6_H75V2BBXX_L3_1.fq.gz Raw_Reads/HV_023_06_forward.fq.gz
mv HV_023_CKDL190143347-1a-7_H75V2BBXX_L3_1.fq.gz Raw_Reads/HV_023_07_forward.fq.gz
mv HV_023_CKDL190143347-1a-8_H75V2BBXX_L3_1.fq.gz Raw_Reads/HV_023_08_forward.fq.gz
mv HV_023_CKDL190143347-1a-9_H75V2BBXX_L3_1.fq.gz Raw_Reads/HV_023_09_forward.fq.gz
mv HV_023_CKDL190143347-1a-10_H75V2BBXX_L3_1.fq.gz Raw_Reads/HV_023_10_forward.fq.gz
mv HV_023_CKDL190143347-1a-11_H75V2BBXX_L3_1.fq.gz Raw_Reads/HV_023_11_forward.fq.gz
mv HV_023_CKDL190143347-1a-12_H75V2BBXX_L3_1.fq.gz Raw_Reads/HV_023_12_forward.fq.gz
mv HV_023_CKDL190143347-1a-13_H75V2BBXX_L3_1.fq.gz Raw_Reads/HV_023_13_forward.fq.gz
mv HV_023_CKDL190143347-1a-14_H75V2BBXX_L3_1.fq.gz Raw_Reads/HV_023_14_forward.fq.gz
mv HV_023_CKDL190143347-1a-15_H75V2BBXX_L3_1.fq.gz Raw_Reads/HV_023_15_forward.fq.gz
mv HV_023_CKDL190143347-1a-16_H75V2BBXX_L3_1.fq.gz Raw_Reads/HV_023_16_forward.fq.gz
mv HV_023_CKDL190143347-1a-18_H75V2BBXX_L3_1.fq.gz Raw_Reads/HV_023_18_forward.fq.gz
mv HV_023_CKDL190143347-1a-19_H75V2BBXX_L3_1.fq.gz Raw_Reads/HV_023_19_forward.fq.gz
mv HV_023_CKDL190143347-1a-20_H75V2BBXX_L3_1.fq.gz Raw_Reads/HV_023_20_forward.fq.gz
mv HV_023_CKDL190143347-1a-21_H75V2BBXX_L3_1.fq.gz Raw_Reads/HV_023_21_forward.fq.gz
mv HV_023_CKDL190143347-1a-22_H75V2BBXX_L3_1.fq.gz Raw_Reads/HV_023_22_forward.fq.gz
mv HV_023_CKDL190143347-1a-23_H75V2BBXX_L3_1.fq.gz Raw_Reads/HV_023_23_forward.fq.gz
mv HV_023_CKDL190143347-1a-25_H75V2BBXX_L3_1.fq.gz Raw_Reads/HV_023_25_forward.fq.gz
mv HV_023_CKDL190143347-1a-27_H75V2BBXX_L3_1.fq.gz Raw_Reads/HV_023_27_forward.fq.gz
mv HV_024_CKDL190143348-1a-1_H75V2BBXX_L4_1.fq.gz Raw_Reads/HV_024_01_forward.fq.gz
mv HV_024_CKDL190143348-1a-2_H75V2BBXX_L4_1.fq.gz Raw_Reads/HV_024_02_forward.fq.gz
mv HV_024_CKDL190143348-1a-3_H75V2BBXX_L4_1.fq.gz Raw_Reads/HV_024_03_forward.fq.gz
mv HV_024_CKDL190143348-1a-4_H75V2BBXX_L4_1.fq.gz Raw_Reads/HV_024_04_forward.fq.gz
mv HV_024_CKDL190143348-1a-5_H75V2BBXX_L4_1.fq.gz Raw_Reads/HV_024_05_forward.fq.gz
mv HV_024_CKDL190143348-1a-6_H75V2BBXX_L4_1.fq.gz Raw_Reads/HV_024_06_forward.fq.gz
mv HV_024_CKDL190143348-1a-7_H75V2BBXX_L4_1.fq.gz Raw_Reads/HV_024_07_forward.fq.gz
mv HV_024_CKDL190143348-1a-8_H75V2BBXX_L4_1.fq.gz Raw_Reads/HV_024_08_forward.fq.gz
mv HV_024_CKDL190143348-1a-9_H75V2BBXX_L4_1.fq.gz Raw_Reads/HV_024_09_forward.fq.gz
mv HV_024_CKDL190143348-1a-10_H75V2BBXX_L4_1.fq.gz Raw_Reads/HV_024_10_forward.fq.gz
mv HV_024_CKDL190143348-1a-11_H75V2BBXX_L4_1.fq.gz Raw_Reads/HV_024_11_forward.fq.gz
mv HV_024_CKDL190143348-1a-12_H75V2BBXX_L4_1.fq.gz Raw_Reads/HV_024_12_forward.fq.gz
mv HV_024_CKDL190143348-1a-13_H75V2BBXX_L4_1.fq.gz Raw_Reads/HV_024_13_forward.fq.gz
mv HV_024_CKDL190143348-1a-14_H75V2BBXX_L4_1.fq.gz Raw_Reads/HV_024_14_forward.fq.gz
mv HV_024_CKDL190143348-1a-15_H75V2BBXX_L4_1.fq.gz Raw_Reads/HV_024_15_forward.fq.gz
mv HV_024_CKDL190143348-1a-16_H75V2BBXX_L4_1.fq.gz Raw_Reads/HV_024_16_forward.fq.gz
mv HV_024_CKDL190143348-1a-18_H75V2BBXX_L4_1.fq.gz Raw_Reads/HV_024_18_forward.fq.gz
mv HV_024_CKDL190143348-1a-19_H75V2BBXX_L4_1.fq.gz Raw_Reads/HV_024_19_forward.fq.gz
mv HV_024_CKDL190143348-1a-20_H75V2BBXX_L4_1.fq.gz Raw_Reads/HV_024_20_forward.fq.gz
mv HV_024_CKDL190143348-1a-21_H75V2BBXX_L4_1.fq.gz Raw_Reads/HV_024_21_forward.fq.gz
mv HV_024_CKDL190143348-1a-22_H75V2BBXX_L4_1.fq.gz Raw_Reads/HV_024_22_forward.fq.gz
mv HV_024_CKDL190143348-1a-23_H75V2BBXX_L4_1.fq.gz Raw_Reads/HV_024_23_forward.fq.gz
mv HV_024_CKDL190143348-1a-25_H75V2BBXX_L4_1.fq.gz Raw_Reads/HV_024_25_forward.fq.gz
mv HV_024_CKDL190143348-1a-27_H75V2BBXX_L4_1.fq.gz Raw_Reads/HV_024_27_forward.fq.gz
mv HV_025_CKDL190143349-1a-1_H75V2BBXX_L5_1.fq.gz Raw_Reads/HV_025_01_forward.fq.gz
mv HV_025_CKDL190143349-1a-2_H75V2BBXX_L5_1.fq.gz Raw_Reads/HV_025_02_forward.fq.gz
mv HV_025_CKDL190143349-1a-3_H75V2BBXX_L5_1.fq.gz Raw_Reads/HV_025_03_forward.fq.gz
mv HV_025_CKDL190143349-1a-4_H75V2BBXX_L5_1.fq.gz Raw_Reads/HV_025_04_forward.fq.gz
mv HV_025_CKDL190143349-1a-5_H75V2BBXX_L5_1.fq.gz Raw_Reads/HV_025_05_forward.fq.gz
mv HV_025_CKDL190143349-1a-6_H75V2BBXX_L5_1.fq.gz Raw_Reads/HV_025_06_forward.fq.gz
mv HV_025_CKDL190143349-1a-7_H75V2BBXX_L5_1.fq.gz Raw_Reads/HV_025_07_forward.fq.gz
mv HV_025_CKDL190143349-1a-8_H75V2BBXX_L5_1.fq.gz Raw_Reads/HV_025_08_forward.fq.gz
mv HV_025_CKDL190143349-1a-9_H75V2BBXX_L5_1.fq.gz Raw_Reads/HV_025_09_forward.fq.gz
mv HV_025_CKDL190143349-1a-10_H75V2BBXX_L5_1.fq.gz Raw_Reads/HV_025_10_forward.fq.gz
mv HV_025_CKDL190143349-1a-11_H75V2BBXX_L5_1.fq.gz Raw_Reads/HV_025_11_forward.fq.gz
mv HV_025_CKDL190143349-1a-12_H75V2BBXX_L5_1.fq.gz Raw_Reads/HV_025_12_forward.fq.gz
mv HV_025_CKDL190143349-1a-13_H75V2BBXX_L5_1.fq.gz Raw_Reads/HV_025_13_forward.fq.gz
mv HV_025_CKDL190143349-1a-14_H75V2BBXX_L5_1.fq.gz Raw_Reads/HV_025_14_forward.fq.gz
mv HV_025_CKDL190143349-1a-15_H75V2BBXX_L5_1.fq.gz Raw_Reads/HV_025_15_forward.fq.gz
mv HV_025_CKDL190143349-1a-16_H75V2BBXX_L5_1.fq.gz Raw_Reads/HV_025_16_forward.fq.gz
mv HV_025_CKDL190143349-1a-18_H75V2BBXX_L5_1.fq.gz Raw_Reads/HV_025_18_forward.fq.gz
mv HV_025_CKDL190143349-1a-19_H75V2BBXX_L5_1.fq.gz Raw_Reads/HV_025_19_forward.fq.gz
mv HV_025_CKDL190143349-1a-20_H75V2BBXX_L5_1.fq.gz Raw_Reads/HV_025_20_forward.fq.gz
mv HV_025_CKDL190143349-1a-21_H75V2BBXX_L5_1.fq.gz Raw_Reads/HV_025_21_forward.fq.gz
mv HV_025_CKDL190143349-1a-22_H75V2BBXX_L5_1.fq.gz Raw_Reads/HV_025_22_forward.fq.gz
mv HV_025_CKDL190143349-1a-23_H75V2BBXX_L5_1.fq.gz Raw_Reads/HV_025_23_forward.fq.gz
mv HV_025_CKDL190143349-1a-25_H75V2BBXX_L5_1.fq.gz Raw_Reads/HV_025_25_forward.fq.gz
mv HV_025_CKDL190143349-1a-27_H75V2BBXX_L5_1.fq.gz Raw_Reads/HV_025_27_forward.fq.gz
mv HV_026_CKDL190144758-1a-1_H75W3BBXX_L4_1.fq.gz Raw_Reads/HV_026_01_forward.fq.gz
mv HV_026_CKDL190144758-1a-2_H75W3BBXX_L4_1.fq.gz Raw_Reads/HV_026_02_forward.fq.gz
mv HV_026_CKDL190144758-1a-3_H75W3BBXX_L4_1.fq.gz Raw_Reads/HV_026_03_forward.fq.gz
mv HV_026_CKDL190144758-1a-4_H75W3BBXX_L4_1.fq.gz Raw_Reads/HV_026_04_forward.fq.gz
mv HV_026_CKDL190144758-1a-5_H75W3BBXX_L4_1.fq.gz Raw_Reads/HV_026_05_forward.fq.gz
mv HV_026_CKDL190144758-1a-6_H75W3BBXX_L4_1.fq.gz Raw_Reads/HV_026_06_forward.fq.gz
mv HV_026_CKDL190144758-1a-7_H75W3BBXX_L4_1.fq.gz Raw_Reads/HV_026_07_forward.fq.gz
mv HV_026_CKDL190144758-1a-8_H75W3BBXX_L4_1.fq.gz Raw_Reads/HV_026_08_forward.fq.gz
mv HV_026_CKDL190144758-1a-9_H75W3BBXX_L4_1.fq.gz Raw_Reads/HV_026_09_forward.fq.gz
mv HV_026_CKDL190144758-1a-10_H75W3BBXX_L4_1.fq.gz Raw_Reads/HV_026_10_forward.fq.gz
mv HV_026_CKDL190144758-1a-11_H75W3BBXX_L4_1.fq.gz Raw_Reads/HV_026_11_forward.fq.gz
mv HV_026_CKDL190144758-1a-12_H75W3BBXX_L4_1.fq.gz Raw_Reads/HV_026_12_forward.fq.gz
mv HV_026_CKDL190144758-1a-13_H75W3BBXX_L4_1.fq.gz Raw_Reads/HV_026_13_forward.fq.gz
mv HV_026_CKDL190144758-1a-14_H75W3BBXX_L4_1.fq.gz Raw_Reads/HV_026_14_forward.fq.gz
mv HV_026_CKDL190144758-1a-15_H75W3BBXX_L4_1.fq.gz Raw_Reads/HV_026_15_forward.fq.gz
mv HV_026_CKDL190144758-1a-16_H75W3BBXX_L4_1.fq.gz Raw_Reads/HV_026_16_forward.fq.gz
mv HV_026_CKDL190144758-1a-18_H75W3BBXX_L4_1.fq.gz Raw_Reads/HV_026_18_forward.fq.gz
mv HV_026_CKDL190144758-1a-19_H75W3BBXX_L4_1.fq.gz Raw_Reads/HV_026_19_forward.fq.gz
mv HV_026_CKDL190144758-1a-20_H75W3BBXX_L4_1.fq.gz Raw_Reads/HV_026_20_forward.fq.gz
mv HV_026_CKDL190144758-1a-21_H75W3BBXX_L4_1.fq.gz Raw_Reads/HV_026_21_forward.fq.gz
mv HV_026_CKDL190144758-1a-22_H75W3BBXX_L4_1.fq.gz Raw_Reads/HV_026_22_forward.fq.gz
mv HV_026_CKDL190144758-1a-23_H75W3BBXX_L4_1.fq.gz Raw_Reads/HV_026_23_forward.fq.gz
mv HV_026_CKDL190144758-1a-25_H75W3BBXX_L4_1.fq.gz Raw_Reads/HV_026_25_forward.fq.gz
mv HV_026_CKDL190144758-1a-27_H75W3BBXX_L4_1.fq.gz Raw_Reads/HV_026_27_forward.fq.gz
mv HV_027_CKDL190144759-1a-1_H75TMBBXX_L8_1.fq.gz Raw_Reads/HV_027_01_forward.fq.gz
mv HV_027_CKDL190144759-1a-2_H75TMBBXX_L8_1.fq.gz Raw_Reads/HV_027_02_forward.fq.gz
mv HV_027_CKDL190144759-1a-3_H75TMBBXX_L8_1.fq.gz Raw_Reads/HV_027_03_forward.fq.gz
mv HV_027_CKDL190144759-1a-4_H75TMBBXX_L8_1.fq.gz Raw_Reads/HV_027_04_forward.fq.gz
mv HV_027_CKDL190144759-1a-5_H75TMBBXX_L8_1.fq.gz Raw_Reads/HV_027_05_forward.fq.gz
mv HV_027_CKDL190144759-1a-6_H75TMBBXX_L8_1.fq.gz Raw_Reads/HV_027_06_forward.fq.gz
mv HV_027_CKDL190144759-1a-7_H75TMBBXX_L8_1.fq.gz Raw_Reads/HV_027_07_forward.fq.gz
mv HV_027_CKDL190144759-1a-8_H75TMBBXX_L8_1.fq.gz Raw_Reads/HV_027_08_forward.fq.gz
mv HV_027_CKDL190144759-1a-9_H75TMBBXX_L8_1.fq.gz Raw_Reads/HV_027_09_forward.fq.gz
mv HV_027_CKDL190144759-1a-10_H75TMBBXX_L8_1.fq.gz Raw_Reads/HV_027_10_forward.fq.gz
mv HV_027_CKDL190144759-1a-11_H75TMBBXX_L8_1.fq.gz Raw_Reads/HV_027_11_forward.fq.gz
mv HV_027_CKDL190144759-1a-12_H75TMBBXX_L8_1.fq.gz Raw_Reads/HV_027_12_forward.fq.gz
mv HV_027_CKDL190144759-1a-13_H75TMBBXX_L8_1.fq.gz Raw_Reads/HV_027_13_forward.fq.gz
mv HV_027_CKDL190144759-1a-14_H75TMBBXX_L8_1.fq.gz Raw_Reads/HV_027_14_forward.fq.gz
mv HV_027_CKDL190144759-1a-15_H75TMBBXX_L8_1.fq.gz Raw_Reads/HV_027_15_forward.fq.gz
mv HV_027_CKDL190144759-1a-16_H75TMBBXX_L8_1.fq.gz Raw_Reads/HV_027_16_forward.fq.gz
mv HV_027_CKDL190144759-1a-18_H75TMBBXX_L8_1.fq.gz Raw_Reads/HV_027_18_forward.fq.gz
mv HV_027_CKDL190144759-1a-19_H75TMBBXX_L8_1.fq.gz Raw_Reads/HV_027_19_forward.fq.gz
mv HV_027_CKDL190144759-1a-20_H75TMBBXX_L8_1.fq.gz Raw_Reads/HV_027_20_forward.fq.gz
mv HV_027_CKDL190144759-1a-21_H75TMBBXX_L8_1.fq.gz Raw_Reads/HV_027_21_forward.fq.gz
mv HV_027_CKDL190144759-1a-22_H75TMBBXX_L8_1.fq.gz Raw_Reads/HV_027_22_forward.fq.gz
mv HV_027_CKDL190144759-1a-23_H75TMBBXX_L8_1.fq.gz Raw_Reads/HV_027_23_forward.fq.gz
mv HV_027_CKDL190144759-1a-25_H75TMBBXX_L8_1.fq.gz Raw_Reads/HV_027_25_forward.fq.gz
mv HV_027_CKDL190144759-1a-27_H75TMBBXX_L8_1.fq.gz Raw_Reads/HV_027_27_forward.fq.gz
mv HV028_1_CKDL200148804-1a-1_H7HJWBBXX_L3_1.fq.gz Raw_Reads/HV_028_01_forward.fq.gz
mv HV028_2_CKDL200148804-1a-2_H7HJWBBXX_L3_1.fq.gz Raw_Reads/HV_028_02_forward.fq.gz
mv HV028_3_CKDL200148804-1a-3_H7HJWBBXX_L3_1.fq.gz Raw_Reads/HV_028_03_forward.fq.gz
mv HV028_4_CKDL200148804-1a-4_H7HJWBBXX_L3_1.fq.gz Raw_Reads/HV_028_04_forward.fq.gz
mv HV028_5_CKDL200148804-1a-5_H7HJWBBXX_L3_1.fq.gz Raw_Reads/HV_028_05_forward.fq.gz
mv HV028_6_CKDL200148804-1a-6_H7HJWBBXX_L3_1.fq.gz Raw_Reads/HV_028_06_forward.fq.gz
mv HV028_7_CKDL200148804-1a-7_H7HJWBBXX_L3_1.fq.gz Raw_Reads/HV_028_07_forward.fq.gz
mv HV028_8_CKDL200148804-1a-8_H7HJWBBXX_L3_1.fq.gz Raw_Reads/HV_028_08_forward.fq.gz
mv HV028_9_CKDL200148804-1a-9_H7HJWBBXX_L3_1.fq.gz Raw_Reads/HV_028_09_forward.fq.gz
mv HV028_10_CKDL200148804-1a-10_H7HJWBBXX_L3_1.fq.gz Raw_Reads/HV_028_10_forward.fq.gz
mv HV028_11_CKDL200148804-1a-11_H7HJWBBXX_L3_1.fq.gz Raw_Reads/HV_028_11_forward.fq.gz
mv HV028_12_CKDL200148804-1a-12_H7HJWBBXX_L3_1.fq.gz Raw_Reads/HV_028_12_forward.fq.gz
mv HV028_13_CKDL200148804-1a-13_H7HJWBBXX_L3_1.fq.gz Raw_Reads/HV_028_13_forward.fq.gz
mv HV028_14_CKDL200148804-1a-14_H7HJWBBXX_L3_1.fq.gz Raw_Reads/HV_028_14_forward.fq.gz
mv HV028_15_CKDL200148804-1a-15_H7HJWBBXX_L3_1.fq.gz Raw_Reads/HV_028_15_forward.fq.gz
mv HV028_16_CKDL200148804-1a-16_H7HJWBBXX_L3_1.fq.gz Raw_Reads/HV_028_16_forward.fq.gz
mv HV028_18_CKDL200148804-1a-18_H7HJWBBXX_L3_1.fq.gz Raw_Reads/HV_028_18_forward.fq.gz
mv HV028_19_CKDL200148804-1a-19_H7HJWBBXX_L3_1.fq.gz Raw_Reads/HV_028_19_forward.fq.gz
mv HV028_20_CKDL200148804-1a-20_H7HJWBBXX_L3_1.fq.gz Raw_Reads/HV_028_20_forward.fq.gz
mv HV028_21_CKDL200148804-1a-21_H7HJWBBXX_L3_1.fq.gz Raw_Reads/HV_028_21_forward.fq.gz
mv HV028_22_CKDL200148804-1a-22_H7HJWBBXX_L3_1.fq.gz Raw_Reads/HV_028_22_forward.fq.gz
mv HV028_23_CKDL200148804-1a-23_H7HJWBBXX_L3_1.fq.gz Raw_Reads/HV_028_23_forward.fq.gz
mv HV028_25_CKDL200148804-1a-25_H7HJWBBXX_L3_1.fq.gz Raw_Reads/HV_028_25_forward.fq.gz
mv HV028_27_CKDL200148804-1a-27_H7HJWBBXX_L3_1.fq.gz Raw_Reads/HV_028_27_forward.fq.gz
mv HV029_1_CKDL200148805-1a-1_H7HJWBBXX_L4_1.fq.gz Raw_Reads/HV_029_01_forward.fq.gz
mv HV029_2_CKDL200148805-1a-2_H7HJWBBXX_L4_1.fq.gz Raw_Reads/HV_029_02_forward.fq.gz
mv HV029_3_CKDL200148805-1a-3_H7HJWBBXX_L4_1.fq.gz Raw_Reads/HV_029_03_forward.fq.gz
mv HV029_4_CKDL200148805-1a-4_H7HJWBBXX_L4_1.fq.gz Raw_Reads/HV_029_04_forward.fq.gz
mv HV029_5_CKDL200148805-1a-5_H7HJWBBXX_L4_1.fq.gz Raw_Reads/HV_029_05_forward.fq.gz
mv HV029_6_CKDL200148805-1a-6_H7HJWBBXX_L4_1.fq.gz Raw_Reads/HV_029_06_forward.fq.gz
mv HV029_7_CKDL200148805-1a-7_H7HJWBBXX_L4_1.fq.gz Raw_Reads/HV_029_07_forward.fq.gz
mv HV029_8_CKDL200148805-1a-8_H7HJWBBXX_L4_1.fq.gz Raw_Reads/HV_029_08_forward.fq.gz
mv HV029_9_CKDL200148805-1a-9_H7HJWBBXX_L4_1.fq.gz Raw_Reads/HV_029_09_forward.fq.gz
mv HV029_10_CKDL200148805-1a-10_H7HJWBBXX_L4_1.fq.gz Raw_Reads/HV_029_10_forward.fq.gz
mv HV029_11_CKDL200148805-1a-11_H7HJWBBXX_L4_1.fq.gz Raw_Reads/HV_029_11_forward.fq.gz
mv HV029_12_CKDL200148805-1a-12_H7HJWBBXX_L4_1.fq.gz Raw_Reads/HV_029_12_forward.fq.gz
mv HV029_13_CKDL200148805-1a-13_H7HJWBBXX_L4_1.fq.gz Raw_Reads/HV_029_13_forward.fq.gz
mv HV029_14_CKDL200148805-1a-14_H7HJWBBXX_L4_1.fq.gz Raw_Reads/HV_029_14_forward.fq.gz
mv HV029_15_CKDL200148805-1a-15_H7HJWBBXX_L4_1.fq.gz Raw_Reads/HV_029_15_forward.fq.gz
mv HV029_16_CKDL200148805-1a-16_H7HJWBBXX_L4_1.fq.gz Raw_Reads/HV_029_16_forward.fq.gz
mv HV029_18_CKDL200148805-1a-18_H7HJWBBXX_L4_1.fq.gz Raw_Reads/HV_029_18_forward.fq.gz
mv HV029_19_CKDL200148805-1a-19_H7HJWBBXX_L4_1.fq.gz Raw_Reads/HV_029_19_forward.fq.gz
mv HV029_20_CKDL200148805-1a-20_H7HJWBBXX_L4_1.fq.gz Raw_Reads/HV_029_20_forward.fq.gz
mv HV029_21_CKDL200148805-1a-21_H7HJWBBXX_L4_1.fq.gz Raw_Reads/HV_029_21_forward.fq.gz
mv HV029_22_CKDL200148805-1a-22_H7HJWBBXX_L4_1.fq.gz Raw_Reads/HV_029_22_forward.fq.gz
mv HV029_23_CKDL200148805-1a-23_H7HJWBBXX_L4_1.fq.gz Raw_Reads/HV_029_23_forward.fq.gz
mv HV029_25_CKDL200148805-1a-25_H7HJWBBXX_L4_1.fq.gz Raw_Reads/HV_029_25_forward.fq.gz
mv HV029_27_CKDL200148805-1a-27_H7HJWBBXX_L4_1.fq.gz Raw_Reads/HV_029_27_forward.fq.gz
mv HV030_1_CKDL200148806-1a-1_H7HJWBBXX_L5_1.fq.gz Raw_Reads/HV_030_01_forward.fq.gz
mv HV030_2_CKDL200148806-1a-2_H7HJWBBXX_L5_1.fq.gz Raw_Reads/HV_030_02_forward.fq.gz
mv HV030_3_CKDL200148806-1a-3_H7HJWBBXX_L5_1.fq.gz Raw_Reads/HV_030_03_forward.fq.gz
mv HV030_4_CKDL200148806-1a-4_H7HJWBBXX_L5_1.fq.gz Raw_Reads/HV_030_04_forward.fq.gz
mv HV030_5_CKDL200148806-1a-5_H7HJWBBXX_L5_1.fq.gz Raw_Reads/HV_030_05_forward.fq.gz
mv HV030_6_CKDL200148806-1a-6_H7HJWBBXX_L5_1.fq.gz Raw_Reads/HV_030_06_forward.fq.gz
mv HV030_7_CKDL200148806-1a-7_H7HJWBBXX_L5_1.fq.gz Raw_Reads/HV_030_07_forward.fq.gz
mv HV030_8_CKDL200148806-1a-8_H7HJWBBXX_L5_1.fq.gz Raw_Reads/HV_030_08_forward.fq.gz
mv HV030_9_CKDL200148806-1a-9_H7HJWBBXX_L5_1.fq.gz Raw_Reads/HV_030_09_forward.fq.gz
mv HV030_10_CKDL200148806-1a-10_H7HJWBBXX_L5_1.fq.gz Raw_Reads/HV_030_10_forward.fq.gz
mv HV030_11_CKDL200148806-1a-11_H7HJWBBXX_L5_1.fq.gz Raw_Reads/HV_030_11_forward.fq.gz
mv HV030_12_CKDL200148806-1a-12_H7HJWBBXX_L5_1.fq.gz Raw_Reads/HV_030_12_forward.fq.gz
mv HV030_13_CKDL200148806-1a-13_H7HJWBBXX_L5_1.fq.gz Raw_Reads/HV_030_13_forward.fq.gz
mv HV030_14_CKDL200148806-1a-14_H7HJWBBXX_L5_1.fq.gz Raw_Reads/HV_030_14_forward.fq.gz
mv HV030_15_CKDL200148806-1a-15_H7HJWBBXX_L5_1.fq.gz Raw_Reads/HV_030_15_forward.fq.gz
mv HV030_16_CKDL200148806-1a-16_H7HJWBBXX_L5_1.fq.gz Raw_Reads/HV_030_16_forward.fq.gz
mv HV030_18_CKDL200148806-1a-18_H7HJWBBXX_L5_1.fq.gz Raw_Reads/HV_030_18_forward.fq.gz
mv HV030_19_CKDL200148806-1a-19_H7HJWBBXX_L5_1.fq.gz Raw_Reads/HV_030_19_forward.fq.gz
mv HV030_20_CKDL200148806-1a-20_H7HJWBBXX_L5_1.fq.gz Raw_Reads/HV_030_20_forward.fq.gz
mv HV030_21_CKDL200148806-1a-21_H7HJWBBXX_L5_1.fq.gz Raw_Reads/HV_030_21_forward.fq.gz
mv HV030_22_CKDL200148806-1a-22_H7HJWBBXX_L5_1.fq.gz Raw_Reads/HV_030_22_forward.fq.gz
mv HV030_23_CKDL200148806-1a-23_H7HJWBBXX_L5_1.fq.gz Raw_Reads/HV_030_23_forward.fq.gz
mv HV030_25_CKDL200148806-1a-25_H7HJWBBXX_L5_1.fq.gz Raw_Reads/HV_030_25_forward.fq.gz
mv HV030_27_CKDL200148806-1a-27_H7HJWBBXX_L5_1.fq.gz Raw_Reads/HV_030_27_forward.fq.gz

############################################
## -------------------------------------- ##
## ------------ Reverse Reads ----------- ##
## -------------------------------------- ##
############################################

mv EG_Sflab_241_1_USPD16092227-1_HV5CHBBXX_L2_2.fq.gz Raw_Reads/HV_001_01_reverse.fq.gz
mv EG_Sflab_241_2_USPD16092227-2_HV5CHBBXX_L2_2.fq.gz Raw_Reads/HV_001_02_reverse.fq.gz
mv EG_Sflab_241_3_USPD16092227-3_HV5CHBBXX_L2_2.fq.gz Raw_Reads/HV_001_03_reverse.fq.gz
mv EG_Sflab_241_4_USPD16092227-4_HV5CHBBXX_L2_2.fq.gz Raw_Reads/HV_001_04_reverse.fq.gz
mv EG_Sflab_241_5_USPD16092227-5_HV5CHBBXX_L2_2.fq.gz Raw_Reads/HV_001_05_reverse.fq.gz
mv EG_Sflab_241_6_USPD16092227-6_HV5CHBBXX_L2_2.fq.gz Raw_Reads/HV_001_06_reverse.fq.gz
mv EG_Sflab_241_7_USPD16092227-7_HV5CHBBXX_L2_2.fq.gz Raw_Reads/HV_001_07_reverse.fq.gz
mv EG_Sflab_241_8_USPD16092227-8_HV5CHBBXX_L2_2.fq.gz Raw_Reads/HV_001_08_reverse.fq.gz
mv EG_Sflab_241_9_USPD16092227-9_HV5CHBBXX_L2_2.fq.gz Raw_Reads/HV_001_09_reverse.fq.gz
mv EG_Sflab_241_10_USPD16092227-10_HV5CHBBXX_L2_2.fq.gz Raw_Reads/HV_001_10_reverse.fq.gz
mv EG_Sflab_241_11_USPD16092227-11_HV5CHBBXX_L2_2.fq.gz Raw_Reads/HV_001_11_reverse.fq.gz
mv EG_Sflab_241_12_USPD16092227-12_HV5CHBBXX_L2_2.fq.gz Raw_Reads/HV_001_12_reverse.fq.gz
mv EG_Sflab_241_13_USPD16092227-13_HV5CHBBXX_L2_2.fq.gz Raw_Reads/HV_001_13_reverse.fq.gz
mv EG_Sflab_241_14_USPD16092227-14_HV5CHBBXX_L2_2.fq.gz Raw_Reads/HV_001_14_reverse.fq.gz
mv EG_Sflab_241_15_USPD16092227-15_HV5CHBBXX_L2_2.fq.gz Raw_Reads/HV_001_15_reverse.fq.gz
mv EG_Sflab_241_16_USPD16092227-16_HV5CHBBXX_L2_2.fq.gz Raw_Reads/HV_001_16_reverse.fq.gz
mv EG_Sflab_241_18_USPD16092227-18_HV5CHBBXX_L2_2.fq.gz Raw_Reads/HV_001_18_reverse.fq.gz
mv EG_Sflab_241_19_USPD16092227-19_HV5CHBBXX_L2_2.fq.gz Raw_Reads/HV_001_19_reverse.fq.gz
mv EG_Sflab_241_20_USPD16092227-20_HV5CHBBXX_L2_2.fq.gz Raw_Reads/HV_001_20_reverse.fq.gz
mv EG_Sflab_241_21_USPD16092227-21_HV5CHBBXX_L2_2.fq.gz Raw_Reads/HV_001_21_reverse.fq.gz
mv EG_Sflab_241_22_USPD16092227-22_HV5CHBBXX_L2_2.fq.gz Raw_Reads/HV_001_22_reverse.fq.gz
mv EG_Sflab_241_23_USPD16092227-23_HV5CHBBXX_L2_2.fq.gz Raw_Reads/HV_001_23_reverse.fq.gz
mv EG_Sflab_241_25_USPD16092227-25_HV5CHBBXX_L2_2.fq.gz Raw_Reads/HV_001_25_reverse.fq.gz
mv EG_Sflab_241_27_USPD16092227-27_HV5CHBBXX_L2_2.fq.gz Raw_Reads/HV_001_27_reverse.fq.gz
mv HV002_1_USPD16094887-1_HWJJHBBXX_L2_2.fq.gz Raw_Reads/HV_002_01_reverse.fq.gz
mv HV002_2_USPD16094887-2_HWJJHBBXX_L2_2.fq.gz Raw_Reads/HV_002_02_reverse.fq.gz
mv HV002_3_USPD16094887-3_HWJJHBBXX_L2_2.fq.gz Raw_Reads/HV_002_03_reverse.fq.gz
mv HV002_4_USPD16094887-4_HWJJHBBXX_L2_2.fq.gz Raw_Reads/HV_002_04_reverse.fq.gz
mv HV002_5_USPD16094887-5_HWJJHBBXX_L2_2.fq.gz Raw_Reads/HV_002_05_reverse.fq.gz
mv HV002_6_USPD16094887-6_HWJJHBBXX_L2_2.fq.gz Raw_Reads/HV_002_06_reverse.fq.gz
mv HV002_7_USPD16094887-7_HWJJHBBXX_L2_2.fq.gz Raw_Reads/HV_002_07_reverse.fq.gz
mv HV002_8_USPD16094887-8_HWJJHBBXX_L2_2.fq.gz Raw_Reads/HV_002_08_reverse.fq.gz
mv HV002_9_USPD16094887-9_HWJJHBBXX_L2_2.fq.gz Raw_Reads/HV_002_09_reverse.fq.gz
mv HV002_10_USPD16094887-10_HWJJHBBXX_L2_2.fq.gz Raw_Reads/HV_002_10_reverse.fq.gz
mv HV002_11_USPD16094887-11_HWJJHBBXX_L2_2.fq.gz Raw_Reads/HV_002_11_reverse.fq.gz
mv HV002_12_USPD16094887-12_HWJJHBBXX_L2_2.fq.gz Raw_Reads/HV_002_12_reverse.fq.gz
mv HV002_13_USPD16094887-13_HWJJHBBXX_L2_2.fq.gz Raw_Reads/HV_002_13_reverse.fq.gz
mv HV002_14_USPD16094887-14_HWJJHBBXX_L2_2.fq.gz Raw_Reads/HV_002_14_reverse.fq.gz
mv HV002_15_USPD16094887-15_HWJJHBBXX_L2_2.fq.gz Raw_Reads/HV_002_15_reverse.fq.gz
mv HV002_18_USPD16094887-18_HWJJHBBXX_L2_2.fq.gz Raw_Reads/HV_002_18_reverse.fq.gz
mv HV002_19_USPD16094887-19_HWJJHBBXX_L2_2.fq.gz Raw_Reads/HV_002_19_reverse.fq.gz
mv HV002_20_USPD16094887-20_HWJJHBBXX_L2_2.fq.gz Raw_Reads/HV_002_20_reverse.fq.gz
mv HV002_21_USPD16094887-21_HWJJHBBXX_L2_2.fq.gz Raw_Reads/HV_002_21_reverse.fq.gz
mv HV002_22_USPD16094887-22_HWJJHBBXX_L2_2.fq.gz Raw_Reads/HV_002_22_reverse.fq.gz
mv HV002_23_USPD16094887-23_HWJJHBBXX_L2_2.fq.gz Raw_Reads/HV_002_23_reverse.fq.gz
mv HV003_1_USPD16094888-1_HWJJHBBXX_L1_2.fq.gz Raw_Reads/HV_003_01_reverse.fq.gz
mv HV003_2_USPD16094888-2_HWJJHBBXX_L1_2.fq.gz Raw_Reads/HV_003_02_reverse.fq.gz
mv HV003_3_USPD16094888-3_HWJJHBBXX_L1_2.fq.gz Raw_Reads/HV_003_03_reverse.fq.gz
mv HV003_4_USPD16094888-4_HWJJHBBXX_L1_2.fq.gz Raw_Reads/HV_003_04_reverse.fq.gz
mv HV003_5_USPD16094888-5_HWJJHBBXX_L1_2.fq.gz Raw_Reads/HV_003_05_reverse.fq.gz
mv HV003_6_USPD16094888-6_HWJJHBBXX_L1_2.fq.gz Raw_Reads/HV_003_06_reverse.fq.gz
mv HV003_7_USPD16094888-7_HWJJHBBXX_L1_2.fq.gz Raw_Reads/HV_003_07_reverse.fq.gz
mv HV003_8_USPD16094888-8_HWJJHBBXX_L1_2.fq.gz Raw_Reads/HV_003_08_reverse.fq.gz
mv HV003_9_USPD16094888-9_HWJJHBBXX_L1_2.fq.gz Raw_Reads/HV_003_09_reverse.fq.gz
mv HV003_10_USPD16094888-10_HWJJHBBXX_L1_2.fq.gz Raw_Reads/HV_003_10_reverse.fq.gz
mv HV003_11_USPD16094888-11_HWJJHBBXX_L1_2.fq.gz Raw_Reads/HV_003_11_reverse.fq.gz
mv HV003_12_USPD16094888-12_HWJJHBBXX_L1_2.fq.gz Raw_Reads/HV_003_12_reverse.fq.gz
mv HV003_13_USPD16094888-13_HWJJHBBXX_L1_2.fq.gz Raw_Reads/HV_003_13_reverse.fq.gz
mv HV003_14_USPD16094888-14_HWJJHBBXX_L1_2.fq.gz Raw_Reads/HV_003_14_reverse.fq.gz
mv HV003_15_USPD16094888-15_HWJJHBBXX_L1_2.fq.gz Raw_Reads/HV_003_15_reverse.fq.gz
mv HV003_16_USPD16094888-16_HWJJHBBXX_L1_2.fq.gz Raw_Reads/HV_003_16_reverse.fq.gz
mv HV003_18_USPD16094888-18_HWJJHBBXX_L1_2.fq.gz Raw_Reads/HV_003_18_reverse.fq.gz
mv HV003_19_USPD16094888-19_HWJJHBBXX_L1_2.fq.gz Raw_Reads/HV_003_19_reverse.fq.gz
mv HV003_20_USPD16094888-20_HWJJHBBXX_L1_2.fq.gz Raw_Reads/HV_003_20_reverse.fq.gz
mv HV003_21_USPD16094888-21_HWJJHBBXX_L1_2.fq.gz Raw_Reads/HV_003_21_reverse.fq.gz
mv HV003_22_USPD16094888-22_HWJJHBBXX_L1_2.fq.gz Raw_Reads/HV_003_22_reverse.fq.gz
mv HV004_1_USPD16094889-1_HWKHHBBXX_L8_2.fq.gz Raw_Reads/HV_004_01_reverse.fq.gz
mv HV004_2_USPD16094889-2_HWKHHBBXX_L8_2.fq.gz Raw_Reads/HV_004_02_reverse.fq.gz
mv HV004_3_USPD16094889-3_HWKHHBBXX_L8_2.fq.gz Raw_Reads/HV_004_03_reverse.fq.gz
mv HV004_4_USPD16094889-4_HWKHHBBXX_L8_2.fq.gz Raw_Reads/HV_004_04_reverse.fq.gz
mv HV004_5_USPD16094889-5_HWKHHBBXX_L8_2.fq.gz Raw_Reads/HV_004_05_reverse.fq.gz
mv HV004_6_USPD16094889-6_HWKHHBBXX_L8_2.fq.gz Raw_Reads/HV_004_06_reverse.fq.gz
mv HV004_7_USPD16094889-7_HWKHHBBXX_L8_2.fq.gz Raw_Reads/HV_004_07_reverse.fq.gz
mv HV004_8_USPD16094889-8_HWKHHBBXX_L8_2.fq.gz Raw_Reads/HV_004_08_reverse.fq.gz
mv HV004_9_USPD16094889-9_HWKHHBBXX_L8_2.fq.gz Raw_Reads/HV_004_09_reverse.fq.gz
mv HV004_10_USPD16094889-10_HWKHHBBXX_L8_2.fq.gz Raw_Reads/HV_004_10_reverse.fq.gz
mv HV004_11_USPD16094889-11_HWKHHBBXX_L8_2.fq.gz Raw_Reads/HV_004_11_reverse.fq.gz
mv HV004_12_USPD16094889-12_HWKHHBBXX_L8_2.fq.gz Raw_Reads/HV_004_12_reverse.fq.gz
mv HV004_13_USPD16094889-13_HWKHHBBXX_L8_2.fq.gz Raw_Reads/HV_004_13_reverse.fq.gz
mv HV004_14_USPD16094889-14_HWKHHBBXX_L8_2.fq.gz Raw_Reads/HV_004_14_reverse.fq.gz
mv HV004_15_USPD16094889-15_HWKHHBBXX_L8_2.fq.gz Raw_Reads/HV_004_15_reverse.fq.gz
mv HV004_16_USPD16094889-16_HWKHHBBXX_L8_2.fq.gz Raw_Reads/HV_004_16_reverse.fq.gz
mv HV004_18_USPD16094889-18_HWKHHBBXX_L8_2.fq.gz Raw_Reads/HV_004_18_reverse.fq.gz
mv HV004_19_USPD16094889-19_HWKHHBBXX_L8_2.fq.gz Raw_Reads/HV_004_19_reverse.fq.gz
mv HV004_20_USPD16094889-20_HWKHHBBXX_L8_2.fq.gz Raw_Reads/HV_004_20_reverse.fq.gz
mv HV004_21_USPD16094889-21_HWKHHBBXX_L8_2.fq.gz Raw_Reads/HV_004_21_reverse.fq.gz
mv HV004_22_USPD16094889-22_HWKHHBBXX_L8_2.fq.gz Raw_Reads/HV_004_22_reverse.fq.gz
mv HV004_23_USPD16094889-23_HWKHHBBXX_L8_2.fq.gz Raw_Reads/HV_004_23_reverse.fq.gz
mv HV004_25_USPD16094889-25_HWKHHBBXX_L8_2.fq.gz Raw_Reads/HV_004_25_reverse.fq.gz
mv HV004_27_USPD16094889-27_HWKHHBBXX_L8_2.fq.gz Raw_Reads/HV_004_27_reverse.fq.gz
mv HV005_1_USPD16094891-1_HWKHHBBXX_L7_2.fq.gz Raw_Reads/HV_005_01_reverse.fq.gz
mv HV005_3_USPD16094891-3_HWKHHBBXX_L7_2.fq.gz Raw_Reads/HV_005_03_reverse.fq.gz
mv HV005_4_USPD16094891-4_HWKHHBBXX_L7_2.fq.gz Raw_Reads/HV_005_04_reverse.fq.gz
mv HV005_5_USPD16094891-5_HWKHHBBXX_L7_2.fq.gz Raw_Reads/HV_005_05_reverse.fq.gz
mv HV005_6_USPD16094891-6_HWKHHBBXX_L7_2.fq.gz Raw_Reads/HV_005_06_reverse.fq.gz
mv HV005_7_USPD16094891-7_HWKHHBBXX_L7_2.fq.gz Raw_Reads/HV_005_07_reverse.fq.gz
mv HV005_8_USPD16094891-8_HWKHHBBXX_L7_2.fq.gz Raw_Reads/HV_005_08_reverse.fq.gz
mv HV005_9_USPD16094891-9_HWKHHBBXX_L7_2.fq.gz Raw_Reads/HV_005_09_reverse.fq.gz
mv HV005_10_USPD16094891-10_HWKHHBBXX_L7_2.fq.gz Raw_Reads/HV_005_10_reverse.fq.gz
mv HV005_11_USPD16094891-11_HWKHHBBXX_L7_2.fq.gz Raw_Reads/HV_005_11_reverse.fq.gz
mv HV005_12_USPD16094891-12_HWKHHBBXX_L7_2.fq.gz Raw_Reads/HV_005_12_reverse.fq.gz
mv HV005_13_USPD16094891-13_HWKHHBBXX_L7_2.fq.gz Raw_Reads/HV_005_13_reverse.fq.gz
mv HV005_14_USPD16094891-14_HWKHHBBXX_L7_2.fq.gz Raw_Reads/HV_005_14_reverse.fq.gz
mv HV005_15_USPD16094891-15_HWKHHBBXX_L7_2.fq.gz Raw_Reads/HV_005_15_reverse.fq.gz
mv HV005_16_USPD16094891-16_HWKHHBBXX_L7_2.fq.gz Raw_Reads/HV_005_16_reverse.fq.gz
mv HV005_18_USPD16094891-18_HWKHHBBXX_L7_2.fq.gz Raw_Reads/HV_005_18_reverse.fq.gz
mv HV005_19_USPD16094891-19_HWKHHBBXX_L7_2.fq.gz Raw_Reads/HV_005_19_reverse.fq.gz
mv HV005_20_USPD16094891-20_HWKHHBBXX_L7_2.fq.gz Raw_Reads/HV_005_20_reverse.fq.gz
mv HV005_22_USPD16094891-22_HWKHHBBXX_L7_2.fq.gz Raw_Reads/HV_005_22_reverse.fq.gz
mv HV005_23_USPD16094891-23_HWKHHBBXX_L7_2.fq.gz Raw_Reads/HV_005_23_reverse.fq.gz
mv HV005_25_USPD16094891-25_HWKHHBBXX_L7_2.fq.gz Raw_Reads/HV_005_25_reverse.fq.gz
mv HV005_27_USPD16094891-27_HWKHHBBXX_L7_2.fq.gz Raw_Reads/HV_005_27_reverse.fq.gz
mv HV_006_1_USPD16099387-1_H2FJ5BBXX_L6_2.fq.gz Raw_Reads/HV_006_01_reverse.fq.gz
mv HV_006_2_USPD16099387-2_H2FJ5BBXX_L6_2.fq.gz Raw_Reads/HV_006_02_reverse.fq.gz
mv HV_006_3_USPD16099387-3_H2FJ5BBXX_L6_2.fq.gz Raw_Reads/HV_006_03_reverse.fq.gz
mv HV_006_4_USPD16099387-4_H2FJ5BBXX_L6_2.fq.gz Raw_Reads/HV_006_04_reverse.fq.gz
mv HV_006_5_USPD16099387-5_H2FJ5BBXX_L6_2.fq.gz Raw_Reads/HV_006_05_reverse.fq.gz
mv HV_006_6_USPD16099387-6_H2FJ5BBXX_L6_2.fq.gz Raw_Reads/HV_006_06_reverse.fq.gz
mv HV_006_7_USPD16099387-7_H2FJ5BBXX_L6_2.fq.gz Raw_Reads/HV_006_07_reverse.fq.gz
mv HV_006_8_USPD16099387-8_H2FJ5BBXX_L6_2.fq.gz Raw_Reads/HV_006_08_reverse.fq.gz
mv HV_006_9_USPD16099387-9_H2FJ5BBXX_L6_2.fq.gz Raw_Reads/HV_006_09_reverse.fq.gz
mv HV_006_10_USPD16099387-10_H2FJ5BBXX_L6_2.fq.gz Raw_Reads/HV_006_10_reverse.fq.gz
mv HV_006_11_USPD16099387-11_H2FJ5BBXX_L6_2.fq.gz Raw_Reads/HV_006_11_reverse.fq.gz
mv HV_006_12_USPD16099387-12_H2FJ5BBXX_L6_2.fq.gz Raw_Reads/HV_006_12_reverse.fq.gz
mv HV_006_13_USPD16099387-13_H2FJ5BBXX_L6_2.fq.gz Raw_Reads/HV_006_13_reverse.fq.gz
mv HV_006_14_USPD16099387-14_H2FJ5BBXX_L6_2.fq.gz Raw_Reads/HV_006_14_reverse.fq.gz
mv HV_006_15_USPD16099387-15_H2FJ5BBXX_L6_2.fq.gz Raw_Reads/HV_006_15_reverse.fq.gz
mv HV_006_16_USPD16099387-16_H2FJ5BBXX_L6_2.fq.gz Raw_Reads/HV_006_16_reverse.fq.gz
mv HV_006_18_USPD16099387-18_H2FJ5BBXX_L6_2.fq.gz Raw_Reads/HV_006_18_reverse.fq.gz
mv HV_006_19_USPD16099387-19_H2FJ5BBXX_L6_2.fq.gz Raw_Reads/HV_006_19_reverse.fq.gz
mv HV_006_20_USPD16099387-20_H2FJ5BBXX_L6_2.fq.gz Raw_Reads/HV_006_20_reverse.fq.gz
mv HV_006_21_USPD16099387-21_H2FJ5BBXX_L6_2.fq.gz Raw_Reads/HV_006_21_reverse.fq.gz
mv HV_006_22_USPD16099387-22_H2FJ5BBXX_L6_2.fq.gz Raw_Reads/HV_006_22_reverse.fq.gz
mv HV_006_23_USPD16099387-23_H2FJ5BBXX_L6_2.fq.gz Raw_Reads/HV_006_23_reverse.fq.gz
mv HV_006_25_USPD16099387-25_H2FJ5BBXX_L6_2.fq.gz Raw_Reads/HV_006_25_reverse.fq.gz
mv HV_006_27_USPD16099387-27_H2FJ5BBXX_L6_2.fq.gz Raw_Reads/HV_006_27_reverse.fq.gz
mv HV_007_1_USPD16099388-1_H2FJ5BBXX_L7_2.fq.gz Raw_Reads/HV_007_01_reverse.fq.gz
mv HV_007_2_USPD16099388-2_H2FJ5BBXX_L7_2.fq.gz Raw_Reads/HV_007_02_reverse.fq.gz
mv HV_007_3_USPD16099388-3_H2FJ5BBXX_L7_2.fq.gz Raw_Reads/HV_007_03_reverse.fq.gz
mv HV_007_4_USPD16099388-4_H2FJ5BBXX_L7_2.fq.gz Raw_Reads/HV_007_04_reverse.fq.gz
mv HV_007_5_USPD16099388-5_H2FJ5BBXX_L7_2.fq.gz Raw_Reads/HV_007_05_reverse.fq.gz
mv HV_007_6_USPD16099388-6_H2FJ5BBXX_L7_2.fq.gz Raw_Reads/HV_007_06_reverse.fq.gz
mv HV_007_7_USPD16099388-7_H2FJ5BBXX_L7_2.fq.gz Raw_Reads/HV_007_07_reverse.fq.gz
mv HV_007_8_USPD16099388-8_H2FJ5BBXX_L7_2.fq.gz Raw_Reads/HV_007_08_reverse.fq.gz
mv HV_007_9_USPD16099388-9_H2FJ5BBXX_L7_2.fq.gz Raw_Reads/HV_007_09_reverse.fq.gz
mv HV_007_10_USPD16099388-10_H2FJ5BBXX_L7_2.fq.gz Raw_Reads/HV_007_10_reverse.fq.gz
mv HV_007_11_USPD16099388-11_H2FJ5BBXX_L7_2.fq.gz Raw_Reads/HV_007_11_reverse.fq.gz
mv HV_007_12_USPD16099388-12_H2FJ5BBXX_L7_2.fq.gz Raw_Reads/HV_007_12_reverse.fq.gz
mv HV_007_13_USPD16099388-13_H2FJ5BBXX_L7_2.fq.gz Raw_Reads/HV_007_13_reverse.fq.gz
mv HV_007_14_USPD16099388-14_H2FJ5BBXX_L7_2.fq.gz Raw_Reads/HV_007_14_reverse.fq.gz
mv HV_007_15_USPD16099388-15_H2FJ5BBXX_L7_2.fq.gz Raw_Reads/HV_007_15_reverse.fq.gz
mv HV_007_16_USPD16099388-16_H2FJ5BBXX_L7_2.fq.gz Raw_Reads/HV_007_16_reverse.fq.gz
mv HV_007_18_USPD16099388-18_H2FJ5BBXX_L7_2.fq.gz Raw_Reads/HV_007_18_reverse.fq.gz
mv HV_007_19_USPD16099388-19_H2FJ5BBXX_L7_2.fq.gz Raw_Reads/HV_007_19_reverse.fq.gz
mv HV_007_20_USPD16099388-20_H2FJ5BBXX_L7_2.fq.gz Raw_Reads/HV_007_20_reverse.fq.gz
mv HV_007_21_USPD16099388-21_H2FJ5BBXX_L7_2.fq.gz Raw_Reads/HV_007_21_reverse.fq.gz
mv HV_007_22_USPD16099388-22_H2FJ5BBXX_L7_2.fq.gz Raw_Reads/HV_007_22_reverse.fq.gz
mv HV_007_23_USPD16099388-23_H2FJ5BBXX_L7_2.fq.gz Raw_Reads/HV_007_23_reverse.fq.gz
mv HV_007_25_USPD16099388-25_H2FJ5BBXX_L7_2.fq.gz Raw_Reads/HV_007_25_reverse.fq.gz
mv HV_007_27_USPD16099388-27_H2FJ5BBXX_L7_2.fq.gz Raw_Reads/HV_007_27_reverse.fq.gz
mv HV_008_1_USPD16099385-1_H2FJ5BBXX_L4_2.fq.gz Raw_Reads/HV_008_01_reverse.fq.gz
mv HV_008_2_USPD16099385-2_H2FJ5BBXX_L4_2.fq.gz Raw_Reads/HV_008_02_reverse.fq.gz
mv HV_008_3_USPD16099385-3_H2FJ5BBXX_L4_2.fq.gz Raw_Reads/HV_008_03_reverse.fq.gz
mv HV_008_4_USPD16099385-4_H2FJ5BBXX_L4_2.fq.gz Raw_Reads/HV_008_04_reverse.fq.gz
mv HV_008_5_USPD16099385-5_H2FJ5BBXX_L4_2.fq.gz Raw_Reads/HV_008_05_reverse.fq.gz
mv HV_008_6_USPD16099385-6_H2FJ5BBXX_L4_2.fq.gz Raw_Reads/HV_008_06_reverse.fq.gz
mv HV_008_7_USPD16099385-7_H2FJ5BBXX_L4_2.fq.gz Raw_Reads/HV_008_07_reverse.fq.gz
mv HV_008_8_USPD16099385-8_H2FJ5BBXX_L4_2.fq.gz Raw_Reads/HV_008_08_reverse.fq.gz
mv HV_008_9_USPD16099385-9_H2FJ5BBXX_L4_2.fq.gz Raw_Reads/HV_008_09_reverse.fq.gz
mv HV_008_10_USPD16099385-10_H2FJ5BBXX_L4_2.fq.gz Raw_Reads/HV_008_10_reverse.fq.gz
mv HV_008_11_USPD16099385-11_H2FJ5BBXX_L4_2.fq.gz Raw_Reads/HV_008_11_reverse.fq.gz
mv HV_008_12_USPD16099385-12_H2FJ5BBXX_L4_2.fq.gz Raw_Reads/HV_008_12_reverse.fq.gz
mv HV_008_13_USPD16099385-13_H2FJ5BBXX_L4_2.fq.gz Raw_Reads/HV_008_13_reverse.fq.gz
mv HV_008_14_USPD16099385-14_H2FJ5BBXX_L4_2.fq.gz Raw_Reads/HV_008_14_reverse.fq.gz
mv HV_008_15_USPD16099385-15_H2FJ5BBXX_L4_2.fq.gz Raw_Reads/HV_008_15_reverse.fq.gz
mv HV_008_16_USPD16099385-16_H2FJ5BBXX_L4_2.fq.gz Raw_Reads/HV_008_16_reverse.fq.gz
mv HV_008_18_USPD16099385-18_H2FJ5BBXX_L4_2.fq.gz Raw_Reads/HV_008_18_reverse.fq.gz
mv HV_008_19_USPD16099385-19_H2FJ5BBXX_L4_2.fq.gz Raw_Reads/HV_008_19_reverse.fq.gz
mv HV_008_20_USPD16099385-20_H2FJ5BBXX_L4_2.fq.gz Raw_Reads/HV_008_20_reverse.fq.gz
mv HV_008_21_USPD16099385-21_H2FJ5BBXX_L4_2.fq.gz Raw_Reads/HV_008_21_reverse.fq.gz
mv HV_008_22_USPD16099385-22_H2FJ5BBXX_L4_2.fq.gz Raw_Reads/HV_008_22_reverse.fq.gz
mv HV_008_23_USPD16099385-23_H2FJ5BBXX_L4_2.fq.gz Raw_Reads/HV_008_23_reverse.fq.gz
mv HV_008_25_USPD16099385-25_H2FJ5BBXX_L4_2.fq.gz Raw_Reads/HV_008_25_reverse.fq.gz
mv HV_008_27_USPD16099385-27_H2FJ5BBXX_L4_2.fq.gz Raw_Reads/HV_008_27_reverse.fq.gz
mv HV_009_1_USPD16099386-1_H2FJ5BBXX_L5_2.fq.gz Raw_Reads/HV_009_01_reverse.fq.gz
mv HV_009_3_USPD16099386-3_H2FJ5BBXX_L5_2.fq.gz Raw_Reads/HV_009_03_reverse.fq.gz
mv HV_009_4_USPD16099386-4_H2FJ5BBXX_L5_2.fq.gz Raw_Reads/HV_009_04_reverse.fq.gz
mv HV_009_5_USPD16099386-5_H2FJ5BBXX_L5_2.fq.gz Raw_Reads/HV_009_05_reverse.fq.gz
mv HV_009_6_USPD16099386-6_H2FJ5BBXX_L5_2.fq.gz Raw_Reads/HV_009_06_reverse.fq.gz
mv HV_009_7_USPD16099386-7_H2FJ5BBXX_L5_2.fq.gz Raw_Reads/HV_009_07_reverse.fq.gz
mv HV_009_8_USPD16099386-8_H2FJ5BBXX_L5_2.fq.gz Raw_Reads/HV_009_08_reverse.fq.gz
mv HV_009_9_USPD16099386-9_H2FJ5BBXX_L5_2.fq.gz Raw_Reads/HV_009_09_reverse.fq.gz
mv HV_009_10_USPD16099386-10_H2FJ5BBXX_L5_2.fq.gz Raw_Reads/HV_009_10_reverse.fq.gz
mv HV_009_11_USPD16099386-11_H2FJ5BBXX_L5_2.fq.gz Raw_Reads/HV_009_11_reverse.fq.gz
mv HV_009_12_USPD16099386-12_H2FJ5BBXX_L5_2.fq.gz Raw_Reads/HV_009_12_reverse.fq.gz
mv HV_009_13_USPD16099386-13_H2FJ5BBXX_L5_2.fq.gz Raw_Reads/HV_009_13_reverse.fq.gz
mv HV_009_14_USPD16099386-14_H2FJ5BBXX_L5_2.fq.gz Raw_Reads/HV_009_14_reverse.fq.gz
mv HV_009_15_USPD16099386-15_H2FJ5BBXX_L5_2.fq.gz Raw_Reads/HV_009_15_reverse.fq.gz
mv HV_009_16_USPD16099386-16_H2FJ5BBXX_L5_2.fq.gz Raw_Reads/HV_009_16_reverse.fq.gz
mv HV_009_18_USPD16099386-18_H2FJ5BBXX_L5_2.fq.gz Raw_Reads/HV_009_18_reverse.fq.gz
mv HV_009_19_USPD16099386-19_H2FJ5BBXX_L5_2.fq.gz Raw_Reads/HV_009_19_reverse.fq.gz
mv HV_009_20_USPD16099386-20_H2FJ5BBXX_L5_2.fq.gz Raw_Reads/HV_009_20_reverse.fq.gz
mv HV_009_21_USPD16099386-21_H2FJ5BBXX_L5_2.fq.gz Raw_Reads/HV_009_21_reverse.fq.gz
mv HV_009_22_USPD16099386-22_H2FJ5BBXX_L5_2.fq.gz Raw_Reads/HV_009_22_reverse.fq.gz
mv HV_009_23_USPD16099386-23_H2FJ5BBXX_L5_2.fq.gz Raw_Reads/HV_009_23_reverse.fq.gz
mv HV_009_25_USPD16099386-25_H2FJ5BBXX_L5_2.fq.gz Raw_Reads/HV_009_25_reverse.fq.gz
mv HV_009_27_USPD16099386-27_H2FJ5BBXX_L5_2.fq.gz Raw_Reads/HV_009_27_reverse.fq.gz
mv HV_010_1_USPD16099384-1_H2FJ5BBXX_L3_2.fq.gz Raw_Reads/HV_010_01_reverse.fq.gz
mv HV_010_2_USPD16099384-2_H2FJ5BBXX_L3_2.fq.gz Raw_Reads/HV_010_02_reverse.fq.gz
mv HV_010_3_USPD16099384-3_H2FJ5BBXX_L3_2.fq.gz Raw_Reads/HV_010_03_reverse.fq.gz
mv HV_010_4_USPD16099384-4_H2FJ5BBXX_L3_2.fq.gz Raw_Reads/HV_010_04_reverse.fq.gz
mv HV_010_5_USPD16099384-5_H2FJ5BBXX_L3_2.fq.gz Raw_Reads/HV_010_05_reverse.fq.gz
mv HV_010_6_USPD16099384-6_H2FJ5BBXX_L3_2.fq.gz Raw_Reads/HV_010_06_reverse.fq.gz
mv HV_010_7_USPD16099384-7_H2FJ5BBXX_L3_2.fq.gz Raw_Reads/HV_010_07_reverse.fq.gz
mv HV_010_8_USPD16099384-8_H2FJ5BBXX_L3_2.fq.gz Raw_Reads/HV_010_08_reverse.fq.gz
mv HV_010_9_USPD16099384-9_H2FJ5BBXX_L3_2.fq.gz Raw_Reads/HV_010_09_reverse.fq.gz
mv HV_010_10_USPD16099384-10_H2FJ5BBXX_L3_2.fq.gz Raw_Reads/HV_010_10_reverse.fq.gz
mv HV_010_11_USPD16099384-11_H2FJ5BBXX_L3_2.fq.gz Raw_Reads/HV_010_11_reverse.fq.gz
mv HV_010_12_USPD16099384-12_H2FJ5BBXX_L3_2.fq.gz Raw_Reads/HV_010_12_reverse.fq.gz
mv HV_010_13_USPD16099384-13_H2FJ5BBXX_L3_2.fq.gz Raw_Reads/HV_010_13_reverse.fq.gz
mv HV_010_14_USPD16099384-14_H2FJ5BBXX_L3_2.fq.gz Raw_Reads/HV_010_14_reverse.fq.gz
mv HV_010_15_USPD16099384-15_H2FJ5BBXX_L3_2.fq.gz Raw_Reads/HV_010_15_reverse.fq.gz
mv HV_010_16_USPD16099384-16_H2FJ5BBXX_L3_2.fq.gz Raw_Reads/HV_010_16_reverse.fq.gz
mv HV_010_18_USPD16099384-18_H2FJ5BBXX_L3_2.fq.gz Raw_Reads/HV_010_18_reverse.fq.gz
mv HV_010_19_USPD16099384-19_H2FJ5BBXX_L3_2.fq.gz Raw_Reads/HV_010_19_reverse.fq.gz
mv HV_010_20_USPD16099384-20_H2FJ5BBXX_L3_2.fq.gz Raw_Reads/HV_010_20_reverse.fq.gz
mv HV_010_21_USPD16099384-21_H2FJ5BBXX_L3_2.fq.gz Raw_Reads/HV_010_21_reverse.fq.gz
mv HV_010_22_USPD16099384-22_H2FJ5BBXX_L3_2.fq.gz Raw_Reads/HV_010_22_reverse.fq.gz
mv HV_010_25_USPD16099384-25_H2FJ5BBXX_L3_2.fq.gz Raw_Reads/HV_010_25_reverse.fq.gz
mv HV_010_27_USPD16099384-27_H2FJ5BBXX_L3_2.fq.gz Raw_Reads/HV_010_27_reverse.fq.gz
mv HV_011_1_USPD16099389-1_H2FJ5BBXX_L8_2.fq.gz Raw_Reads/HV_011_01_reverse.fq.gz
mv HV_011_2_USPD16099389-2_H2FJ5BBXX_L8_2.fq.gz Raw_Reads/HV_011_02_reverse.fq.gz
mv HV_011_3_USPD16099389-3_H2FJ5BBXX_L8_2.fq.gz Raw_Reads/HV_011_03_reverse.fq.gz
mv HV_011_4_USPD16099389-4_H2FJ5BBXX_L8_2.fq.gz Raw_Reads/HV_011_04_reverse.fq.gz
mv HV_011_5_USPD16099389-5_H2FJ5BBXX_L8_2.fq.gz Raw_Reads/HV_011_05_reverse.fq.gz
mv HV_011_6_USPD16099389-6_H2FJ5BBXX_L8_2.fq.gz Raw_Reads/HV_011_06_reverse.fq.gz
mv HV_011_7_USPD16099389-7_H2FJ5BBXX_L8_2.fq.gz Raw_Reads/HV_011_07_reverse.fq.gz
mv HV_011_8_USPD16099389-8_H2FJ5BBXX_L8_2.fq.gz Raw_Reads/HV_011_08_reverse.fq.gz
mv HV_011_9_USPD16099389-9_H2FJ5BBXX_L8_2.fq.gz Raw_Reads/HV_011_09_reverse.fq.gz
mv HV_011_10_USPD16099389-10_H2FJ5BBXX_L8_2.fq.gz Raw_Reads/HV_011_10_reverse.fq.gz
mv HV_011_11_USPD16099389-11_H2FJ5BBXX_L8_2.fq.gz Raw_Reads/HV_011_11_reverse.fq.gz
mv HV_011_12_USPD16099389-12_H2FJ5BBXX_L8_2.fq.gz Raw_Reads/HV_011_12_reverse.fq.gz
mv HV_011_13_USPD16099389-13_H2FJ5BBXX_L8_2.fq.gz Raw_Reads/HV_011_13_reverse.fq.gz
mv HV_011_14_USPD16099389-14_H2FJ5BBXX_L8_2.fq.gz Raw_Reads/HV_011_14_reverse.fq.gz
mv HV_011_15_USPD16099389-15_H2FJ5BBXX_L8_2.fq.gz Raw_Reads/HV_011_15_reverse.fq.gz
mv HV_011_16_USPD16099389-16_H2FJ5BBXX_L8_2.fq.gz Raw_Reads/HV_011_16_reverse.fq.gz
mv HV_011_18_USPD16099389-18_H2FJ5BBXX_L8_2.fq.gz Raw_Reads/HV_011_18_reverse.fq.gz
mv HV_011_19_USPD16099389-19_H2FJ5BBXX_L8_2.fq.gz Raw_Reads/HV_011_19_reverse.fq.gz
mv HV_011_20_USPD16099389-20_H2FJ5BBXX_L8_2.fq.gz Raw_Reads/HV_011_20_reverse.fq.gz
mv HV_011_21_USPD16099389-21_H2FJ5BBXX_L8_2.fq.gz Raw_Reads/HV_011_21_reverse.fq.gz
mv HV_011_22_USPD16099389-22_H2FJ5BBXX_L8_2.fq.gz Raw_Reads/HV_011_22_reverse.fq.gz
mv HV_011_23_USPD16099389-23_H2FJ5BBXX_L8_2.fq.gz Raw_Reads/HV_011_23_reverse.fq.gz
mv HV_011_25_USPD16099389-25_H2FJ5BBXX_L8_2.fq.gz Raw_Reads/HV_011_25_reverse.fq.gz
mv HV_011_27_USPD16099389-27_H2FJ5BBXX_L8_2.fq.gz Raw_Reads/HV_011_27_reverse.fq.gz
mv HV012_1_USPD16102686-1_H333VBBXX_L2_2.fq.gz Raw_Reads/HV_012_01_reverse.fq.gz
mv HV012_2_USPD16102686-2_H333VBBXX_L2_2.fq.gz Raw_Reads/HV_012_02_reverse.fq.gz
mv HV012_3_USPD16102686-3_H333VBBXX_L2_2.fq.gz Raw_Reads/HV_012_03_reverse.fq.gz
mv HV012_4_USPD16102686-4_H333VBBXX_L2_2.fq.gz Raw_Reads/HV_012_04_reverse.fq.gz
mv HV012_5_USPD16102686-5_H333VBBXX_L2_2.fq.gz Raw_Reads/HV_012_05_reverse.fq.gz
mv HV012_6_USPD16102686-6_H333VBBXX_L2_2.fq.gz Raw_Reads/HV_012_06_reverse.fq.gz
mv HV012_7_USPD16102686-7_H333VBBXX_L2_2.fq.gz Raw_Reads/HV_012_07_reverse.fq.gz
mv HV012_8_USPD16102686-8_H333VBBXX_L2_2.fq.gz Raw_Reads/HV_012_08_reverse.fq.gz
mv HV012_9_USPD16102686-9_H333VBBXX_L2_2.fq.gz Raw_Reads/HV_012_09_reverse.fq.gz
mv HV012_10_USPD16102686-10_H333VBBXX_L2_2.fq.gz Raw_Reads/HV_012_10_reverse.fq.gz
mv HV012_11_USPD16102686-11_H333VBBXX_L2_2.fq.gz Raw_Reads/HV_012_11_reverse.fq.gz
mv HV012_12_USPD16102686-12_H333VBBXX_L2_2.fq.gz Raw_Reads/HV_012_12_reverse.fq.gz
mv HV012_13_USPD16102686-13_H333VBBXX_L2_2.fq.gz Raw_Reads/HV_012_13_reverse.fq.gz
mv HV012_14_USPD16102686-14_H333VBBXX_L2_2.fq.gz Raw_Reads/HV_012_14_reverse.fq.gz
mv HV012_15_USPD16102686-15_H333VBBXX_L2_2.fq.gz Raw_Reads/HV_012_15_reverse.fq.gz
mv HV012_16_USPD16102686-16_H333VBBXX_L2_2.fq.gz Raw_Reads/HV_012_16_reverse.fq.gz
mv HV012_18_USPD16102686-18_H333VBBXX_L2_2.fq.gz Raw_Reads/HV_012_18_reverse.fq.gz
mv HV012_19_USPD16102686-19_H333VBBXX_L2_2.fq.gz Raw_Reads/HV_012_19_reverse.fq.gz
mv HV012_20_USPD16102686-20_H333VBBXX_L2_2.fq.gz Raw_Reads/HV_012_20_reverse.fq.gz
mv HV012_21_USPD16102686-21_H333VBBXX_L2_2.fq.gz Raw_Reads/HV_012_21_reverse.fq.gz
mv HV012_22_USPD16102686-22_H333VBBXX_L2_2.fq.gz Raw_Reads/HV_012_22_reverse.fq.gz
mv HV012_23_USPD16102686-23_H333VBBXX_L2_2.fq.gz Raw_Reads/HV_012_23_reverse.fq.gz
mv HV012_25_USPD16102686-25_H333VBBXX_L2_2.fq.gz Raw_Reads/HV_012_25_reverse.fq.gz
mv HV012_27_USPD16102686-27_H333VBBXX_L2_2.fq.gz Raw_Reads/HV_012_27_reverse.fq.gz
mv HV013_1_USPD16102680-1_H3FVNBBXX_L1_2.fq.gz Raw_Reads/HV_013_01_reverse.fq.gz
mv HV013_2_USPD16102680-2_H3FVNBBXX_L1_2.fq.gz Raw_Reads/HV_013_02_reverse.fq.gz
mv HV013_3_USPD16102680-3_H3FVNBBXX_L1_2.fq.gz Raw_Reads/HV_013_03_reverse.fq.gz
mv HV013_4_USPD16102680-4_H3FVNBBXX_L1_2.fq.gz Raw_Reads/HV_013_04_reverse.fq.gz
mv HV013_5_USPD16102680-5_H3FVNBBXX_L1_2.fq.gz Raw_Reads/HV_013_05_reverse.fq.gz
mv HV013_7_USPD16102680-7_H3FVNBBXX_L1_2.fq.gz Raw_Reads/HV_013_07_reverse.fq.gz
mv HV013_8_USPD16102680-8_H3FVNBBXX_L1_2.fq.gz Raw_Reads/HV_013_08_reverse.fq.gz
mv HV013_9_USPD16102680-9_H3FVNBBXX_L1_2.fq.gz Raw_Reads/HV_013_09_reverse.fq.gz
mv HV013_10_USPD16102680-10_H3FVNBBXX_L1_2.fq.gz Raw_Reads/HV_013_10_reverse.fq.gz
mv HV013_11_USPD16102680-11_H3FVNBBXX_L1_2.fq.gz Raw_Reads/HV_013_11_reverse.fq.gz
mv HV013_12_USPD16102680-12_H3FVNBBXX_L1_2.fq.gz Raw_Reads/HV_013_12_reverse.fq.gz
mv HV013_13_USPD16102680-13_H3FVNBBXX_L1_2.fq.gz Raw_Reads/HV_013_13_reverse.fq.gz
mv HV013_14_USPD16102680-14_H3FVNBBXX_L1_2.fq.gz Raw_Reads/HV_013_14_reverse.fq.gz
mv HV013_15_USPD16102680-15_H3FVNBBXX_L1_2.fq.gz Raw_Reads/HV_013_15_reverse.fq.gz
mv HV013_16_USPD16102680-16_H3FVNBBXX_L1_2.fq.gz Raw_Reads/HV_013_16_reverse.fq.gz
mv HV013_18_USPD16102680-18_H3FVNBBXX_L1_2.fq.gz Raw_Reads/HV_013_18_reverse.fq.gz
mv HV013_19_USPD16102680-19_H3FVNBBXX_L1_2.fq.gz Raw_Reads/HV_013_19_reverse.fq.gz
mv HV013_20_USPD16102680-20_H3FVNBBXX_L1_2.fq.gz Raw_Reads/HV_013_20_reverse.fq.gz
mv HV013_21_USPD16102680-21_H3FVNBBXX_L1_2.fq.gz Raw_Reads/HV_013_21_reverse.fq.gz
mv HV013_22_USPD16102680-22_H3FVNBBXX_L1_2.fq.gz Raw_Reads/HV_013_22_reverse.fq.gz
mv HV013_23_USPD16102680-23_H3FVNBBXX_L1_2.fq.gz Raw_Reads/HV_013_23_reverse.fq.gz
mv HV013_25_USPD16102680-25_H3FVNBBXX_L1_2.fq.gz Raw_Reads/HV_013_25_reverse.fq.gz
mv HV013_27_USPD16102680-27_H3FVNBBXX_L1_2.fq.gz Raw_Reads/HV_013_27_reverse.fq.gz
mv HV014_1_USPD16102681-1_H3FVNBBXX_L2_2.fq.gz Raw_Reads/HV_014_01_reverse.fq.gz
mv HV014_2_USPD16102681-2_H3FVNBBXX_L2_2.fq.gz Raw_Reads/HV_014_02_reverse.fq.gz
mv HV014_3_USPD16102681-3_H3FVNBBXX_L2_2.fq.gz Raw_Reads/HV_014_03_reverse.fq.gz
mv HV014_4_USPD16102681-4_H3FVNBBXX_L2_2.fq.gz Raw_Reads/HV_014_04_reverse.fq.gz
mv HV014_5_USPD16102681-5_H3FVNBBXX_L2_2.fq.gz Raw_Reads/HV_014_05_reverse.fq.gz
mv HV014_6_USPD16102681-6_H3FVNBBXX_L2_2.fq.gz Raw_Reads/HV_014_06_reverse.fq.gz
mv HV014_7_USPD16102681-7_H3FVNBBXX_L2_2.fq.gz Raw_Reads/HV_014_07_reverse.fq.gz
mv HV014_8_USPD16102681-8_H3FVNBBXX_L2_2.fq.gz Raw_Reads/HV_014_08_reverse.fq.gz
mv HV014_9_USPD16102681-9_H3FVNBBXX_L2_2.fq.gz Raw_Reads/HV_014_09_reverse.fq.gz
mv HV014_10_USPD16102681-10_H3FVNBBXX_L2_2.fq.gz Raw_Reads/HV_014_10_reverse.fq.gz
mv HV014_11_USPD16102681-11_H3FVNBBXX_L2_2.fq.gz Raw_Reads/HV_014_11_reverse.fq.gz
mv HV014_12_USPD16102681-12_H3FVNBBXX_L2_2.fq.gz Raw_Reads/HV_014_12_reverse.fq.gz
mv HV014_13_USPD16102681-13_H3FVNBBXX_L2_2.fq.gz Raw_Reads/HV_014_13_reverse.fq.gz
mv HV014_14_USPD16102681-14_H3FVNBBXX_L2_2.fq.gz Raw_Reads/HV_014_14_reverse.fq.gz
mv HV014_15_USPD16102681-15_H3FVNBBXX_L2_2.fq.gz Raw_Reads/HV_014_15_reverse.fq.gz
mv HV014_16_USPD16102681-16_H3FVNBBXX_L2_2.fq.gz Raw_Reads/HV_014_16_reverse.fq.gz
mv HV014_18_USPD16102681-18_H3FVNBBXX_L2_2.fq.gz Raw_Reads/HV_014_18_reverse.fq.gz
mv HV014_19_USPD16102681-19_H3FVNBBXX_L2_2.fq.gz Raw_Reads/HV_014_19_reverse.fq.gz
mv HV014_20_USPD16102681-20_H3FVNBBXX_L2_2.fq.gz Raw_Reads/HV_014_20_reverse.fq.gz
mv HV014_21_USPD16102681-21_H3FVNBBXX_L2_2.fq.gz Raw_Reads/HV_014_21_reverse.fq.gz
mv HV014_22_USPD16102681-22_H3FVNBBXX_L2_2.fq.gz Raw_Reads/HV_014_22_reverse.fq.gz
mv HV014_23_USPD16102681-23_H3FVNBBXX_L2_2.fq.gz Raw_Reads/HV_014_23_reverse.fq.gz
mv HV014_25_USPD16102681-25_H3FVNBBXX_L2_2.fq.gz Raw_Reads/HV_014_25_reverse.fq.gz
mv HV014_27_USPD16102681-27_H3FVNBBXX_L2_2.fq.gz Raw_Reads/HV_014_27_reverse.fq.gz
mv HV015_1_USPD16102683-1_H3FVNBBXX_L4_2.fq.gz Raw_Reads/HV_015_01_reverse.fq.gz
mv HV015_2_USPD16102683-2_H3FVNBBXX_L4_2.fq.gz Raw_Reads/HV_015_02_reverse.fq.gz
mv HV015_3_USPD16102683-3_H3FVNBBXX_L4_2.fq.gz Raw_Reads/HV_015_03_reverse.fq.gz
mv HV015_4_USPD16102683-4_H3FVNBBXX_L4_2.fq.gz Raw_Reads/HV_015_04_reverse.fq.gz
mv HV015_5_USPD16102683-5_H3FVNBBXX_L4_2.fq.gz Raw_Reads/HV_015_05_reverse.fq.gz
mv HV015_6_USPD16102683-6_H3FVNBBXX_L4_2.fq.gz Raw_Reads/HV_015_06_reverse.fq.gz
mv HV015_7_USPD16102683-7_H3FVNBBXX_L4_2.fq.gz Raw_Reads/HV_015_07_reverse.fq.gz
mv HV015_8_USPD16102683-8_H3FVNBBXX_L4_2.fq.gz Raw_Reads/HV_015_08_reverse.fq.gz
mv HV015_9_USPD16102683-9_H3FVNBBXX_L4_2.fq.gz Raw_Reads/HV_015_09_reverse.fq.gz
mv HV015_10_USPD16102683-10_H3FVNBBXX_L4_2.fq.gz Raw_Reads/HV_015_10_reverse.fq.gz
mv HV015_11_USPD16102683-11_H3FVNBBXX_L4_2.fq.gz Raw_Reads/HV_015_11_reverse.fq.gz
mv HV015_12_USPD16102683-12_H3FVNBBXX_L4_2.fq.gz Raw_Reads/HV_015_12_reverse.fq.gz
mv HV015_13_USPD16102683-13_H3FVNBBXX_L4_2.fq.gz Raw_Reads/HV_015_13_reverse.fq.gz
mv HV015_14_USPD16102683-14_H3FVNBBXX_L4_2.fq.gz Raw_Reads/HV_015_14_reverse.fq.gz
mv HV015_15_USPD16102683-15_H3FVNBBXX_L4_2.fq.gz Raw_Reads/HV_015_15_reverse.fq.gz
mv HV015_16_USPD16102683-16_H3FVNBBXX_L4_2.fq.gz Raw_Reads/HV_015_16_reverse.fq.gz
mv HV015_18_USPD16102683-18_H3FVNBBXX_L4_2.fq.gz Raw_Reads/HV_015_18_reverse.fq.gz
mv HV015_19_USPD16102683-19_H3FVNBBXX_L4_2.fq.gz Raw_Reads/HV_015_19_reverse.fq.gz
mv HV015_20_USPD16102683-20_H3FVNBBXX_L4_2.fq.gz Raw_Reads/HV_015_20_reverse.fq.gz
mv HV015_21_USPD16102683-21_H3FVNBBXX_L4_2.fq.gz Raw_Reads/HV_015_21_reverse.fq.gz
mv HV015_23_USPD16102683-23_H3FVNBBXX_L4_2.fq.gz Raw_Reads/HV_015_23_reverse.fq.gz
mv HV015_25_USPD16102683-25_H3FVNBBXX_L4_2.fq.gz Raw_Reads/HV_015_25_reverse.fq.gz
mv HV015_27_USPD16102683-27_H3FVNBBXX_L4_2.fq.gz Raw_Reads/HV_015_27_reverse.fq.gz
mv HV016_1_USPD16102682-1_H3FVNBBXX_L3_2.fq.gz Raw_Reads/HV_016_01_reverse.fq.gz
mv HV016_2_USPD16102682-2_H3FVNBBXX_L3_2.fq.gz Raw_Reads/HV_016_02_reverse.fq.gz
mv HV016_3_USPD16102682-3_H3FVNBBXX_L3_2.fq.gz Raw_Reads/HV_016_03_reverse.fq.gz
mv HV016_4_USPD16102682-4_H3FVNBBXX_L3_2.fq.gz Raw_Reads/HV_016_04_reverse.fq.gz
mv HV016_5_USPD16102682-5_H3FVNBBXX_L3_2.fq.gz Raw_Reads/HV_016_05_reverse.fq.gz
mv HV016_6_USPD16102682-6_H3FVNBBXX_L3_2.fq.gz Raw_Reads/HV_016_06_reverse.fq.gz
mv HV016_7_USPD16102682-7_H3FVNBBXX_L3_2.fq.gz Raw_Reads/HV_016_07_reverse.fq.gz
mv HV016_8_USPD16102682-8_H3FVNBBXX_L3_2.fq.gz Raw_Reads/HV_016_08_reverse.fq.gz
mv HV016_9_USPD16102682-9_H3FVNBBXX_L3_2.fq.gz Raw_Reads/HV_016_09_reverse.fq.gz
mv HV016_10_USPD16102682-10_H3FVNBBXX_L3_2.fq.gz Raw_Reads/HV_016_10_reverse.fq.gz
mv HV016_11_USPD16102682-11_H3FVNBBXX_L3_2.fq.gz Raw_Reads/HV_016_11_reverse.fq.gz
mv HV016_12_USPD16102682-12_H3FVNBBXX_L3_2.fq.gz Raw_Reads/HV_016_12_reverse.fq.gz
mv HV016_13_USPD16102682-13_H3FVNBBXX_L3_2.fq.gz Raw_Reads/HV_016_13_reverse.fq.gz
mv HV016_14_USPD16102682-14_H3FVNBBXX_L3_2.fq.gz Raw_Reads/HV_016_14_reverse.fq.gz
mv HV016_15_USPD16102682-15_H3FVNBBXX_L3_2.fq.gz Raw_Reads/HV_016_15_reverse.fq.gz
mv HV016_16_USPD16102682-16_H3FVNBBXX_L3_2.fq.gz Raw_Reads/HV_016_16_reverse.fq.gz
mv HV016_18_USPD16102682-18_H3FVNBBXX_L3_2.fq.gz Raw_Reads/HV_016_18_reverse.fq.gz
mv HV016_19_USPD16102682-19_H3FVNBBXX_L3_2.fq.gz Raw_Reads/HV_016_19_reverse.fq.gz
mv HV016_20_USPD16102682-20_H3FVNBBXX_L3_2.fq.gz Raw_Reads/HV_016_20_reverse.fq.gz
mv HV016_21_USPD16102682-21_H3FVNBBXX_L3_2.fq.gz Raw_Reads/HV_016_21_reverse.fq.gz
mv HV016_22_USPD16102682-22_H3FVNBBXX_L3_2.fq.gz Raw_Reads/HV_016_22_reverse.fq.gz
mv HV016_23_USPD16102682-23_H3FVNBBXX_L3_2.fq.gz Raw_Reads/HV_016_23_reverse.fq.gz
mv HV016_25_USPD16102682-25_H3FVNBBXX_L3_2.fq.gz Raw_Reads/HV_016_25_reverse.fq.gz
mv HV016_27_USPD16102682-27_H3FVNBBXX_L3_2.fq.gz Raw_Reads/HV_016_27_reverse.fq.gz
mv HV017_1_USPD16102684-1_H3FVNBBXX_L5_2.fq.gz Raw_Reads/HV_017_01_reverse.fq.gz
mv HV017_2_USPD16102684-2_H3FVNBBXX_L5_2.fq.gz Raw_Reads/HV_017_02_reverse.fq.gz
mv HV017_3_USPD16102684-3_H3FVNBBXX_L5_2.fq.gz Raw_Reads/HV_017_03_reverse.fq.gz
mv HV017_4_USPD16102684-4_H3FVNBBXX_L5_2.fq.gz Raw_Reads/HV_017_04_reverse.fq.gz
mv HV017_5_USPD16102684-5_H3FVNBBXX_L5_2.fq.gz Raw_Reads/HV_017_05_reverse.fq.gz
mv HV017_6_USPD16102684-6_H3FVNBBXX_L5_2.fq.gz Raw_Reads/HV_017_06_reverse.fq.gz
mv HV017_7_USPD16102684-7_H3FVNBBXX_L5_2.fq.gz Raw_Reads/HV_017_07_reverse.fq.gz
mv HV017_8_USPD16102684-8_H3FVNBBXX_L5_2.fq.gz Raw_Reads/HV_017_08_reverse.fq.gz
mv HV017_9_USPD16102684-9_H3FVNBBXX_L5_2.fq.gz Raw_Reads/HV_017_09_reverse.fq.gz
mv HV017_10_USPD16102684-10_H3FVNBBXX_L5_2.fq.gz Raw_Reads/HV_017_10_reverse.fq.gz
mv HV017_11_USPD16102684-11_H3FVNBBXX_L5_2.fq.gz Raw_Reads/HV_017_11_reverse.fq.gz
mv HV017_12_USPD16102684-12_H3FVNBBXX_L5_2.fq.gz Raw_Reads/HV_017_12_reverse.fq.gz
mv HV017_13_USPD16102684-13_H3FVNBBXX_L5_2.fq.gz Raw_Reads/HV_017_13_reverse.fq.gz
mv HV017_14_USPD16102684-14_H3FVNBBXX_L5_2.fq.gz Raw_Reads/HV_017_14_reverse.fq.gz
mv HV017_15_USPD16102684-15_H3FVNBBXX_L5_2.fq.gz Raw_Reads/HV_017_15_reverse.fq.gz
mv HV017_16_USPD16102684-16_H3FVNBBXX_L5_2.fq.gz Raw_Reads/HV_017_16_reverse.fq.gz
mv HV017_18_USPD16102684-18_H3FVNBBXX_L5_2.fq.gz Raw_Reads/HV_017_18_reverse.fq.gz
mv HV017_19_USPD16102684-19_H3FVNBBXX_L5_2.fq.gz Raw_Reads/HV_017_19_reverse.fq.gz
mv HV017_20_USPD16102684-20_H3FVNBBXX_L5_2.fq.gz Raw_Reads/HV_017_20_reverse.fq.gz
mv HV017_21_USPD16102684-21_H3FVNBBXX_L5_2.fq.gz Raw_Reads/HV_017_21_reverse.fq.gz
mv HV017_22_USPD16102684-22_H3FVNBBXX_L5_2.fq.gz Raw_Reads/HV_017_22_reverse.fq.gz
mv HV017_23_USPD16102684-23_H3FVNBBXX_L5_2.fq.gz Raw_Reads/HV_017_23_reverse.fq.gz
mv HV017_25_USPD16102684-25_H3FVNBBXX_L5_2.fq.gz Raw_Reads/HV_017_25_reverse.fq.gz
mv HV017_27_USPD16102684-27_H3FVNBBXX_L5_2.fq.gz Raw_Reads/HV_017_27_reverse.fq.gz
mv HV018_01_2.fq.gz Raw_Reads/HV_018_01_reverse.fq.gz
mv HV018_02_2.fq.gz Raw_Reads/HV_018_02_reverse.fq.gz
mv HV018_03_2.fq.gz Raw_Reads/HV_018_03_reverse.fq.gz
mv HV018_04_2.fq.gz Raw_Reads/HV_018_04_reverse.fq.gz
mv HV018_05_2.fq.gz Raw_Reads/HV_018_05_reverse.fq.gz
mv HV018_06_2.fq.gz Raw_Reads/HV_018_06_reverse.fq.gz
mv HV018_07_2.fq.gz Raw_Reads/HV_018_07_reverse.fq.gz
mv HV018_08_2.fq.gz Raw_Reads/HV_018_08_reverse.fq.gz
mv HV018_09_2.fq.gz Raw_Reads/HV_018_09_reverse.fq.gz
mv HV018_10_2.fq.gz Raw_Reads/HV_018_10_reverse.fq.gz
mv HV018_11_2.fq.gz Raw_Reads/HV_018_11_reverse.fq.gz
mv HV018_12_2.fq.gz Raw_Reads/HV_018_12_reverse.fq.gz
mv HV018_13_2.fq.gz Raw_Reads/HV_018_13_reverse.fq.gz
mv HV018_14_2.fq.gz Raw_Reads/HV_018_14_reverse.fq.gz
mv HV018_15_2.fq.gz Raw_Reads/HV_018_15_reverse.fq.gz
mv HV018_16_2.fq.gz Raw_Reads/HV_018_16_reverse.fq.gz
mv HV018_18_2.fq.gz Raw_Reads/HV_018_18_reverse.fq.gz
mv HV018_19_2.fq.gz Raw_Reads/HV_018_19_reverse.fq.gz
mv HV018_20_2.fq.gz Raw_Reads/HV_018_20_reverse.fq.gz
mv HV018_21_2.fq.gz Raw_Reads/HV_018_21_reverse.fq.gz
mv HV018_22_2.fq.gz Raw_Reads/HV_018_22_reverse.fq.gz
mv HV018_23_2.fq.gz Raw_Reads/HV_018_23_reverse.fq.gz
mv HV018_25_2.fq.gz Raw_Reads/HV_018_25_reverse.fq.gz
mv HV018_27_2.fq.gz Raw_Reads/HV_018_27_reverse.fq.gz
mv HV019_01_2.fq.gz Raw_Reads/HV_019_01_reverse.fq.gz
mv HV019_02_2.fq.gz Raw_Reads/HV_019_02_reverse.fq.gz
mv HV019_03_2.fq.gz Raw_Reads/HV_019_03_reverse.fq.gz
mv HV019_04_2.fq.gz Raw_Reads/HV_019_04_reverse.fq.gz
mv HV019_05_2.fq.gz Raw_Reads/HV_019_05_reverse.fq.gz
mv HV019_06_2.fq.gz Raw_Reads/HV_019_06_reverse.fq.gz
mv HV019_07_2.fq.gz Raw_Reads/HV_019_07_reverse.fq.gz
mv HV019_08_2.fq.gz Raw_Reads/HV_019_08_reverse.fq.gz
mv HV019_09_2.fq.gz Raw_Reads/HV_019_09_reverse.fq.gz
mv HV019_10_2.fq.gz Raw_Reads/HV_019_10_reverse.fq.gz
mv HV019_11_2.fq.gz Raw_Reads/HV_019_11_reverse.fq.gz
mv HV019_12_2.fq.gz Raw_Reads/HV_019_12_reverse.fq.gz
mv HV019_13_2.fq.gz Raw_Reads/HV_019_13_reverse.fq.gz
mv HV019_14_2.fq.gz Raw_Reads/HV_019_14_reverse.fq.gz
mv HV019_15_2.fq.gz Raw_Reads/HV_019_15_reverse.fq.gz
mv HV019_16_2.fq.gz Raw_Reads/HV_019_16_reverse.fq.gz
mv HV019_18_2.fq.gz Raw_Reads/HV_019_18_reverse.fq.gz
mv HV019_19_2.fq.gz Raw_Reads/HV_019_19_reverse.fq.gz
mv HV019_20_2.fq.gz Raw_Reads/HV_019_20_reverse.fq.gz
mv HV019_21_2.fq.gz Raw_Reads/HV_019_21_reverse.fq.gz
mv HV019_22_2.fq.gz Raw_Reads/HV_019_22_reverse.fq.gz
mv HV019_23_2.fq.gz Raw_Reads/HV_019_23_reverse.fq.gz
mv HV019_25_2.fq.gz Raw_Reads/HV_019_25_reverse.fq.gz
mv HV019_27_2.fq.gz Raw_Reads/HV_019_27_reverse.fq.gz
mv HV020_01_2.fq.gz Raw_Reads/HV_020_01_reverse.fq.gz
mv HV020_02_2.fq.gz Raw_Reads/HV_020_02_reverse.fq.gz
mv HV020_03_2.fq.gz Raw_Reads/HV_020_03_reverse.fq.gz
mv HV020_04_2.fq.gz Raw_Reads/HV_020_04_reverse.fq.gz
mv HV020_05_2.fq.gz Raw_Reads/HV_020_05_reverse.fq.gz
mv HV020_06_2.fq.gz Raw_Reads/HV_020_06_reverse.fq.gz
mv HV020_07_2.fq.gz Raw_Reads/HV_020_07_reverse.fq.gz
mv HV020_08_2.fq.gz Raw_Reads/HV_020_08_reverse.fq.gz
mv HV020_09_2.fq.gz Raw_Reads/HV_020_09_reverse.fq.gz
mv HV020_10_2.fq.gz Raw_Reads/HV_020_10_reverse.fq.gz
mv HV020_11_2.fq.gz Raw_Reads/HV_020_11_reverse.fq.gz
mv HV020_12_2.fq.gz Raw_Reads/HV_020_12_reverse.fq.gz
mv HV020_13_2.fq.gz Raw_Reads/HV_020_13_reverse.fq.gz
mv HV020_14_2.fq.gz Raw_Reads/HV_020_14_reverse.fq.gz
mv HV020_15_2.fq.gz Raw_Reads/HV_020_15_reverse.fq.gz
mv HV020_16_2.fq.gz Raw_Reads/HV_020_16_reverse.fq.gz
mv HV020_18_2.fq.gz Raw_Reads/HV_020_18_reverse.fq.gz
mv HV020_19_2.fq.gz Raw_Reads/HV_020_19_reverse.fq.gz
mv HV020_20_2.fq.gz Raw_Reads/HV_020_20_reverse.fq.gz
mv HV020_21_2.fq.gz Raw_Reads/HV_020_21_reverse.fq.gz
mv HV020_22_2.fq.gz Raw_Reads/HV_020_22_reverse.fq.gz
mv HV020_23_2.fq.gz Raw_Reads/HV_020_23_reverse.fq.gz
mv HV020_25_2.fq.gz Raw_Reads/HV_020_25_reverse.fq.gz
mv HV020_27_2.fq.gz Raw_Reads/HV_020_27_reverse.fq.gz
mv HV021_01_2.fq.gz Raw_Reads/HV_021_01_reverse.fq.gz
mv HV021_02_2.fq.gz Raw_Reads/HV_021_02_reverse.fq.gz
mv HV021_03_2.fq.gz Raw_Reads/HV_021_03_reverse.fq.gz
mv HV021_04_2.fq.gz Raw_Reads/HV_021_04_reverse.fq.gz
mv HV021_05_2.fq.gz Raw_Reads/HV_021_05_reverse.fq.gz
mv HV021_06_2.fq.gz Raw_Reads/HV_021_06_reverse.fq.gz
mv HV021_07_2.fq.gz Raw_Reads/HV_021_07_reverse.fq.gz
mv HV021_08_2.fq.gz Raw_Reads/HV_021_08_reverse.fq.gz
mv HV021_09_2.fq.gz Raw_Reads/HV_021_09_reverse.fq.gz
mv HV021_10_2.fq.gz Raw_Reads/HV_021_10_reverse.fq.gz
mv HV021_11_2.fq.gz Raw_Reads/HV_021_11_reverse.fq.gz
mv HV021_12_2.fq.gz Raw_Reads/HV_021_12_reverse.fq.gz
mv HV021_13_2.fq.gz Raw_Reads/HV_021_13_reverse.fq.gz
mv HV021_14_2.fq.gz Raw_Reads/HV_021_14_reverse.fq.gz
mv HV021_15_2.fq.gz Raw_Reads/HV_021_15_reverse.fq.gz
mv HV021_16_2.fq.gz Raw_Reads/HV_021_16_reverse.fq.gz
mv HV021_18_2.fq.gz Raw_Reads/HV_021_18_reverse.fq.gz
mv HV021_19_2.fq.gz Raw_Reads/HV_021_19_reverse.fq.gz
mv HV021_20_2.fq.gz Raw_Reads/HV_021_20_reverse.fq.gz
mv HV021_21_2.fq.gz Raw_Reads/HV_021_21_reverse.fq.gz
mv HV021_22_2.fq.gz Raw_Reads/HV_021_22_reverse.fq.gz
mv HV021_23_2.fq.gz Raw_Reads/HV_021_23_reverse.fq.gz
mv HV021_25_2.fq.gz Raw_Reads/HV_021_25_reverse.fq.gz
mv HV021_27_2.fq.gz Raw_Reads/HV_021_27_reverse.fq.gz
mv HV_022_CKDL190143346-1a-1_H75V2BBXX_L2_2.fq.gz Raw_Reads/HV_022_01_reverse.fq.gz
mv HV_022_CKDL190143346-1a-2_H75V2BBXX_L2_2.fq.gz Raw_Reads/HV_022_02_reverse.fq.gz
mv HV_022_CKDL190143346-1a-3_H75V2BBXX_L2_2.fq.gz Raw_Reads/HV_022_03_reverse.fq.gz
mv HV_022_CKDL190143346-1a-4_H75V2BBXX_L2_2.fq.gz Raw_Reads/HV_022_04_reverse.fq.gz
mv HV_022_CKDL190143346-1a-5_H75V2BBXX_L2_2.fq.gz Raw_Reads/HV_022_05_reverse.fq.gz
mv HV_022_CKDL190143346-1a-6_H75V2BBXX_L2_2.fq.gz Raw_Reads/HV_022_06_reverse.fq.gz
mv HV_022_CKDL190143346-1a-7_H75V2BBXX_L2_2.fq.gz Raw_Reads/HV_022_07_reverse.fq.gz
mv HV_022_CKDL190143346-1a-8_H75V2BBXX_L2_2.fq.gz Raw_Reads/HV_022_08_reverse.fq.gz
mv HV_022_CKDL190143346-1a-9_H75V2BBXX_L2_2.fq.gz Raw_Reads/HV_022_09_reverse.fq.gz
mv HV_022_CKDL190143346-1a-10_H75V2BBXX_L2_2.fq.gz Raw_Reads/HV_022_10_reverse.fq.gz
mv HV_022_CKDL190143346-1a-11_H75V2BBXX_L2_2.fq.gz Raw_Reads/HV_022_11_reverse.fq.gz
mv HV_022_CKDL190143346-1a-12_H75V2BBXX_L2_2.fq.gz Raw_Reads/HV_022_12_reverse.fq.gz
mv HV_022_CKDL190143346-1a-13_H75V2BBXX_L2_2.fq.gz Raw_Reads/HV_022_13_reverse.fq.gz
mv HV_022_CKDL190143346-1a-14_H75V2BBXX_L2_2.fq.gz Raw_Reads/HV_022_14_reverse.fq.gz
mv HV_022_CKDL190143346-1a-15_H75V2BBXX_L2_2.fq.gz Raw_Reads/HV_022_15_reverse.fq.gz
mv HV_022_CKDL190143346-1a-16_H75V2BBXX_L2_2.fq.gz Raw_Reads/HV_022_16_reverse.fq.gz
mv HV_022_CKDL190143346-1a-18_H75V2BBXX_L2_2.fq.gz Raw_Reads/HV_022_18_reverse.fq.gz
mv HV_022_CKDL190143346-1a-19_H75V2BBXX_L2_2.fq.gz Raw_Reads/HV_022_19_reverse.fq.gz
mv HV_022_CKDL190143346-1a-20_H75V2BBXX_L2_2.fq.gz Raw_Reads/HV_022_20_reverse.fq.gz
mv HV_022_CKDL190143346-1a-21_H75V2BBXX_L2_2.fq.gz Raw_Reads/HV_022_21_reverse.fq.gz
mv HV_022_CKDL190143346-1a-22_H75V2BBXX_L2_2.fq.gz Raw_Reads/HV_022_22_reverse.fq.gz
mv HV_022_CKDL190143346-1a-23_H75V2BBXX_L2_2.fq.gz Raw_Reads/HV_022_23_reverse.fq.gz
mv HV_022_CKDL190143346-1a-25_H75V2BBXX_L2_2.fq.gz Raw_Reads/HV_022_25_reverse.fq.gz
mv HV_022_CKDL190143346-1a-27_H75V2BBXX_L2_2.fq.gz Raw_Reads/HV_022_27_reverse.fq.gz
mv HV_023_CKDL190143347-1a-1_H75V2BBXX_L3_2.fq.gz Raw_Reads/HV_023_01_reverse.fq.gz
mv HV_023_CKDL190143347-1a-2_H75V2BBXX_L3_2.fq.gz Raw_Reads/HV_023_02_reverse.fq.gz
mv HV_023_CKDL190143347-1a-3_H75V2BBXX_L3_2.fq.gz Raw_Reads/HV_023_03_reverse.fq.gz
mv HV_023_CKDL190143347-1a-4_H75V2BBXX_L3_2.fq.gz Raw_Reads/HV_023_04_reverse.fq.gz
mv HV_023_CKDL190143347-1a-5_H75V2BBXX_L3_2.fq.gz Raw_Reads/HV_023_05_reverse.fq.gz
mv HV_023_CKDL190143347-1a-6_H75V2BBXX_L3_2.fq.gz Raw_Reads/HV_023_06_reverse.fq.gz
mv HV_023_CKDL190143347-1a-7_H75V2BBXX_L3_2.fq.gz Raw_Reads/HV_023_07_reverse.fq.gz
mv HV_023_CKDL190143347-1a-8_H75V2BBXX_L3_2.fq.gz Raw_Reads/HV_023_08_reverse.fq.gz
mv HV_023_CKDL190143347-1a-9_H75V2BBXX_L3_2.fq.gz Raw_Reads/HV_023_09_reverse.fq.gz
mv HV_023_CKDL190143347-1a-10_H75V2BBXX_L3_2.fq.gz Raw_Reads/HV_023_10_reverse.fq.gz
mv HV_023_CKDL190143347-1a-11_H75V2BBXX_L3_2.fq.gz Raw_Reads/HV_023_11_reverse.fq.gz
mv HV_023_CKDL190143347-1a-12_H75V2BBXX_L3_2.fq.gz Raw_Reads/HV_023_12_reverse.fq.gz
mv HV_023_CKDL190143347-1a-13_H75V2BBXX_L3_2.fq.gz Raw_Reads/HV_023_13_reverse.fq.gz
mv HV_023_CKDL190143347-1a-14_H75V2BBXX_L3_2.fq.gz Raw_Reads/HV_023_14_reverse.fq.gz
mv HV_023_CKDL190143347-1a-15_H75V2BBXX_L3_2.fq.gz Raw_Reads/HV_023_15_reverse.fq.gz
mv HV_023_CKDL190143347-1a-16_H75V2BBXX_L3_2.fq.gz Raw_Reads/HV_023_16_reverse.fq.gz
mv HV_023_CKDL190143347-1a-18_H75V2BBXX_L3_2.fq.gz Raw_Reads/HV_023_18_reverse.fq.gz
mv HV_023_CKDL190143347-1a-19_H75V2BBXX_L3_2.fq.gz Raw_Reads/HV_023_19_reverse.fq.gz
mv HV_023_CKDL190143347-1a-20_H75V2BBXX_L3_2.fq.gz Raw_Reads/HV_023_20_reverse.fq.gz
mv HV_023_CKDL190143347-1a-21_H75V2BBXX_L3_2.fq.gz Raw_Reads/HV_023_21_reverse.fq.gz
mv HV_023_CKDL190143347-1a-22_H75V2BBXX_L3_2.fq.gz Raw_Reads/HV_023_22_reverse.fq.gz
mv HV_023_CKDL190143347-1a-23_H75V2BBXX_L3_2.fq.gz Raw_Reads/HV_023_23_reverse.fq.gz
mv HV_023_CKDL190143347-1a-25_H75V2BBXX_L3_2.fq.gz Raw_Reads/HV_023_25_reverse.fq.gz
mv HV_023_CKDL190143347-1a-27_H75V2BBXX_L3_2.fq.gz Raw_Reads/HV_023_27_reverse.fq.gz
mv HV_024_CKDL190143348-1a-1_H75V2BBXX_L4_2.fq.gz Raw_Reads/HV_024_01_reverse.fq.gz
mv HV_024_CKDL190143348-1a-2_H75V2BBXX_L4_2.fq.gz Raw_Reads/HV_024_02_reverse.fq.gz
mv HV_024_CKDL190143348-1a-3_H75V2BBXX_L4_2.fq.gz Raw_Reads/HV_024_03_reverse.fq.gz
mv HV_024_CKDL190143348-1a-4_H75V2BBXX_L4_2.fq.gz Raw_Reads/HV_024_04_reverse.fq.gz
mv HV_024_CKDL190143348-1a-5_H75V2BBXX_L4_2.fq.gz Raw_Reads/HV_024_05_reverse.fq.gz
mv HV_024_CKDL190143348-1a-6_H75V2BBXX_L4_2.fq.gz Raw_Reads/HV_024_06_reverse.fq.gz
mv HV_024_CKDL190143348-1a-7_H75V2BBXX_L4_2.fq.gz Raw_Reads/HV_024_07_reverse.fq.gz
mv HV_024_CKDL190143348-1a-8_H75V2BBXX_L4_2.fq.gz Raw_Reads/HV_024_08_reverse.fq.gz
mv HV_024_CKDL190143348-1a-9_H75V2BBXX_L4_2.fq.gz Raw_Reads/HV_024_09_reverse.fq.gz
mv HV_024_CKDL190143348-1a-10_H75V2BBXX_L4_2.fq.gz Raw_Reads/HV_024_10_reverse.fq.gz
mv HV_024_CKDL190143348-1a-11_H75V2BBXX_L4_2.fq.gz Raw_Reads/HV_024_11_reverse.fq.gz
mv HV_024_CKDL190143348-1a-12_H75V2BBXX_L4_2.fq.gz Raw_Reads/HV_024_12_reverse.fq.gz
mv HV_024_CKDL190143348-1a-13_H75V2BBXX_L4_2.fq.gz Raw_Reads/HV_024_13_reverse.fq.gz
mv HV_024_CKDL190143348-1a-14_H75V2BBXX_L4_2.fq.gz Raw_Reads/HV_024_14_reverse.fq.gz
mv HV_024_CKDL190143348-1a-15_H75V2BBXX_L4_2.fq.gz Raw_Reads/HV_024_15_reverse.fq.gz
mv HV_024_CKDL190143348-1a-16_H75V2BBXX_L4_2.fq.gz Raw_Reads/HV_024_16_reverse.fq.gz
mv HV_024_CKDL190143348-1a-18_H75V2BBXX_L4_2.fq.gz Raw_Reads/HV_024_18_reverse.fq.gz
mv HV_024_CKDL190143348-1a-19_H75V2BBXX_L4_2.fq.gz Raw_Reads/HV_024_19_reverse.fq.gz
mv HV_024_CKDL190143348-1a-20_H75V2BBXX_L4_2.fq.gz Raw_Reads/HV_024_20_reverse.fq.gz
mv HV_024_CKDL190143348-1a-21_H75V2BBXX_L4_2.fq.gz Raw_Reads/HV_024_21_reverse.fq.gz
mv HV_024_CKDL190143348-1a-22_H75V2BBXX_L4_2.fq.gz Raw_Reads/HV_024_22_reverse.fq.gz
mv HV_024_CKDL190143348-1a-23_H75V2BBXX_L4_2.fq.gz Raw_Reads/HV_024_23_reverse.fq.gz
mv HV_024_CKDL190143348-1a-25_H75V2BBXX_L4_2.fq.gz Raw_Reads/HV_024_25_reverse.fq.gz
mv HV_024_CKDL190143348-1a-27_H75V2BBXX_L4_2.fq.gz Raw_Reads/HV_024_27_reverse.fq.gz
mv HV_025_CKDL190143349-1a-1_H75V2BBXX_L5_2.fq.gz Raw_Reads/HV_025_01_reverse.fq.gz
mv HV_025_CKDL190143349-1a-2_H75V2BBXX_L5_2.fq.gz Raw_Reads/HV_025_02_reverse.fq.gz
mv HV_025_CKDL190143349-1a-3_H75V2BBXX_L5_2.fq.gz Raw_Reads/HV_025_03_reverse.fq.gz
mv HV_025_CKDL190143349-1a-4_H75V2BBXX_L5_2.fq.gz Raw_Reads/HV_025_04_reverse.fq.gz
mv HV_025_CKDL190143349-1a-5_H75V2BBXX_L5_2.fq.gz Raw_Reads/HV_025_05_reverse.fq.gz
mv HV_025_CKDL190143349-1a-6_H75V2BBXX_L5_2.fq.gz Raw_Reads/HV_025_06_reverse.fq.gz
mv HV_025_CKDL190143349-1a-7_H75V2BBXX_L5_2.fq.gz Raw_Reads/HV_025_07_reverse.fq.gz
mv HV_025_CKDL190143349-1a-8_H75V2BBXX_L5_2.fq.gz Raw_Reads/HV_025_08_reverse.fq.gz
mv HV_025_CKDL190143349-1a-9_H75V2BBXX_L5_2.fq.gz Raw_Reads/HV_025_09_reverse.fq.gz
mv HV_025_CKDL190143349-1a-10_H75V2BBXX_L5_2.fq.gz Raw_Reads/HV_025_10_reverse.fq.gz
mv HV_025_CKDL190143349-1a-11_H75V2BBXX_L5_2.fq.gz Raw_Reads/HV_025_11_reverse.fq.gz
mv HV_025_CKDL190143349-1a-12_H75V2BBXX_L5_2.fq.gz Raw_Reads/HV_025_12_reverse.fq.gz
mv HV_025_CKDL190143349-1a-13_H75V2BBXX_L5_2.fq.gz Raw_Reads/HV_025_13_reverse.fq.gz
mv HV_025_CKDL190143349-1a-14_H75V2BBXX_L5_2.fq.gz Raw_Reads/HV_025_14_reverse.fq.gz
mv HV_025_CKDL190143349-1a-15_H75V2BBXX_L5_2.fq.gz Raw_Reads/HV_025_15_reverse.fq.gz
mv HV_025_CKDL190143349-1a-16_H75V2BBXX_L5_2.fq.gz Raw_Reads/HV_025_16_reverse.fq.gz
mv HV_025_CKDL190143349-1a-18_H75V2BBXX_L5_2.fq.gz Raw_Reads/HV_025_18_reverse.fq.gz
mv HV_025_CKDL190143349-1a-19_H75V2BBXX_L5_2.fq.gz Raw_Reads/HV_025_19_reverse.fq.gz
mv HV_025_CKDL190143349-1a-20_H75V2BBXX_L5_2.fq.gz Raw_Reads/HV_025_20_reverse.fq.gz
mv HV_025_CKDL190143349-1a-21_H75V2BBXX_L5_2.fq.gz Raw_Reads/HV_025_21_reverse.fq.gz
mv HV_025_CKDL190143349-1a-22_H75V2BBXX_L5_2.fq.gz Raw_Reads/HV_025_22_reverse.fq.gz
mv HV_025_CKDL190143349-1a-23_H75V2BBXX_L5_2.fq.gz Raw_Reads/HV_025_23_reverse.fq.gz
mv HV_025_CKDL190143349-1a-25_H75V2BBXX_L5_2.fq.gz Raw_Reads/HV_025_25_reverse.fq.gz
mv HV_025_CKDL190143349-1a-27_H75V2BBXX_L5_2.fq.gz Raw_Reads/HV_025_27_reverse.fq.gz
mv HV_026_CKDL190144758-1a-1_H75W3BBXX_L4_2.fq.gz Raw_Reads/HV_026_01_reverse.fq.gz
mv HV_026_CKDL190144758-1a-2_H75W3BBXX_L4_2.fq.gz Raw_Reads/HV_026_02_reverse.fq.gz
mv HV_026_CKDL190144758-1a-3_H75W3BBXX_L4_2.fq.gz Raw_Reads/HV_026_03_reverse.fq.gz
mv HV_026_CKDL190144758-1a-4_H75W3BBXX_L4_2.fq.gz Raw_Reads/HV_026_04_reverse.fq.gz
mv HV_026_CKDL190144758-1a-5_H75W3BBXX_L4_2.fq.gz Raw_Reads/HV_026_05_reverse.fq.gz
mv HV_026_CKDL190144758-1a-6_H75W3BBXX_L4_2.fq.gz Raw_Reads/HV_026_06_reverse.fq.gz
mv HV_026_CKDL190144758-1a-7_H75W3BBXX_L4_2.fq.gz Raw_Reads/HV_026_07_reverse.fq.gz
mv HV_026_CKDL190144758-1a-8_H75W3BBXX_L4_2.fq.gz Raw_Reads/HV_026_08_reverse.fq.gz
mv HV_026_CKDL190144758-1a-9_H75W3BBXX_L4_2.fq.gz Raw_Reads/HV_026_09_reverse.fq.gz
mv HV_026_CKDL190144758-1a-10_H75W3BBXX_L4_2.fq.gz Raw_Reads/HV_026_10_reverse.fq.gz
mv HV_026_CKDL190144758-1a-11_H75W3BBXX_L4_2.fq.gz Raw_Reads/HV_026_11_reverse.fq.gz
mv HV_026_CKDL190144758-1a-12_H75W3BBXX_L4_2.fq.gz Raw_Reads/HV_026_12_reverse.fq.gz
mv HV_026_CKDL190144758-1a-13_H75W3BBXX_L4_2.fq.gz Raw_Reads/HV_026_13_reverse.fq.gz
mv HV_026_CKDL190144758-1a-14_H75W3BBXX_L4_2.fq.gz Raw_Reads/HV_026_14_reverse.fq.gz
mv HV_026_CKDL190144758-1a-15_H75W3BBXX_L4_2.fq.gz Raw_Reads/HV_026_15_reverse.fq.gz
mv HV_026_CKDL190144758-1a-16_H75W3BBXX_L4_2.fq.gz Raw_Reads/HV_026_16_reverse.fq.gz
mv HV_026_CKDL190144758-1a-18_H75W3BBXX_L4_2.fq.gz Raw_Reads/HV_026_18_reverse.fq.gz
mv HV_026_CKDL190144758-1a-19_H75W3BBXX_L4_2.fq.gz Raw_Reads/HV_026_19_reverse.fq.gz
mv HV_026_CKDL190144758-1a-20_H75W3BBXX_L4_2.fq.gz Raw_Reads/HV_026_20_reverse.fq.gz
mv HV_026_CKDL190144758-1a-21_H75W3BBXX_L4_2.fq.gz Raw_Reads/HV_026_21_reverse.fq.gz
mv HV_026_CKDL190144758-1a-22_H75W3BBXX_L4_2.fq.gz Raw_Reads/HV_026_22_reverse.fq.gz
mv HV_026_CKDL190144758-1a-23_H75W3BBXX_L4_2.fq.gz Raw_Reads/HV_026_23_reverse.fq.gz
mv HV_026_CKDL190144758-1a-25_H75W3BBXX_L4_2.fq.gz Raw_Reads/HV_026_25_reverse.fq.gz
mv HV_026_CKDL190144758-1a-27_H75W3BBXX_L4_2.fq.gz Raw_Reads/HV_026_27_reverse.fq.gz
mv HV_027_CKDL190144759-1a-1_H75TMBBXX_L8_2.fq.gz Raw_Reads/HV_027_01_reverse.fq.gz
mv HV_027_CKDL190144759-1a-2_H75TMBBXX_L8_2.fq.gz Raw_Reads/HV_027_02_reverse.fq.gz
mv HV_027_CKDL190144759-1a-3_H75TMBBXX_L8_2.fq.gz Raw_Reads/HV_027_03_reverse.fq.gz
mv HV_027_CKDL190144759-1a-4_H75TMBBXX_L8_2.fq.gz Raw_Reads/HV_027_04_reverse.fq.gz
mv HV_027_CKDL190144759-1a-5_H75TMBBXX_L8_2.fq.gz Raw_Reads/HV_027_05_reverse.fq.gz
mv HV_027_CKDL190144759-1a-6_H75TMBBXX_L8_2.fq.gz Raw_Reads/HV_027_06_reverse.fq.gz
mv HV_027_CKDL190144759-1a-7_H75TMBBXX_L8_2.fq.gz Raw_Reads/HV_027_07_reverse.fq.gz
mv HV_027_CKDL190144759-1a-8_H75TMBBXX_L8_2.fq.gz Raw_Reads/HV_027_08_reverse.fq.gz
mv HV_027_CKDL190144759-1a-9_H75TMBBXX_L8_2.fq.gz Raw_Reads/HV_027_09_reverse.fq.gz
mv HV_027_CKDL190144759-1a-10_H75TMBBXX_L8_2.fq.gz Raw_Reads/HV_027_10_reverse.fq.gz
mv HV_027_CKDL190144759-1a-11_H75TMBBXX_L8_2.fq.gz Raw_Reads/HV_027_11_reverse.fq.gz
mv HV_027_CKDL190144759-1a-12_H75TMBBXX_L8_2.fq.gz Raw_Reads/HV_027_12_reverse.fq.gz
mv HV_027_CKDL190144759-1a-13_H75TMBBXX_L8_2.fq.gz Raw_Reads/HV_027_13_reverse.fq.gz
mv HV_027_CKDL190144759-1a-14_H75TMBBXX_L8_2.fq.gz Raw_Reads/HV_027_14_reverse.fq.gz
mv HV_027_CKDL190144759-1a-15_H75TMBBXX_L8_2.fq.gz Raw_Reads/HV_027_15_reverse.fq.gz
mv HV_027_CKDL190144759-1a-16_H75TMBBXX_L8_2.fq.gz Raw_Reads/HV_027_16_reverse.fq.gz
mv HV_027_CKDL190144759-1a-18_H75TMBBXX_L8_2.fq.gz Raw_Reads/HV_027_18_reverse.fq.gz
mv HV_027_CKDL190144759-1a-19_H75TMBBXX_L8_2.fq.gz Raw_Reads/HV_027_19_reverse.fq.gz
mv HV_027_CKDL190144759-1a-20_H75TMBBXX_L8_2.fq.gz Raw_Reads/HV_027_20_reverse.fq.gz
mv HV_027_CKDL190144759-1a-21_H75TMBBXX_L8_2.fq.gz Raw_Reads/HV_027_21_reverse.fq.gz
mv HV_027_CKDL190144759-1a-22_H75TMBBXX_L8_2.fq.gz Raw_Reads/HV_027_22_reverse.fq.gz
mv HV_027_CKDL190144759-1a-23_H75TMBBXX_L8_2.fq.gz Raw_Reads/HV_027_23_reverse.fq.gz
mv HV_027_CKDL190144759-1a-25_H75TMBBXX_L8_2.fq.gz Raw_Reads/HV_027_25_reverse.fq.gz
mv HV_027_CKDL190144759-1a-27_H75TMBBXX_L8_2.fq.gz Raw_Reads/HV_027_27_reverse.fq.gz
mv HV028_1_CKDL200148804-1a-1_H7HJWBBXX_L3_2.fq.gz Raw_Reads/HV_028_01_reverse.fq.gz
mv HV028_2_CKDL200148804-1a-2_H7HJWBBXX_L3_2.fq.gz Raw_Reads/HV_028_02_reverse.fq.gz
mv HV028_3_CKDL200148804-1a-3_H7HJWBBXX_L3_2.fq.gz Raw_Reads/HV_028_03_reverse.fq.gz
mv HV028_4_CKDL200148804-1a-4_H7HJWBBXX_L3_2.fq.gz Raw_Reads/HV_028_04_reverse.fq.gz
mv HV028_5_CKDL200148804-1a-5_H7HJWBBXX_L3_2.fq.gz Raw_Reads/HV_028_05_reverse.fq.gz
mv HV028_6_CKDL200148804-1a-6_H7HJWBBXX_L3_2.fq.gz Raw_Reads/HV_028_06_reverse.fq.gz
mv HV028_7_CKDL200148804-1a-7_H7HJWBBXX_L3_2.fq.gz Raw_Reads/HV_028_07_reverse.fq.gz
mv HV028_8_CKDL200148804-1a-8_H7HJWBBXX_L3_2.fq.gz Raw_Reads/HV_028_08_reverse.fq.gz
mv HV028_9_CKDL200148804-1a-9_H7HJWBBXX_L3_2.fq.gz Raw_Reads/HV_028_09_reverse.fq.gz
mv HV028_10_CKDL200148804-1a-10_H7HJWBBXX_L3_2.fq.gz Raw_Reads/HV_028_10_reverse.fq.gz
mv HV028_11_CKDL200148804-1a-11_H7HJWBBXX_L3_2.fq.gz Raw_Reads/HV_028_11_reverse.fq.gz
mv HV028_12_CKDL200148804-1a-12_H7HJWBBXX_L3_2.fq.gz Raw_Reads/HV_028_12_reverse.fq.gz
mv HV028_13_CKDL200148804-1a-13_H7HJWBBXX_L3_2.fq.gz Raw_Reads/HV_028_13_reverse.fq.gz
mv HV028_14_CKDL200148804-1a-14_H7HJWBBXX_L3_2.fq.gz Raw_Reads/HV_028_14_reverse.fq.gz
mv HV028_15_CKDL200148804-1a-15_H7HJWBBXX_L3_2.fq.gz Raw_Reads/HV_028_15_reverse.fq.gz
mv HV028_16_CKDL200148804-1a-16_H7HJWBBXX_L3_2.fq.gz Raw_Reads/HV_028_16_reverse.fq.gz
mv HV028_18_CKDL200148804-1a-18_H7HJWBBXX_L3_2.fq.gz Raw_Reads/HV_028_18_reverse.fq.gz
mv HV028_19_CKDL200148804-1a-19_H7HJWBBXX_L3_2.fq.gz Raw_Reads/HV_028_19_reverse.fq.gz
mv HV028_20_CKDL200148804-1a-20_H7HJWBBXX_L3_2.fq.gz Raw_Reads/HV_028_20_reverse.fq.gz
mv HV028_21_CKDL200148804-1a-21_H7HJWBBXX_L3_2.fq.gz Raw_Reads/HV_028_21_reverse.fq.gz
mv HV028_22_CKDL200148804-1a-22_H7HJWBBXX_L3_2.fq.gz Raw_Reads/HV_028_22_reverse.fq.gz
mv HV028_23_CKDL200148804-1a-23_H7HJWBBXX_L3_2.fq.gz Raw_Reads/HV_028_23_reverse.fq.gz
mv HV028_25_CKDL200148804-1a-25_H7HJWBBXX_L3_2.fq.gz Raw_Reads/HV_028_25_reverse.fq.gz
mv HV028_27_CKDL200148804-1a-27_H7HJWBBXX_L3_2.fq.gz Raw_Reads/HV_028_27_reverse.fq.gz
mv HV029_1_CKDL200148805-1a-1_H7HJWBBXX_L4_2.fq.gz Raw_Reads/HV_029_01_reverse.fq.gz
mv HV029_2_CKDL200148805-1a-2_H7HJWBBXX_L4_2.fq.gz Raw_Reads/HV_029_02_reverse.fq.gz
mv HV029_3_CKDL200148805-1a-3_H7HJWBBXX_L4_2.fq.gz Raw_Reads/HV_029_03_reverse.fq.gz
mv HV029_4_CKDL200148805-1a-4_H7HJWBBXX_L4_2.fq.gz Raw_Reads/HV_029_04_reverse.fq.gz
mv HV029_5_CKDL200148805-1a-5_H7HJWBBXX_L4_2.fq.gz Raw_Reads/HV_029_05_reverse.fq.gz
mv HV029_6_CKDL200148805-1a-6_H7HJWBBXX_L4_2.fq.gz Raw_Reads/HV_029_06_reverse.fq.gz
mv HV029_7_CKDL200148805-1a-7_H7HJWBBXX_L4_2.fq.gz Raw_Reads/HV_029_07_reverse.fq.gz
mv HV029_8_CKDL200148805-1a-8_H7HJWBBXX_L4_2.fq.gz Raw_Reads/HV_029_08_reverse.fq.gz
mv HV029_9_CKDL200148805-1a-9_H7HJWBBXX_L4_2.fq.gz Raw_Reads/HV_029_09_reverse.fq.gz
mv HV029_10_CKDL200148805-1a-10_H7HJWBBXX_L4_2.fq.gz Raw_Reads/HV_029_10_reverse.fq.gz
mv HV029_11_CKDL200148805-1a-11_H7HJWBBXX_L4_2.fq.gz Raw_Reads/HV_029_11_reverse.fq.gz
mv HV029_12_CKDL200148805-1a-12_H7HJWBBXX_L4_2.fq.gz Raw_Reads/HV_029_12_reverse.fq.gz
mv HV029_13_CKDL200148805-1a-13_H7HJWBBXX_L4_2.fq.gz Raw_Reads/HV_029_13_reverse.fq.gz
mv HV029_14_CKDL200148805-1a-14_H7HJWBBXX_L4_2.fq.gz Raw_Reads/HV_029_14_reverse.fq.gz
mv HV029_15_CKDL200148805-1a-15_H7HJWBBXX_L4_2.fq.gz Raw_Reads/HV_029_15_reverse.fq.gz
mv HV029_16_CKDL200148805-1a-16_H7HJWBBXX_L4_2.fq.gz Raw_Reads/HV_029_16_reverse.fq.gz
mv HV029_18_CKDL200148805-1a-18_H7HJWBBXX_L4_2.fq.gz Raw_Reads/HV_029_18_reverse.fq.gz
mv HV029_19_CKDL200148805-1a-19_H7HJWBBXX_L4_2.fq.gz Raw_Reads/HV_029_19_reverse.fq.gz
mv HV029_20_CKDL200148805-1a-20_H7HJWBBXX_L4_2.fq.gz Raw_Reads/HV_029_20_reverse.fq.gz
mv HV029_21_CKDL200148805-1a-21_H7HJWBBXX_L4_2.fq.gz Raw_Reads/HV_029_21_reverse.fq.gz
mv HV029_22_CKDL200148805-1a-22_H7HJWBBXX_L4_2.fq.gz Raw_Reads/HV_029_22_reverse.fq.gz
mv HV029_23_CKDL200148805-1a-23_H7HJWBBXX_L4_2.fq.gz Raw_Reads/HV_029_23_reverse.fq.gz
mv HV029_25_CKDL200148805-1a-25_H7HJWBBXX_L4_2.fq.gz Raw_Reads/HV_029_25_reverse.fq.gz
mv HV029_27_CKDL200148805-1a-27_H7HJWBBXX_L4_2.fq.gz Raw_Reads/HV_029_27_reverse.fq.gz
mv HV030_1_CKDL200148806-1a-1_H7HJWBBXX_L5_2.fq.gz Raw_Reads/HV_030_01_reverse.fq.gz
mv HV030_2_CKDL200148806-1a-2_H7HJWBBXX_L5_2.fq.gz Raw_Reads/HV_030_02_reverse.fq.gz
mv HV030_3_CKDL200148806-1a-3_H7HJWBBXX_L5_2.fq.gz Raw_Reads/HV_030_03_reverse.fq.gz
mv HV030_4_CKDL200148806-1a-4_H7HJWBBXX_L5_2.fq.gz Raw_Reads/HV_030_04_reverse.fq.gz
mv HV030_5_CKDL200148806-1a-5_H7HJWBBXX_L5_2.fq.gz Raw_Reads/HV_030_05_reverse.fq.gz
mv HV030_6_CKDL200148806-1a-6_H7HJWBBXX_L5_2.fq.gz Raw_Reads/HV_030_06_reverse.fq.gz
mv HV030_7_CKDL200148806-1a-7_H7HJWBBXX_L5_2.fq.gz Raw_Reads/HV_030_07_reverse.fq.gz
mv HV030_8_CKDL200148806-1a-8_H7HJWBBXX_L5_2.fq.gz Raw_Reads/HV_030_08_reverse.fq.gz
mv HV030_9_CKDL200148806-1a-9_H7HJWBBXX_L5_2.fq.gz Raw_Reads/HV_030_09_reverse.fq.gz
mv HV030_10_CKDL200148806-1a-10_H7HJWBBXX_L5_2.fq.gz Raw_Reads/HV_030_10_reverse.fq.gz
mv HV030_11_CKDL200148806-1a-11_H7HJWBBXX_L5_2.fq.gz Raw_Reads/HV_030_11_reverse.fq.gz
mv HV030_12_CKDL200148806-1a-12_H7HJWBBXX_L5_2.fq.gz Raw_Reads/HV_030_12_reverse.fq.gz
mv HV030_13_CKDL200148806-1a-13_H7HJWBBXX_L5_2.fq.gz Raw_Reads/HV_030_13_reverse.fq.gz
mv HV030_14_CKDL200148806-1a-14_H7HJWBBXX_L5_2.fq.gz Raw_Reads/HV_030_14_reverse.fq.gz
mv HV030_15_CKDL200148806-1a-15_H7HJWBBXX_L5_2.fq.gz Raw_Reads/HV_030_15_reverse.fq.gz
mv HV030_16_CKDL200148806-1a-16_H7HJWBBXX_L5_2.fq.gz Raw_Reads/HV_030_16_reverse.fq.gz
mv HV030_18_CKDL200148806-1a-18_H7HJWBBXX_L5_2.fq.gz Raw_Reads/HV_030_18_reverse.fq.gz
mv HV030_19_CKDL200148806-1a-19_H7HJWBBXX_L5_2.fq.gz Raw_Reads/HV_030_19_reverse.fq.gz
mv HV030_20_CKDL200148806-1a-20_H7HJWBBXX_L5_2.fq.gz Raw_Reads/HV_030_20_reverse.fq.gz
mv HV030_21_CKDL200148806-1a-21_H7HJWBBXX_L5_2.fq.gz Raw_Reads/HV_030_21_reverse.fq.gz
mv HV030_22_CKDL200148806-1a-22_H7HJWBBXX_L5_2.fq.gz Raw_Reads/HV_030_22_reverse.fq.gz
mv HV030_23_CKDL200148806-1a-23_H7HJWBBXX_L5_2.fq.gz Raw_Reads/HV_030_23_reverse.fq.gz
mv HV030_25_CKDL200148806-1a-25_H7HJWBBXX_L5_2.fq.gz Raw_Reads/HV_030_25_reverse.fq.gz
mv HV030_27_CKDL200148806-1a-27_H7HJWBBXX_L5_2.fq.gz Raw_Reads/HV_030_27_reverse.fq.gz

############################################
## -------------------------------------- ##
## ------------ Control Reads ----------- ##
## -------------------------------------- ##
############################################

mv Raw_Reads/*23_forward.fq.gz Raw_Control_Reads/
mv Raw_Reads/*23_reverse.fq.gz Raw_Control_Reads/
mv Raw_Reads/*25_forward.fq.gz Raw_Control_Reads/
mv Raw_Reads/*25_reverse.fq.gz Raw_Control_Reads/

