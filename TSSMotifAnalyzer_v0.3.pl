#!/usr/bin/perl -w

#====================================================================================================================================================#
#<use>
$|++; #---turn on the auto flush for the progress bar
use strict;
use File::Path;
use Time::HiRes qw( time );
use Math::Random;
use Math::Round;
use Storable;
use Getopt::Long;
use File::Basename;
use File::Spec::Functions qw(rel2abs);
use List::Util qw (sum shuffle min max);
use threads;
use threads::shared;
use Statistics::Descriptive;
use URI::Escape;
use Cwd 'abs_path';
#<\use>
#====================================================================================================================================================#

#====================================================================================================================================================#
#<doc>
#	Description
#		This is a Perl script to analyze a predefined set of motifs around the TSS sites.
#
#	Input
#
#		--fullReadDefineNATPileupIndxPath=		path [compulsory]; path of the full read pileup index.hsh.pls which will be used to define genes with or without NAT ;
#		--validTSSCntgPlsPath=					path [compulsory]; a hash sotroable contains gene based TSS information, generated in TEXPolIIDataMerger;;
#		--fastaPath=							file path [compulsory]; the path fasta file contains the genome sequence, for generating blank perl storables;
#		--gffPath=								path[compulsory]; path of the reference GFF for gene annotation;
#		--outDir=								output directory
#	
#	Usage
#		/Volumes/A_MPro2TB/softwareForNGS/myPerlScripts/transcriptModel/TSSMotifAnalyzer/v0.1/TSSMotifAnalyzer_v0.1.pl --gffPath=/Volumes/A_MPro2TB/softwareForNGS/myPerlScripts/pipeLines/marathon/resources/genome/EHI_v2_allFearues.forPileupCounter.gff --fullReadDefineNATPileupIndxPath=/Volumes/A_MPro2TB/NGS/analyses/EHI_NAT/artefactCorrectedPileup/EHI_Iluminia_standard/modelFit/W200.S100.I1.R0.10000000/modelFitData.CL0.999/correctedCovPls.upr/index.hsh.pls --fastaPath=/Volumes/A_MPro2TB/softwareForNGS/myPerlScripts/pipeLines/marathon/resources/genome/EHI_v13.fa --mRNARefPtHshPath=/Volumes/A_MPro2TB/softwareForNGS/myPerlScripts/transcriptModel/TSSFinder/v0.1/cleanedScriptTest/result/MC3.PD75.BF_peak.CR99.9.CC90/storable/mRNARefPtHsh.pls  --gffPath=/Volumes/A_MPro2TB/softwareForNGS/myPerlScripts/pipeLines/marathon/resources/genome/EHI_v2_allFearues.forPileupCounter.gff --outDir=/Volumes/A_MPro2TB/softwareForNGS/myPerlScripts/transcriptModel/TSSMotifAnalyzer/v0.1/betaTest/
#
#	v0.2
#		[13/10/2013 15:34] added validTSSCntgPlsPath option, for analyzing the motifs of all valid TSS;
#		[13/10/2013 15:37] mRNARefPtHshPath option was removed, will used geneBasedTSSInfoHshPlsPath to generate the same information;
#		[13/10/2013 16:53] update the readGff subroutine
#
#	v0.3
#		[13/10/2013 23:17] removed geneBasedTSSInfoHshPlsPath option, directlt get the mRNA reference point from validTSSCntgPlsPath
#<\doc>
#====================================================================================================================================================#

#====================================================================================================================================================#
#<lastCmdCalled>
#
#	[2013-10-13 23:16]	/Volumes/A_MPro2TB/softwareForNGS/myPerlScripts/transcriptModel/TSSMotifAnalyzer/v0.3/TSSMotifAnalyzer_v0.3.pl --fullReadDefineNATPileupIndxPath=/Volumes/A_MPro2TB/NGS/analyses/EHI_NAT/mergedPileup/EHI_Standard/storable/merged/index.hsh.pls --fastaPath=/Volumes/A_MPro2TB/softwareForNGS/resources/genome/AmoebaDB/inUse/AmoebaDB-3.0/EHI/fasta/genome.sorted.fa --validTSSCntgPlsPath=/Volumes/A_MPro2TB/softwareForNGS/myPerlScripts/transcriptModel/TEXPolIIDataMerger/v0.1/TEXPolIIDataMerger/HTr99.9_HPr99.9_HTc95_HPc90_LTr95_lPr90_LTc75_LPc75.UD100.ID.100/storable/TSSCov/index.hsh.pls --gffPath=/Volumes/A_MPro2TB/softwareForNGS/resources/genome/AmoebaDB/inUse/AmoebaDB-3.0/EHI/gff/forPileupCounter.gff --outDir=/Volumes/A_MPro2TB/softwareForNGS/myPerlScripts/transcriptModel/TSSMotifAnalyzer/v0.3/mergedPolIITEXTSS/
#
#	/Volumes/A_MPro2TB/softwareForNGS/myPerlScripts/transcriptModel/TSSMotifAnalyzer/v0.3/TSSMotifAnalyzer_v0.3.pl
#	--fullReadDefineNATPileupIndxPath=/Volumes/A_MPro2TB/NGS/analyses/EHI_NAT/mergedPileup/EHI_Standard/storable/merged/index.hsh.pls
#	--fastaPath=/Volumes/A_MPro2TB/softwareForNGS/resources/genome/AmoebaDB/inUse/AmoebaDB-3.0/EHI/fasta/genome.sorted.fa
#	--validTSSCntgPlsPath=/Volumes/A_MPro2TB/softwareForNGS/myPerlScripts/transcriptModel/TEXPolIIDataMerger/v0.1/TEXPolIIDataMerger/HTr99.999_HPr95_HTc95_HPc90_LTr90_lPr90_LTc50_LPc50.UD100.ID.100/storable/TSSCov/index.hsh.pls
#	--gffPath=/Volumes/A_MPro2TB/softwareForNGS/resources/genome/AmoebaDB/inUse/AmoebaDB-3.0/EHI/gff/forPileupCounter.gff
#	--outDir=/Volumes/A_MPro2TB/softwareForNGS/myPerlScripts/transcriptModel/TSSMotifAnalyzer/v0.3/HTr99.999_HPr95_HTc95_HPc90_LTr90_lPr90_LTc50_LPc50.UD100.ID.100/
#
#<\lastCmdCalled>
#====================================================================================================================================================#

#====================================================================================================================================================#
#<global>
my $scriptDirPath = dirname(rel2abs($0));
my $globalTmpLogPath = "$scriptDirPath/tmp.log.txt";
open TMPLOG, ">", $globalTmpLogPath;
#<\global>
#====================================================================================================================================================#

#====================================================================================================================================================#
{	#Main sections lexical scope starts
#====================================================================================================================================================#

#====================================================================================================================================================#
#	section 0_startingTasks
#	primaryDependOnSub: printCMDLogOrFinishMessage|2737, readParameters|3102
#	secondaryDependOnSub: currentTime|694
#
#<section ID="startingTasks" num="0">
#-----print start log
&printCMDLogOrFinishMessage("CMDLog");#->2737

#-----read parameters
my ($fullReadDefineNATPileupIndxPath, $validTSSCntgPlsPath, $fastaPath, $gffPath, $outDir) = &readParameters();#->3102
#<\section>
#====================================================================================================================================================#

#====================================================================================================================================================#
#	section 1_defineHardCodedParam
#	primaryDependOnSub: defineHardCodedParameters|841
#	secondaryDependOnSub: >none
#
#<section ID="defineHardCodedParam" num="1">
my ($forceDefineEnd3NAT, $siteSeqRng, $cmltvePprtnLimit, $maxPolymerSize, $extraMargin, $abvBkgrdFactor, $end3NATDefineParamHsh_ref, $predictTSSMotifInfoHsh_ref, $minStreakHsh_ref, $resultDirTag, $TSSScoreResultTag, $logTransformPVal, $skipEndNt, $TSScoreRegSizeHsh_ref, $forceCalculateTSSScoreIndivmRNA) = &defineHardCodedParameters();#->841
#---used in defineEnd3NATInmRNA

#<\section>
#====================================================================================================================================================#

#====================================================================================================================================================#
#	section 2_defineHardCodedPath
#	primaryDependOnSub: >none
#	secondaryDependOnSub: >none
#
#<section ID="defineHardCodedPath" num="2">
my $coreMotifPath = "/Volumes/A_MPro2TB/softwareForNGS/myPerlScripts/transcriptModel/TSSMotifAnalyzer/v0.2/core_dreme.txt";
my $TATAMotifPath = "/Volumes/A_MPro2TB/softwareForNGS/myPerlScripts/transcriptModel/TSSMotifAnalyzer/v0.2/TATA_dreme.txt";
#<\section>
#====================================================================================================================================================#

#====================================================================================================================================================#
#	section 3_defineOutDirPath
#	primaryDependOnSub: >none
#	secondaryDependOnSub: >none
#
#<section ID="defineOutDirPath" num="3">
my @mkDirAry;
my $resultStorableDir = "$outDir/$resultDirTag/storable/"; push @mkDirAry, $resultStorableDir;
my $resultLogDir = "$outDir/$resultDirTag/log/"; push @mkDirAry, $resultLogDir;
my $resultMotifDir = "$outDir/$resultDirTag/motifUsed/"; push @mkDirAry, $resultMotifDir;
my $resultGFFDir = "$outDir/$resultDirTag/GFF/"; push @mkDirAry, $resultGFFDir;
my $resultWigDir = "$outDir/$resultDirTag/wig/"; push @mkDirAry, $resultWigDir;
my $resultXMLDir = "$outDir/$resultDirTag/XML/"; push @mkDirAry, $resultXMLDir;
my $resultFastaDir = "$outDir/$resultDirTag/fasta/"; push @mkDirAry, $resultFastaDir;
my $mastDir = "$outDir/$resultDirTag/MAST/"; push @mkDirAry, $mastDir;
my $mastBackgroundDir = "$mastDir/background/"; push @mkDirAry, $mastBackgroundDir;
my $mastMotifDir = "$mastDir/motif/"; push @mkDirAry, $mastMotifDir;
my $mastRunDir = "$mastDir/run/"; push @mkDirAry, $mastRunDir;

my $TSSScoreDir = "$outDir/$resultDirTag/TSSScore/$TSSScoreResultTag/"; push @mkDirAry, $TSSScoreDir;
my $TSSScoreStorableDir = "$TSSScoreDir/storable"; push @mkDirAry, $TSSScoreStorableDir;
my $TSSScoreWigDir = "$TSSScoreDir/wig"; push @mkDirAry, $TSSScoreWigDir;
my $TSSScoreLogDir = "$TSSScoreDir/log"; push @mkDirAry, $TSSScoreLogDir;

my $generalggplotDirHsh_ref = {};
my $TSSScoreggplotDirHsh_ref = {};
foreach my $fileType (qw /dat pdf R log/) {
	$generalggplotDirHsh_ref->{$fileType} = "$outDir/$resultDirTag/ggplot/$fileType"; push @mkDirAry, $generalggplotDirHsh_ref->{$fileType};
	$TSSScoreggplotDirHsh_ref->{$fileType} = "$TSSScoreDir/ggplot/$fileType"; push @mkDirAry, $TSSScoreggplotDirHsh_ref->{$fileType};
}

my $weblogoDirHsh_ref = {};
foreach my $fileType (qw /pdf fasta cmd/) {$weblogoDirHsh_ref->{$fileType} = "$outDir/$resultDirTag/weblogo/$fileType"; push @mkDirAry, $weblogoDirHsh_ref->{$fileType};}

system ("mkdir -pm 777 $_") foreach @mkDirAry;
#<\section>
#====================================================================================================================================================#

#====================================================================================================================================================#
#	section 4_defineOutFilePath
#	primaryDependOnSub: >none
#	secondaryDependOnSub: >none
#
#<section ID="defineOutFilePath" num="4">
my $mRNANATInfoLogPath = "$resultLogDir/mRNANATInfoLog.xls";
my $mRNANATInfoHshPath = "$resultStorableDir/mRNANATInfoHsh.pls";
my $bkgdNtFreqHshPath = "$resultStorableDir/bkgdNtFreqHsh.pls";
my $InrMotifFilePath = "$mastMotifDir/InrMotifFile.txt";
my $motifFilePathHsh_ref = {};
$motifFilePathHsh_ref->{'Inr'} = $InrMotifFilePath;
$motifFilePathHsh_ref->{'core'} = $coreMotifPath;
$motifFilePathHsh_ref->{'TATA'} = $TATAMotifPath;

system ("cp -f $coreMotifPath $resultMotifDir/core.txt");
system ("cp -f $TATAMotifPath $resultMotifDir/TATA.txt");

#<\section>
#====================================================================================================================================================#

#====================================================================================================================================================#
#	section 5_processInputData
#	primaryDependOnSub: checkGeneInfo|496, generateFastaLengthHsh|981, getCtgryGeneInfo|1246, getIndivCntgCovPlsPath|1331, readGFF_oneRNAPerGene|2945, readMultiFasta|3044, reverseComplementRefFasta|3158, zipUnzipCntgCovInPlsPathHsh|3349
#	secondaryDependOnSub: currentTime|694, reportStatus|3137
#
#<section ID="processInputData" num="5">

#----------Read fasta
my $fastaHsh_ref = &readMultiFasta($fastaPath);#->3044
my ($fastaWithRevComHsh_ref, $revComFastaPath) = &reverseComplementRefFasta($fastaHsh_ref, $resultFastaDir);#->3158
my ($fastaLengthHsh_ref) = &generateFastaLengthHsh($fastaHsh_ref);#->981

#----------Read Gff
my ($geneInfoHsh_ref) = &readGFF_oneRNAPerGene($gffPath);#->2945
&checkGeneInfo($geneInfoHsh_ref);#->496

#----------Get bfmRNA and mRNA ranges
my @mRNAAry = qw/bfmRNA/;
my ($mRNAInfoHsh_ref, $mRNAByCntgHsh_ref)= &getCtgryGeneInfo($geneInfoHsh_ref, \@mRNAAry);#->1246

#----------Get the paths of the fullReadDefineNATPileup
my $fullReadDefineNATPileupPathHsh_ref = &getIndivCntgCovPlsPath($fullReadDefineNATPileupIndxPath);#->1331
&zipUnzipCntgCovInPlsPathHsh('unzip', $fullReadDefineNATPileupPathHsh_ref);#->3349
#<\section>
#====================================================================================================================================================#

#====================================================================================================================================================#
#	section 6_defineGeneWithEnd3NAT
#	primaryDependOnSub: defineEnd3NATInmRNA|746, plotmRNANATInfo|2552, printmRNANATInfoLog|2908
#	secondaryDependOnSub: ggplotCulmulativeFrequency|1977, reportStatus|3137
#
#<section ID="defineGeneWithEnd3NAT" num="6">
my ($mRNANATInfoHsh_ref) = &defineEnd3NATInmRNA($mRNAInfoHsh_ref, $mRNAByCntgHsh_ref, $fullReadDefineNATPileupPathHsh_ref, $end3NATDefineParamHsh_ref, $mRNANATInfoHshPath, $forceDefineEnd3NAT);#->746
&plotmRNANATInfo($mRNANATInfoHsh_ref, $generalggplotDirHsh_ref);#->2552
&printmRNANATInfoLog($mRNANATInfoHsh_ref, $mRNANATInfoLogPath);#->2908
#<\section>
#====================================================================================================================================================#

#====================================================================================================================================================#
#	section 7_retrieveSequenceSurroundingPredefinedTSS
#	primaryDependOnSub: calculateBackgroundNucleotideFrequency|402, generateShuffleSeq|1045, getSequenceAroundmRNAReferencePoint|1589, getmRNAReferencePoints|1943
#	secondaryDependOnSub: generateMASTBackgroundFile|1000, getBaseAtTSSAndExon|1182, getGeneWithValidTSSAtGeneEnd|1283, getIndivCntgCovPlsPath|1331, investigateTSSRelativeToATGAndTAA|2169, reportStatus|3137
#
#<section ID="retrieveSequenceSurroundingPredefinedTSS" num="7">
my ($mRNARefPtHsh_ref) = &getmRNAReferencePoints($mRNAInfoHsh_ref, $validTSSCntgPlsPath);#->1943

my ($seqAroundSiteHsh_ref, $seqAroundSiteInfoHsh_ref) = &getSequenceAroundmRNAReferencePoint($mRNAInfoHsh_ref, $mRNARefPtHsh_ref, $fastaHsh_ref, $siteSeqRng, $resultFastaDir, $resultStorableDir);#->1589

my ($shuffleSeqAroundSiteHsh_ref, $shuffleSeqAroundSiteInfoHsh_ref) = &generateShuffleSeq($seqAroundSiteHsh_ref, $resultFastaDir, $resultStorableDir);#->1045

my ($bkgdNtFreqHsh_ref) = &calculateBackgroundNucleotideFrequency($seqAroundSiteHsh_ref, $maxPolymerSize, $mastBackgroundDir, $bkgdNtFreqHshPath);#->402
#<\section>
#====================================================================================================================================================#

#====================================================================================================================================================#
#	section 8_getInitiatorMotif
#	primaryDependOnSub: calculateBaseCompositionInAlignments|448, createInrMotifFile|593, defineDegenerateNucleotideHsh|712, getInrMotif|1364, plotInrMotifProportion|2312
#	secondaryDependOnSub: createWeblogo|667, ggplotXYLinesMultipleSamples|2134, reportStatus|3137
#
#<section ID="getInitiatorMotif" num="8">
my ($degenNtHsh_ref) = &defineDegenerateNucleotideHsh();#->712

my ($InrMotifFreqHsh_ref, $InrSeqAlignHsh_ref) = &getInrMotif($seqAroundSiteHsh_ref, $weblogoDirHsh_ref);#->1364

&plotInrMotifProportion($InrMotifFreqHsh_ref, $generalggplotDirHsh_ref, $cmltvePprtnLimit);#->2312

my (undef, $InrBaseComHsh_ref) = &calculateBaseCompositionInAlignments($InrSeqAlignHsh_ref);#->448

&createInrMotifFile($bkgdNtFreqHsh_ref, $InrBaseComHsh_ref, $InrMotifFilePath);#->593
#<\section>
#====================================================================================================================================================#

#====================================================================================================================================================#
#	section 9_scanMotifOccurence
#	primaryDependOnSub: scanMotifAroundSiteWithMAST|3210, scanMotifWholeGenomeWithMAST|3277
#	secondaryDependOnSub: createMotifHitHshStorable|628, generateMASTBackgroundFile|1000, getMastGenomeBothStrandHit|1415, getSingleMotifMASTLogPostionalData|1657, reportStatus|3137, runMAST|3189, storeMotifHitToHshStorable|3324
#
#<section ID="scanMotifOccurence" num="9">
my ($cntgMotifHitPlsPathHsh_ref) = &scanMotifWholeGenomeWithMAST($revComFastaPath, $fastaHsh_ref, $motifFilePathHsh_ref, $mastRunDir, $maxPolymerSize, $mastBackgroundDir, $fastaLengthHsh_ref, $resultStorableDir);#->3277

my ($mastAroundSiteResultHsh_ref) = &scanMotifAroundSiteWithMAST($seqAroundSiteInfoHsh_ref, $shuffleSeqAroundSiteInfoHsh_ref, $motifFilePathHsh_ref, $bkgdNtFreqHsh_ref, $mastRunDir, $predictTSSMotifInfoHsh_ref, $resultStorableDir);#->3210

#<\section>
#====================================================================================================================================================#

#====================================================================================================================================================#
#	section 10_defineMotifBoundary
#	primaryDependOnSub: defineMotifPositionBoundsUsingShuffleBackground|926, plotMotifPostionFactor|2354
#	secondaryDependOnSub: getMotifBoundCutoff|1452, ggplotXYLinesMultipleSamples|2134
#
#<section ID="defineMotifBoundary" num="10">
my ($motifPostionFactorHsh_ref) = &defineMotifPositionBoundsUsingShuffleBackground($predictTSSMotifInfoHsh_ref, $mastAroundSiteResultHsh_ref, $generalggplotDirHsh_ref, $minStreakHsh_ref, $abvBkgrdFactor, $extraMargin, $siteSeqRng);#->926
&plotMotifPostionFactor($motifPostionFactorHsh_ref, $generalggplotDirHsh_ref);#->2354
#<\section>
#====================================================================================================================================================#

#====================================================================================================================================================#
#	section 11_destroyBigUnusedVar
#	primaryDependOnSub: >none
#	secondaryDependOnSub: >none
#
#<section ID="destroyBigUnusedVar" num="11">
undef $mastAroundSiteResultHsh_ref;
undef $seqAroundSiteHsh_ref;
#<\section>
#====================================================================================================================================================#

#====================================================================================================================================================#
#	section 12_predictTSS
#	primaryDependOnSub: predictGenomeWideTSS|2591, printTSSScoreWiggle|2834
#	secondaryDependOnSub: checkRunningThreadAndWaitToJoin|522, createEmptyStorableForGenowideTSSPredictionData|550, generateThreadHshWithRandomCntg|1107, printWiggleSingleTrackFromCntgCovPlsPathHsh|2870, reportStatus|3137
#
#<section ID="predictTSS" num="12">
my ($genomeWideTSSPlsPathHsh_ref) = &predictGenomeWideTSS($motifPostionFactorHsh_ref, $cntgMotifHitPlsPathHsh_ref, $predictTSSMotifInfoHsh_ref, $fastaHsh_ref, $TSSScoreStorableDir, $logTransformPVal);#->2591
&printTSSScoreWiggle($TSSScoreWigDir, $genomeWideTSSPlsPathHsh_ref);#->2834
#<\section>
#====================================================================================================================================================#

#====================================================================================================================================================#
#	section 13_analyzeTSSScore
#	primaryDependOnSub: getAverageTSSScoreForEnd3NATGenes|1134, getProportionOfGenesWithTSSScoreAboveCutOff|1530, getTSSScoreAroundmRNAReferencePoint|1692, getTSSScoreForAllGenes|1784, getTSSScoreInExonAndTSS|1884, plotAverageTSSScoreForEnd3NATGenes|2269, plotTSSScoreAroundmRNAReferencePoint|2391, plotTSSScoreInExonAndTSS|2435, plotTSSScoreInGeneWithAndWithoutEnd3NAT|2493, printGeneAboveTSSScoreLog|2770
#	secondaryDependOnSub: ggplotMultiSampleBoxWhisker|2016, ggplotTwoSampleHistogram|2060, ggplotXYLinesMultipleSamples|2134, reportStatus|3137
#
#<section ID="analyzeTSSScore" num="13">
my ($TSSScoreExonTSSHshNonZero_ref, $TSSScoreExonTSSHshWithZero_ref) = &getTSSScoreInExonAndTSS($genomeWideTSSPlsPathHsh_ref, $mRNAByCntgHsh_ref, $mRNARefPtHsh_ref, $siteSeqRng, $mRNAInfoHsh_ref);#->1884
&plotTSSScoreInExonAndTSS($TSSScoreExonTSSHshNonZero_ref, $TSSScoreExonTSSHshWithZero_ref, $TSSScoreggplotDirHsh_ref);#->2435

my ($TSSScoreRefPtPlotHsh_ref, $TSSScoreRefPtBymRNAHsh_ref) = &getTSSScoreAroundmRNAReferencePoint($mRNAInfoHsh_ref, $mRNARefPtHsh_ref, $genomeWideTSSPlsPathHsh_ref, $siteSeqRng, $mRNAByCntgHsh_ref);#->1692
&plotTSSScoreAroundmRNAReferencePoint($TSSScoreRefPtPlotHsh_ref, $TSSScoreggplotDirHsh_ref, $mRNARefPtHsh_ref);#->2391

my ($TSSScoreIndivmRNAHsh_ref) = &getTSSScoreForAllGenes($mRNAInfoHsh_ref, $genomeWideTSSPlsPathHsh_ref, $mRNAByCntgHsh_ref, $skipEndNt, $TSScoreRegSizeHsh_ref, $fastaLengthHsh_ref, $TSSScoreStorableDir, $forceCalculateTSSScoreIndivmRNA);#->1784

my ($end3NATGeneAvgTSSScoreHsh_ref) = &getAverageTSSScoreForEnd3NATGenes($TSSScoreIndivmRNAHsh_ref, $mRNANATInfoHsh_ref, $mRNARefPtHsh_ref);#->1134

&plotAverageTSSScoreForEnd3NATGenes($end3NATGeneAvgTSSScoreHsh_ref, $TSSScoreggplotDirHsh_ref);#->2269

my ($geneAboveTSSScoreHsh_ref) = &getProportionOfGenesWithTSSScoreAboveCutOff($TSSScoreIndivmRNAHsh_ref, $mRNANATInfoHsh_ref, $mRNARefPtHsh_ref);#->1530

&printGeneAboveTSSScoreLog($geneAboveTSSScoreHsh_ref, $TSSScoreLogDir);#->2770

&plotTSSScoreInGeneWithAndWithoutEnd3NAT($TSSScoreIndivmRNAHsh_ref, $TSSScoreggplotDirHsh_ref, $mRNARefPtHsh_ref, $mRNANATInfoHsh_ref);#->2493

#<\section>
#====================================================================================================================================================#

#====================================================================================================================================================#
#	section 14_finishingTasks
#	primaryDependOnSub: printCMDLogOrFinishMessage|2737, zipUnzipCntgCovInPlsPathHsh|3349
#	secondaryDependOnSub: currentTime|694, reportStatus|3137
#
#<section ID="finishingTasks" num="14">
#-----open the XML in finder
#&zipUnzipCntgCovInPlsPathHsh('zip', $fullReadDefineNATPileupPathHsh_ref);

#-----print the log
&printCMDLogOrFinishMessage("finishMessage");#->2737
system ("open $TSSScoreggplotDirHsh_ref->{'pdf'}");
close TMPLOG;
#<\section>
#====================================================================================================================================================#

#====================================================================================================================================================#
}	#Main sections lexical scope ends
#====================================================================================================================================================#

#====================================================================================================================================================#
#List of subroutines by category
#
#	alignment [n=1]:
#		calculateBaseCompositionInAlignments
#
#	baseComposition [n=2]:
#		calculateBaseCompositionInAlignments, defineDegenerateNucleotideHsh
#
#	fasta [n=1]:
#		readMultiFasta
#
#	general [n=8]:
#		checkGeneInfo, currentTime, getCtgryGeneInfo
#		printCMDLogOrFinishMessage, readGFF_oneRNAPerGene, readMultiFasta
#		readParameters, reportStatus
#
#	gff [n=3]:
#		checkGeneInfo, getCtgryGeneInfo, readGFF_oneRNAPerGene
#
#	ggplot [n=2]:
#		ggplotTwoSampleHistogram, ggplotXYLinesMultipleSamples
#
#	multithread [n=2]:
#		checkRunningThreadAndWaitToJoin, generateThreadHshWithRandomCntg
#
#	plotInR [n=2]:
#		ggplotTwoSampleHistogram, ggplotXYLinesMultipleSamples
#
#	reporting [n=1]:
#		currentTime
#
#	specific [n=1]:
#		getmRNAReferencePoints
#
#	storable [n=2]:
#		getIndivCntgCovPlsPath, zipUnzipCntgCovInPlsPathHsh
#
#	thridPartyApp [n=2]:
#		createWeblogo, runMAST
#
#	unassigned [n=42]:
#		calculateBackgroundNucleotideFrequency, createEmptyStorableForGenowideTSSPredictionData, createInrMotifFile
#		createMotifHitHshStorable, defineEnd3NATInmRNA, defineHardCodedParameters
#		defineMotifPositionBoundsUsingShuffleBackground, generateFastaLengthHsh, generateMASTBackgroundFile
#		generateShuffleSeq, getAverageTSSScoreForEnd3NATGenes, getBaseAtTSSAndExon
#		getGeneWithValidTSSAtGeneEnd, getInrMotif, getMastGenomeBothStrandHit
#		getMotifBoundCutoff, getProportionOfGenesWithTSSScoreAboveCutOff, getSequenceAroundmRNAReferencePoint
#		getSingleMotifMASTLogPostionalData, getTSSScoreAroundmRNAReferencePoint, getTSSScoreForAllGenes
#		getTSSScoreInExonAndTSS, ggplotCulmulativeFrequency, ggplotMultiSampleBoxWhisker
#		investigateTSSRelativeToATGAndTAA, plotAverageTSSScoreForEnd3NATGenes, plotInrMotifProportion
#		plotMotifPostionFactor, plotTSSScoreAroundmRNAReferencePoint, plotTSSScoreInExonAndTSS
#		plotTSSScoreInGeneWithAndWithoutEnd3NAT, plotmRNANATInfo, predictGenomeWideTSS
#		printAllMotifWiggle, printGeneAboveTSSScoreLog, printMotifWiggle
#		printTSSScoreWiggle, printmRNANATInfoLog, reverseComplementRefFasta
#		scanMotifAroundSiteWithMAST, scanMotifWholeGenomeWithMAST, storeMotifHitToHshStorable
#
#	wiggle [n=1]:
#		printWiggleSingleTrackFromCntgCovPlsPathHsh
#
#====================================================================================================================================================#

sub calculateBackgroundNucleotideFrequency {
#....................................................................................................................................................#
#	subroutineCategory: unassigned
#	dependOnSub: generateMASTBackgroundFile|1000, reportStatus|3137
#	appearInSub: >none
#	primaryAppearInSection: 7_retrieveSequenceSurroundingPredefinedTSS|212
#	secondaryAppearInSection: >none
#	input: $bkgdNtFreqHshPath, $mastBackgroundDir, $maxPolymerSize, $seqAroundSiteHsh_ref
#	output: $bkgdNtFreqHsh_ref
#	toCall: my ($bkgdNtFreqHsh_ref) = &calculateBackgroundNucleotideFrequency($seqAroundSiteHsh_ref, $maxPolymerSize, $mastBackgroundDir, $bkgdNtFreqHshPath);
#	calledInLine: 223
#....................................................................................................................................................#
	my ($seqAroundSiteHsh_ref, $maxPolymerSize, $mastBackgroundDir, $bkgdNtFreqHshPath) = @_;
	
	my $bkgdNtFreqHsh_ref = {};
	
	if (-s $bkgdNtFreqHshPath) {
	
		$bkgdNtFreqHsh_ref = retrieve($bkgdNtFreqHshPath);
		&reportStatus("Retrieving nucleotide frequencies around all siteTypes", 0,"\n");#->3137
	
	} else {
		foreach my $siteType (keys %{$seqAroundSiteHsh_ref}) {
			my $siteTypeFullSeqHsh_ref = {};
			&reportStatus("Getting nucleotide frequencies around $siteType", 0,"\n");#->3137
		
			foreach my $mRNAID (keys %{$seqAroundSiteHsh_ref->{$siteType}}) {
				$siteTypeFullSeqHsh_ref->{'full'}{$mRNAID} = $seqAroundSiteHsh_ref->{$siteType}{$mRNAID}{'upStrm'}.$seqAroundSiteHsh_ref->{$siteType}{$mRNAID}{'dnStrm'};
				$siteTypeFullSeqHsh_ref->{'upStrm'}{$mRNAID} = $seqAroundSiteHsh_ref->{$siteType}{$mRNAID}{'upStrm'};
				$siteTypeFullSeqHsh_ref->{'dnStrm'}{$mRNAID} = $seqAroundSiteHsh_ref->{$siteType}{$mRNAID}{'dnStrm'};
			}

			foreach my $fullOrUpStrmOrDnStrm (keys %{$siteTypeFullSeqHsh_ref}) {
				my $seqHsh_ref = \%{$siteTypeFullSeqHsh_ref->{$fullOrUpStrmOrDnStrm}};
				my $bfilePath = "$mastBackgroundDir/$siteType.$fullOrUpStrmOrDnStrm.mast.bfile.txt";
				my ($freqHsh_ref, $proportionHsh_ref) = &generateMASTBackgroundFile($seqHsh_ref, $bfilePath, $maxPolymerSize);#->1000
				$bkgdNtFreqHsh_ref->{$siteType}{$fullOrUpStrmOrDnStrm}{'proportionHsh_ref'} = $proportionHsh_ref;
				$bkgdNtFreqHsh_ref->{$siteType}{$fullOrUpStrmOrDnStrm}{'freqHsh_ref'} = $freqHsh_ref;
				$bkgdNtFreqHsh_ref->{$siteType}{$fullOrUpStrmOrDnStrm}{'bfilePath'} = $bfilePath;
			}
		}
		store($bkgdNtFreqHsh_ref, $bkgdNtFreqHshPath);
	}

	return ($bkgdNtFreqHsh_ref);
}
sub calculateBaseCompositionInAlignments {
#....................................................................................................................................................#
#	subroutineCategory: alignment, baseComposition
#	dependOnSub: >none
#	appearInSub: >none
#	primaryAppearInSection: 8_getInitiatorMotif|228
#	secondaryAppearInSection: >none
#	input: $seqAlignHsh_ref
#	output: $baseCompByBaseHsh_ref, $baseCompByPosHsh_ref
#	toCall: my ($baseCompByBaseHsh_ref, $baseCompByPosHsh_ref) = &calculateBaseCompositionInAlignments($seqAlignHsh_ref);
#	calledInLine: 239
#....................................................................................................................................................#

	my ($seqAlignHsh_ref) = @_;
	
	my $baseCountByPosHsh_ref = {};
	my $baseCompByBaseHsh_ref = {};
	my $baseCompByPosHsh_ref = {};
	my $validSeqNum = 0;
	my $tmpLengthHsh_ref = {};
	
	foreach my $seqName (keys %{$seqAlignHsh_ref}) {
		next if $seqAlignHsh_ref->{$seqName} =~ m/[^ATGCatgc]/;
		$validSeqNum++;
		my @seqAry = split //, $seqAlignHsh_ref->{$seqName};
		$tmpLengthHsh_ref->{@seqAry}++;
		for my $pos (0..$#seqAry) {
			my $base = $seqAry[$pos];
			$base =~ tr/atgc/ATGC/;
			$baseCountByPosHsh_ref->{$pos}{$base}++;
		}
	}
	
	my $lengthNum = keys %{$tmpLengthHsh_ref};
	die 'Length of the sequences in the alignment is not uniform\n' if $lengthNum > 1;
	
	foreach my $base (qw/A T G C/) {
		foreach my $pos (sort {$a <=> $b} keys %{$baseCountByPosHsh_ref}) {
			$baseCountByPosHsh_ref->{$pos}{$base} = 0 if not $baseCountByPosHsh_ref->{$pos}{$base};
			my $proportion = $baseCountByPosHsh_ref->{$pos}{$base}/$validSeqNum;
			$baseCompByBaseHsh_ref->{$base}{$pos} = $proportion;
			$baseCompByPosHsh_ref->{$pos}{$base} = $proportion;
		}
	}
	
	return ($baseCompByBaseHsh_ref, $baseCompByPosHsh_ref);
	
}
sub checkGeneInfo {
#....................................................................................................................................................#
#	subroutineCategory: general, gff
#	dependOnSub: reportStatus|3137
#	appearInSub: >none
#	primaryAppearInSection: 5_processInputData|174
#	secondaryAppearInSection: >none
#	input: $geneInfoHsh_ref
#	output: none
#	toCall: &checkGeneInfo($geneInfoHsh_ref);
#	calledInLine: 187
#....................................................................................................................................................#
	
	my ($geneInfoHsh_ref) = @_;
	
	&reportStatus("Checking gene categories", 0, "\n");#->3137
	my $ctrgyCountHsh_ref = {};
	foreach my $geneID (keys %{$geneInfoHsh_ref}) {
		my $ctgry = $geneInfoHsh_ref->{$geneID}{'ctgry'};
		$ctrgyCountHsh_ref->{$ctgry}++;
	}
	
	foreach my $ctgry (sort keys %{$ctrgyCountHsh_ref}) {
		&reportStatus("Item in $ctgry = $ctrgyCountHsh_ref->{$ctgry}", 0, "\n");#->3137
	}
}
sub checkRunningThreadAndWaitToJoin {
#....................................................................................................................................................#
#	subroutineCategory: multithread
#	dependOnSub: reportStatus|3137
#	appearInSub: predictGenomeWideTSS|2591, printTSSScoreWiggle|2834
#	primaryAppearInSection: >none
#	secondaryAppearInSection: 12_predictTSS|282
#	input: $sleepTime, $verbose
#	output: none
#	toCall: &checkRunningThreadAndWaitToJoin($verbose, $sleepTime);
#	calledInLine: 2701, 2865
#....................................................................................................................................................#
	
	my ($verbose, $sleepTime) = @_;
	
	my @runningThrAry = threads->list(threads::running);
	my @joinableThrAry = threads->list(threads::joinable);
	while (@runningThrAry or @joinableThrAry) {
		@runningThrAry = threads->list(threads::running);
		@joinableThrAry = threads->list(threads::joinable);
		foreach my $joinableThr (@joinableThrAry) {
			$joinableThr->detach() if not $joinableThr->is_running();
		}
		my $numThreadRunning = scalar @runningThrAry;
		&reportStatus("The last $numThreadRunning threads are still running", 20, "\r") if $verbose eq 'yes';#->3137
		sleep $sleepTime;
	}
}
sub createEmptyStorableForGenowideTSSPredictionData {
#....................................................................................................................................................#
#	subroutineCategory: unassigned
#	dependOnSub: reportStatus|3137
#	appearInSub: predictGenomeWideTSS|2591
#	primaryAppearInSection: >none
#	secondaryAppearInSection: 12_predictTSS|282
#	input: $fastaHsh_ref, $resultStorableDir
#	output: $allGenomeWideTSSPlsPathExist, $genomeWideTSSPlsPathHsh_ref
#	toCall: my ($genomeWideTSSPlsPathHsh_ref, $allGenomeWideTSSPlsPathExist) = &createEmptyStorableForGenowideTSSPredictionData($fastaHsh_ref, $resultStorableDir);
#	calledInLine: 2604
#....................................................................................................................................................#
	my ($fastaHsh_ref, $resultStorableDir) = @_;
	
	my $genomeWideTSSPlsDir = "$resultStorableDir/genomeWideTSSPrediction/";
	system ("mkdir -pm 777 $genomeWideTSSPlsDir");

	my $genomeWideTSSPlsIdxHsh_ref = {};
	my $genomeWideTSSPlsPathHsh_ref = {};
	my $allGenomeWideTSSPlsPathExist = 'yes';
	
	foreach my $cntg (keys %{$fastaHsh_ref}) {
		my $cntgPlsName = "$cntg.ary.pls";
		my $cntgPlsPath = "$genomeWideTSSPlsDir/$cntgPlsName";
		$genomeWideTSSPlsIdxHsh_ref->{$cntg} = $cntgPlsName;
		$genomeWideTSSPlsPathHsh_ref->{$cntg} = $cntgPlsPath;
		if (not -s $cntgPlsPath) {

			&reportStatus("Creating empty storabe for $cntg", 20, "\r");#->3137

			my $cntgPosAry_ref = ();
			push @{$cntgPosAry_ref}, undef foreach (1..length($fastaHsh_ref->{$cntg}));

			store($cntgPosAry_ref, $cntgPlsPath);
			$allGenomeWideTSSPlsPathExist = 'no';
		}
	}

	my $genomeWideTSSPlsIdxPath = "$genomeWideTSSPlsDir/index.hsh.pls";
	store($genomeWideTSSPlsIdxHsh_ref, $genomeWideTSSPlsIdxPath);

	return ($genomeWideTSSPlsPathHsh_ref, $allGenomeWideTSSPlsPathExist);
}
sub createInrMotifFile {
#....................................................................................................................................................#
#	subroutineCategory: unassigned
#	dependOnSub: >none
#	appearInSub: >none
#	primaryAppearInSection: 8_getInitiatorMotif|228
#	secondaryAppearInSection: >none
#	input: $InrBaseComHsh_ref, $InrMotifFilePath, $bkgdNtFreqHsh_ref
#	output: 
#	toCall: &createInrMotifFile($bkgdNtFreqHsh_ref, $InrBaseComHsh_ref, $InrMotifFilePath);
#	calledInLine: 241
#....................................................................................................................................................#
	my ($bkgdNtFreqHsh_ref, $InrBaseComHsh_ref, $InrMotifFilePath) = @_;
	
	my $proportionHsh_ref = $bkgdNtFreqHsh_ref->{"mRNA_TSS"}{"full"}{'proportionHsh_ref'};
	my $totalPos = keys %{$InrBaseComHsh_ref};
	
	open INRMOTIF, ">", $InrMotifFilePath;
	print INRMOTIF "MEME version 4\n";
	print INRMOTIF "\n";
	print INRMOTIF "ALPHABET= ACGT\n";
	print INRMOTIF "\n";
	print INRMOTIF "strands: +\n";
	print INRMOTIF "\n";
	print INRMOTIF "Background letter frequencies\n";
	print INRMOTIF "A $proportionHsh_ref->{1}{A} C $proportionHsh_ref->{1}{C} G $proportionHsh_ref->{1}{G} T $proportionHsh_ref->{1}{T}\n";
	print INRMOTIF "\n";
	print INRMOTIF "MOTIF Initiator\n";
	print INRMOTIF "letter-probability matrix: alength=4 w=$totalPos nsites=1000 E=0\n";
	foreach my $pos (sort {$a <=> $b} keys %{$InrBaseComHsh_ref}) {
		print INRMOTIF join " ", ($InrBaseComHsh_ref->{$pos}{'A'}, $InrBaseComHsh_ref->{$pos}{'C'}, $InrBaseComHsh_ref->{$pos}{'G'}, $InrBaseComHsh_ref->{$pos}{'T'}."\n");
	}
	close INRMOTIF;
	return ();
}
sub createMotifHitHshStorable {
#....................................................................................................................................................#
#	subroutineCategory: unassigned
#	dependOnSub: reportStatus|3137
#	appearInSub: scanMotifWholeGenomeWithMAST|3277
#	primaryAppearInSection: >none
#	secondaryAppearInSection: 9_scanMotifOccurence|246
#	input: $fastaHsh_ref, $resultStorableDir
#	output: $allCntgMotifHitPlsExist, $cntgMotifHitIdxPlsPath, $cntgMotifHitPlsPathHsh_ref
#	toCall: my ($cntgMotifHitIdxPlsPath, $cntgMotifHitPlsPathHsh_ref, $allCntgMotifHitPlsExist) = &createMotifHitHshStorable($fastaHsh_ref, $resultStorableDir);
#	calledInLine: 3290
#....................................................................................................................................................#
	my ($fastaHsh_ref, $resultStorableDir) = @_;
	
	my $cntgMotifHitHshPlsDir = "$resultStorableDir/motifHitHsh/";
	my $cntgMotifHitIdxPlsPath = "$cntgMotifHitHshPlsDir/index.hsh.pls";
	my $allCntgMotifHitPlsExist = 'yes';
	
	&reportStatus("Creating empty motif hit hsh storable for motifHits", 0,"\n");#->3137

	system ("mkdir -pm 777 $cntgMotifHitHshPlsDir");
	my $cntgMotifHitPlsIdxHsh_ref = {};
	my $cntgMotifHitPlsPathHsh_ref = {};
	foreach my $cntg (keys %{$fastaHsh_ref}) {
		my $cntgMotifHitPlsName = "$cntg.motifHitHsh.pls";
		my $cntgMotifHitPlsPath = "$cntgMotifHitHshPlsDir/$cntgMotifHitPlsName";
		$cntgMotifHitPlsPathHsh_ref->{$cntg} = $cntgMotifHitPlsPath;
		$cntgMotifHitPlsIdxHsh_ref->{$cntg} = $cntgMotifHitPlsName;

		if (not -s $cntgMotifHitPlsPath) {
			my $cntgMotifHitHsh_ref = {};
			store($cntgMotifHitHsh_ref, $cntgMotifHitPlsPath);
			$allCntgMotifHitPlsExist = 'no';
		}
	}
	store($cntgMotifHitPlsIdxHsh_ref, $cntgMotifHitIdxPlsPath);
	
	return ($cntgMotifHitIdxPlsPath, $cntgMotifHitPlsPathHsh_ref, $allCntgMotifHitPlsExist);
}
sub createWeblogo {
#....................................................................................................................................................#
#	subroutineCategory: thridPartyApp
#	dependOnSub: >none
#	appearInSub: getInrMotif|1364
#	primaryAppearInSection: >none
#	secondaryAppearInSection: 8_getInitiatorMotif|228
#	input: $cmdPath, $fastaPath, $pdfPath, $seqAlignHsh_ref, $seqType, $title
#	output: none
#	toCall: &createWeblogo($seqAlignHsh_ref, $pdfPath, $fastaPath, $cmdPath, $seqType, $title);
#	calledInLine: 1399
#....................................................................................................................................................#

	my ($seqAlignHsh_ref, $pdfPath, $fastaPath, $cmdPath, $seqType, $title) = @_;
	
	#---print the fasta
	open (FASTA, ">", $fastaPath);
	foreach my $seqName (keys %{$seqAlignHsh_ref}) {
		print FASTA ">$seqName\n";
		print FASTA "$seqAlignHsh_ref->{$seqName}\n";
	}
	close FASTA;
	my $cmd = "weblogo --fin $fastaPath --datatype fasta --format pdf --fout $pdfPath --sequence-type $seqType --title \"$title\" --color-scheme classic";
	system ("echo \"$cmd\" >$cmdPath");
	system ($cmd);
	
}
sub currentTime {
#....................................................................................................................................................#
#	subroutineCategory: general, reporting
#	dependOnSub: >none
#	appearInSub: printCMDLogOrFinishMessage|2737, readGFF_oneRNAPerGene|2945, reportStatus|3137
#	primaryAppearInSection: >none
#	secondaryAppearInSection: 0_startingTasks|78, 14_finishingTasks|320, 5_processInputData|174
#	input: none
#	output: $runTime
#	toCall: my ($runTime) = &currentTime();
#	calledInLine: 2757, 2760, 2765, 2965, 3153
#....................................................................................................................................................#
	
	my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst)=localtime(time);
	my $runTime = sprintf "%04d-%02d-%02d %02d:%02d", $year+1900, $mon+1,$mday,$hour,$min;	
	
	return $runTime;
}
sub defineDegenerateNucleotideHsh {
#....................................................................................................................................................#
#	subroutineCategory: baseComposition
#	dependOnSub: >none
#	appearInSub: >none
#	primaryAppearInSection: 8_getInitiatorMotif|228
#	secondaryAppearInSection: >none
#	input: none
#	output: $degenNtHsh_ref
#	toCall: my ($degenNtHsh_ref) = &defineDegenerateNucleotideHsh();
#	calledInLine: 233
#....................................................................................................................................................#

	my $degenNtHsh_ref = {};
	
	$degenNtHsh_ref->{'A'} = 'A';
	$degenNtHsh_ref->{'C'} = 'C';
	$degenNtHsh_ref->{'G'} = 'G';
	$degenNtHsh_ref->{'T'} = 'T';
	$degenNtHsh_ref->{'U'} = 'U';
	$degenNtHsh_ref->{'R'} = 'AG';
	$degenNtHsh_ref->{'Y'} = 'CT';
	$degenNtHsh_ref->{'S'} = 'CG';
	$degenNtHsh_ref->{'W'} = 'AT';
	$degenNtHsh_ref->{'K'} = 'GT';
	$degenNtHsh_ref->{'M'} = 'AC';
	$degenNtHsh_ref->{'B'} = 'CGT';
	$degenNtHsh_ref->{'D'} = 'AGT';
	$degenNtHsh_ref->{'H'} = 'ACT';
	$degenNtHsh_ref->{'V'} = 'ACG';
	$degenNtHsh_ref->{'N'} = 'ACGT';
	
	return $degenNtHsh_ref;
}
sub defineEnd3NATInmRNA {
#....................................................................................................................................................#
#	subroutineCategory: unassigned
#	dependOnSub: reportStatus|3137
#	appearInSub: >none
#	primaryAppearInSection: 6_defineGeneWithEnd3NAT|200
#	secondaryAppearInSection: >none
#	input: $end3NATDefineParamHsh_ref, $forceDefineEnd3NAT, $fullReadDefineNATPileupPathHsh_ref, $mRNAByCntgHsh_ref, $mRNAInfoHsh_ref, $mRNANATInfoHshPath
#	output: $mRNANATInfoHsh_ref
#	toCall: my ($mRNANATInfoHsh_ref) = &defineEnd3NATInmRNA($mRNAInfoHsh_ref, $mRNAByCntgHsh_ref, $fullReadDefineNATPileupPathHsh_ref, $end3NATDefineParamHsh_ref, $mRNANATInfoHshPath, $forceDefineEnd3NAT);
#	calledInLine: 205
#....................................................................................................................................................#
	
	my ($mRNAInfoHsh_ref, $mRNAByCntgHsh_ref, $fullReadDefineNATPileupPathHsh_ref, $end3NATDefineParamHsh_ref, $mRNANATInfoHshPath, $forceDefineEnd3NAT) = @_;
	
	my $mRNANATInfoHsh_ref = {};
	my $mRNAWithNAT = my $mRNAWithoutNAT = 0;

	if (-s $mRNANATInfoHshPath and $forceDefineEnd3NAT eq 'no') {
		&reportStatus("mRNANATInfoHshPath found. Skip defining end3NAT", 0, "\n");#->3137
		
		$mRNANATInfoHsh_ref = retrieve($mRNANATInfoHshPath);
		
		foreach my $mRNAID (keys %{$mRNANATInfoHsh_ref}) {
			$mRNAWithNAT++ if $mRNANATInfoHsh_ref->{$mRNAID}{'end3NAT'} eq 'yes';
			$mRNAWithoutNAT++ if $mRNANATInfoHsh_ref->{$mRNAID}{'end3NAT'} eq 'no';
		}
		
	} else {

		my $upStrmRngEnd3NAT = $end3NATDefineParamHsh_ref->{'upStrmRngEnd3NAT'};
		my $dnStrmRngEnd3NAT = $end3NATDefineParamHsh_ref->{'dnStrmRngEnd3NAT'};
		my $minCovPosHitEnd3NAT = $end3NATDefineParamHsh_ref->{'minCovPosHitEnd3NAT'};
		my $minCovPerNtEnd3NAT = $end3NATDefineParamHsh_ref->{'minCovPerNtEnd3NAT'};
		my $minPosHitPctEnd3NAT = $end3NATDefineParamHsh_ref->{'minPosHitPctEnd3NAT'};
		my $maxCovPerNtEnd3NAT = $end3NATDefineParamHsh_ref->{'maxCovPerNtEnd3NAT'};
		my $maxPosHitPctEnd3NAT = $end3NATDefineParamHsh_ref->{'maxPosHitPctEnd3NAT'};
		
		my %tmpAntisenseStrndHsh = ('+'=>'-', '-'=>'+');
		
		foreach my $cntg (keys %{$fullReadDefineNATPileupPathHsh_ref}) {
			&reportStatus("Calculating antisense coverage at gene 3End in $cntg", 20, "\r");#->3137
			my $cntgCovAry_ref = retrieve($fullReadDefineNATPileupPathHsh_ref->{$cntg});
			foreach my $mRNAID (keys %{$mRNAByCntgHsh_ref->{$cntg}}) {
				#--- get gene info
				my $strnd = $mRNAInfoHsh_ref->{$mRNAID}{'strnd'};
				my @CDSRngAry = sort {$a <=> $b} @{$mRNAInfoHsh_ref->{$mRNAID}{'CDSRng'}};
				my ($CDSStart, $CDSEnd) = ($CDSRngAry[0], $CDSRngAry[-1]);
				my @tmpSrchRngAry = ();
				if ($strnd eq '+') {
					@tmpSrchRngAry = ($CDSEnd-$dnStrmRngEnd3NAT..$CDSEnd+$upStrmRngEnd3NAT);
				} else {
					@tmpSrchRngAry = ($CDSStart-$upStrmRngEnd3NAT..$CDSStart+$dnStrmRngEnd3NAT);
				}
			
				#---define all items as zero
				foreach my $item (qw/covSum numPosHit posHitPct covPerNt/) {
					$mRNANATInfoHsh_ref->{$mRNAID}{$item} = 0;
				}

				#---go thr all search positions
				foreach my $pos (@tmpSrchRngAry) {
					my %tmpCovHsh = ('+'=>0, '-'=>0);
					($tmpCovHsh{'+'}, $tmpCovHsh{'-'}) = split /,/, $cntgCovAry_ref->[$pos-1] if $cntgCovAry_ref->[$pos-1];
					$mRNANATInfoHsh_ref->{$mRNAID}{'covSum'} += $tmpCovHsh{$tmpAntisenseStrndHsh{$strnd}};
					$mRNANATInfoHsh_ref->{$mRNAID}{'numPosHit'}++ if ($tmpCovHsh{$tmpAntisenseStrndHsh{$strnd}} >= $minCovPosHitEnd3NAT);
				}
	
				$mRNANATInfoHsh_ref->{$mRNAID}{'posHitPct'} = 100*$mRNANATInfoHsh_ref->{$mRNAID}{'numPosHit'}/@tmpSrchRngAry;
				$mRNANATInfoHsh_ref->{$mRNAID}{'covPerNt'} = $mRNANATInfoHsh_ref->{$mRNAID}{'covSum'}/@tmpSrchRngAry;
				$mRNANATInfoHsh_ref->{$mRNAID}{'cntg'} = $cntg;
				$mRNANATInfoHsh_ref->{$mRNAID}{'srchStart'} = $tmpSrchRngAry[0];
				$mRNANATInfoHsh_ref->{$mRNAID}{'srchEnd'} = $tmpSrchRngAry[-1];
				$mRNANATInfoHsh_ref->{$mRNAID}{'strnd'} = $strnd;

				#---decide whether the gene has 3'End NAT or not;
				if ($mRNANATInfoHsh_ref->{$mRNAID}{'posHitPct'} >= $minPosHitPctEnd3NAT and $mRNANATInfoHsh_ref->{$mRNAID}{'covPerNt'} >= $minCovPerNtEnd3NAT) {
					$mRNANATInfoHsh_ref->{$mRNAID}{'end3NAT'} = 'yes';
					$mRNAWithNAT++;
				} elsif ($mRNANATInfoHsh_ref->{$mRNAID}{'posHitPct'} <= $maxPosHitPctEnd3NAT and $mRNANATInfoHsh_ref->{$mRNAID}{'covPerNt'} <= $maxCovPerNtEnd3NAT) {
					$mRNANATInfoHsh_ref->{$mRNAID}{'end3NAT'} = 'no';
					$mRNAWithoutNAT++;
				} else {
					$mRNANATInfoHsh_ref->{$mRNAID}{'end3NAT'} = 'middle';
				}
			}
		}
	
		store($mRNANATInfoHsh_ref, $mRNANATInfoHshPath) 
	}
	
	&reportStatus("$mRNAWithNAT mRNA NAT+ and $mRNAWithoutNAT mRNA NAT-", 0, "\n");#->3137
	
	return $mRNANATInfoHsh_ref;
}
sub defineHardCodedParameters {
#....................................................................................................................................................#
#	subroutineCategory: unassigned
#	dependOnSub: >none
#	appearInSub: >none
#	primaryAppearInSection: 1_defineHardCodedParam|92
#	secondaryAppearInSection: >none
#	input: none
#	output: $TSSScoreResultTag, $TSScoreRegSizeHsh_ref, $abvBkgrdFactor, $cmltvePprtnLimit, $end3NATDefineParamHsh_ref, $extraMargin, $forceCalculateTSSScoreIndivmRNA, $forceDefineEnd3NAT, $logTransformPVal, $maxPolymerSize, $minStreakHsh_ref, $predictTSSMotifInfoHsh_ref, $resultDirTag, $siteSeqRng, $skipEndNt
#	toCall: my ($forceDefineEnd3NAT, $siteSeqRng, $cmltvePprtnLimit, $maxPolymerSize, $extraMargin, $abvBkgrdFactor, $end3NATDefineParamHsh_ref, $predictTSSMotifInfoHsh_ref, $minStreakHsh_ref, $resultDirTag, $TSSScoreResultTag, $logTransformPVal, $skipEndNt, $TSScoreRegSizeHsh_ref, $forceCalculateTSSScoreIndivmRNA) = &defineHardCodedParameters();
#	calledInLine: 97
#....................................................................................................................................................#
	
	my $forceDefineEnd3NAT = 'yes'; #---force to redine presence or absence of NAT whenever storable found at mRNANATInfoHshPath
	my $forceCalculateTSSScoreIndivmRNA = 'yes'; 
	my $siteSeqRng = 150;
	my $cmltvePprtnLimit = 0.85;
	my $maxPolymerSize = 6; #---to be used to generate background file in mast
	my $extraMargin = 3; #---extra margin added to the motif boundarys defined
	my $abvBkgrdFactor = 1.2; #----the background will be scaled up using this factor during comparison. The purpose is to indentify the signal thay is truely above background
	my $logTransformPVal = "yes"; #-----use 1/pValue or log(1/pValue);
	#my $logTransformPVal = "no"; #-----use 1/pValue or log(1/pValue);
	
	#---ranges for sampling TSS score form all genes
	my $skipEndNt = 0; #---number of nucleotide to be skipped the positions around ATG and TAA
	my $TSScoreRegSizeHsh_ref = (); #---size of the region for sampling TSS Score around ATG and TAA in all gene
	$TSScoreRegSizeHsh_ref->{'upStrm'} = 40;
	$TSScoreRegSizeHsh_ref->{'dnStrm'} = 0;
	
	my $end3NATDefineParamHsh_ref = {};

	#---3'end range
	$end3NATDefineParamHsh_ref->{'upStrmRngEnd3NAT'} = 0;#---upstream of TAA range (ie outside CDS in UTR3) to be search to define whether there're NAT at 3'end of genes
	$end3NATDefineParamHsh_ref->{'dnStrmRngEnd3NAT'} = 100;#---dnstream of TAA range (ie in CDS) to be search to define whether there're NAT at 3'end of genes
	$end3NATDefineParamHsh_ref->{'minCovPosHitEnd3NAT'} = 1;#---minimum coverage at a postion to be regarded as a postion hit

	#---3'end NAT positive
	$end3NATDefineParamHsh_ref->{'minCovPerNtEnd3NAT'} = 5;#---minimum coverage per base as end3NAT postive
	$end3NATDefineParamHsh_ref->{'minPosHitPctEnd3NAT'} = 80;#---minimum percentage of postion hit in the search rng to be regarded as end3NAT postive
	#---3'end NAT negative
	$end3NATDefineParamHsh_ref->{'maxCovPerNtEnd3NAT'} = 0;#---maximum coverage per base as end3NAT negative
	$end3NATDefineParamHsh_ref->{'maxPosHitPctEnd3NAT'} = 0;#---maximum percentage of postion hit in the search rng to be regarded as end3NAT negative

	my $predictTSSMotifInfoHsh_ref = {};
	$predictTSSMotifInfoHsh_ref->{'Inr'}{'scoreFactor'} = 1; #---the factor use to scale the score of the motif
	$predictTSSMotifInfoHsh_ref->{'core'}{'scoreFactor'} = 1; #---the factor use to scale the score of the motif
	$predictTSSMotifInfoHsh_ref->{'TATA'}{'scoreFactor'} = 1; #---the factor use to scale the score of the motif

	$predictTSSMotifInfoHsh_ref->{'Inr'}{'maxPValGenomePredict'} = 0.05; #---the max pval of motif to be valid
	$predictTSSMotifInfoHsh_ref->{'core'}{'maxPValGenomePredict'} = 0.05; #---the max pval of motif to be valid
	$predictTSSMotifInfoHsh_ref->{'TATA'}{'maxPValGenomePredict'} = 0.0005; #---the max pval of motif to be valid

	$predictTSSMotifInfoHsh_ref->{'Inr'}{'maxPValDefineBound'} = 0.05; #---the max pval of motif to be valid for defining the positional bound
	$predictTSSMotifInfoHsh_ref->{'core'}{'maxPValDefineBound'} = 0.05; #---the max pval of motif to be valid for defining the positional bound
	$predictTSSMotifInfoHsh_ref->{'TATA'}{'maxPValDefineBound'} = 0.0005; #---the max pval of motif to be valid for defining the positional bound

	$predictTSSMotifInfoHsh_ref->{'Inr'}{'mustValid'} = 'yes'; #---the motif must be valid or not
	$predictTSSMotifInfoHsh_ref->{'core'}{'mustValid'} = 'yes'; #---the motif must be valid or not
	$predictTSSMotifInfoHsh_ref->{'TATA'}{'mustValid'} = 'yes'; #---the motif must be valid or not

	my $minStreakHsh_ref = {}; #---minimum nt streak to start and end the motif boundary 
	$minStreakHsh_ref->{'positive'} = 3;
	$minStreakHsh_ref->{'negative'} = 2; #----use 2 because for the TATA box and core there're peaks very close to the main peak, use a small number to prevent merging the peaks

	#---generate resultTag
	my %TSSScoreResultTagHsh = ();
	foreach my $motif (keys %{$predictTSSMotifInfoHsh_ref}) {
		my @tmpAry = ();
		foreach my $param (qw/scoreFactor maxPValGenomePredict mustValid/) {
			push @tmpAry, $predictTSSMotifInfoHsh_ref->{$motif}{$param};
		}
		$TSSScoreResultTagHsh{$motif} = join "_", @tmpAry;
	}
	
	my $TSSScoreResultTag = "Inr.$TSSScoreResultTagHsh{'Inr'}.core.$TSSScoreResultTagHsh{'core'}.TATA.$TSSScoreResultTagHsh{'TATA'}.logPVal.$logTransformPVal";

	my @resultDirTagAry = ();
	foreach my $motif (qw/Inr core TATA/) {
		push @resultDirTagAry, "$motif.$predictTSSMotifInfoHsh_ref->{$motif}{'maxPValDefineBound'}";
	}
	
	my $resultDirTag = join ".", @resultDirTagAry;

	return ($forceDefineEnd3NAT, $siteSeqRng, $cmltvePprtnLimit, $maxPolymerSize, $extraMargin, $abvBkgrdFactor, $end3NATDefineParamHsh_ref, $predictTSSMotifInfoHsh_ref, $minStreakHsh_ref, $resultDirTag, $TSSScoreResultTag, $logTransformPVal, $skipEndNt, $TSScoreRegSizeHsh_ref, $forceCalculateTSSScoreIndivmRNA);
}
sub defineMotifPositionBoundsUsingShuffleBackground {
#....................................................................................................................................................#
#	subroutineCategory: unassigned
#	dependOnSub: getMotifBoundCutoff|1452, ggplotXYLinesMultipleSamples|2134
#	appearInSub: >none
#	primaryAppearInSection: 10_defineMotifBoundary|259
#	secondaryAppearInSection: >none
#	input: $abvBkgrdFactor, $extraMargin, $generalggplotDirHsh_ref, $mastAroundSiteResultHsh_ref, $minStreakHsh_ref, $predictTSSMotifInfoHsh_ref, $siteSeqRng
#	output: $motifPostionFactorHsh_ref
#	toCall: my ($motifPostionFactorHsh_ref) = &defineMotifPositionBoundsUsingShuffleBackground($predictTSSMotifInfoHsh_ref, $mastAroundSiteResultHsh_ref, $generalggplotDirHsh_ref, $minStreakHsh_ref, $abvBkgrdFactor, $extraMargin, $siteSeqRng);
#	calledInLine: 264
#....................................................................................................................................................#
	my ($predictTSSMotifInfoHsh_ref, $mastAroundSiteResultHsh_ref, $generalggplotDirHsh_ref, $minStreakHsh_ref, $abvBkgrdFactor, $extraMargin, $siteSeqRng) = @_;

	my $motifPostionFactorHsh_ref = {};
	foreach my $siteType (keys %{$mastAroundSiteResultHsh_ref}) {
		foreach my $motif (keys %{$mastAroundSiteResultHsh_ref->{$siteType}}) {
			my $bothMotifPctHsh_ref = {};
			my $posPctHsh_ref = {};
			foreach my $queryOrShuffle (keys %{$mastAroundSiteResultHsh_ref->{$siteType}{$motif}}) {
				my $motifPctHsh_ref = $mastAroundSiteResultHsh_ref->{$siteType}{$motif}{$queryOrShuffle}{'motifPctHsh_ref'};#---take the has out of the lexcial scope
				foreach my $pos (keys %{$motifPctHsh_ref}) {
					$bothMotifPctHsh_ref->{$queryOrShuffle}{$pos} = $motifPctHsh_ref->{$pos};
					$posPctHsh_ref->{$pos}{$queryOrShuffle} = $motifPctHsh_ref->{$pos};
				}
			}
			my ($leftBound, $rightBound, $postionFactorHsh_ref) = &getMotifBoundCutoff($posPctHsh_ref, $minStreakHsh_ref, $abvBkgrdFactor, $extraMargin, $siteSeqRng, $motif);#->1452
			if ($siteType eq 'mRNA_TSS') {#---get the bound only if siteType if mRNA_TSS and not the initiator
				$motifPostionFactorHsh_ref->{$motif} = $postionFactorHsh_ref;
			}
			#print $siteType."\t".$motif."\t".$leftBound."\t".$rightBound."\n";
			my @posAry = (sort {$a <=> $b} keys %{$posPctHsh_ref});
			my $maxPos = $posAry[-1];
			my $maxHitPVal = $predictTSSMotifInfoHsh_ref->{$motif}{'maxPValDefineBound'};
			my $nameTag = "mast.$motif.p$maxHitPVal.in.$siteType";
			my $plotDataHsh_ref = $bothMotifPctHsh_ref;
			my $pdfPath = $generalggplotDirHsh_ref->{'pdf'}."/$nameTag.pdf";
			my $dataPath = $generalggplotDirHsh_ref->{'dat'}."/$nameTag.dat";
			my $RScriptPath = $generalggplotDirHsh_ref->{'R'}."/$nameTag.R";
			my $logPath = $generalggplotDirHsh_ref->{'log'}."/$nameTag.log";
			my $xAxis = 'relativePositon';
			my $YAxis = 'proportion';
			my $YVariable = 'motif';
			#my $extraArg = '+ ylim(0, 100)';
			my $extraArg = " + scale_x_continuous(breaks=seq(0, $maxPos, by=10))";
			$extraArg .= " + geom_vline(xintercept=c($leftBound), linetype=\"dotted\") + annotate(\"text\", x=$leftBound, y=0, label=\"$leftBound\", vjust=-0.2, hjust=-0.1, angle=90)";
			$extraArg .= " + geom_vline(xintercept=c($rightBound), linetype=\"dotted\") + annotate(\"text\", x=$rightBound, y=0, label=\"$rightBound\", vjust=-0.2, hjust=-0.1, angle=90)";
			my $height = 6;
			my $width = 14;
			&ggplotXYLinesMultipleSamples($plotDataHsh_ref, $dataPath, $pdfPath, $RScriptPath, $logPath, $xAxis, $YAxis, $YVariable, $extraArg, $height, $width);#->2134
		}
	}

	return ($motifPostionFactorHsh_ref);
}
sub generateFastaLengthHsh {
#....................................................................................................................................................#
#	subroutineCategory: unassigned
#	dependOnSub: >none
#	appearInSub: >none
#	primaryAppearInSection: 5_processInputData|174
#	secondaryAppearInSection: >none
#	input: $fastaHsh_ref
#	output: $fastaLengthHsh_ref
#	toCall: my ($fastaLengthHsh_ref) = &generateFastaLengthHsh($fastaHsh_ref);
#	calledInLine: 183
#....................................................................................................................................................#
	my ($fastaHsh_ref) = @_;

	my $fastaLengthHsh_ref = {};
	$fastaLengthHsh_ref->{$_} = length ($fastaHsh_ref->{$_}) foreach (keys %{$fastaHsh_ref});

	return ($fastaLengthHsh_ref);
}
sub generateMASTBackgroundFile {
#....................................................................................................................................................#
#	subroutineCategory: unassigned
#	dependOnSub: >none
#	appearInSub: calculateBackgroundNucleotideFrequency|402, scanMotifWholeGenomeWithMAST|3277
#	primaryAppearInSection: >none
#	secondaryAppearInSection: 7_retrieveSequenceSurroundingPredefinedTSS|212, 9_scanMotifOccurence|246
#	input: $bfilePath, $maxPolymerSize, $seqHsh_ref
#	output: $freqHsh_ref, $proportionHsh_ref
#	toCall: my ($freqHsh_ref, $proportionHsh_ref) = &generateMASTBackgroundFile($seqHsh_ref, $bfilePath, $maxPolymerSize);
#	calledInLine: 436, 3297
#....................................................................................................................................................#
	my ($seqHsh_ref, $bfilePath, $maxPolymerSize) = @_;
	
	my $freqHsh_ref = {};
	my $proportionHsh_ref = {};

	foreach my $polymerSize (1..$maxPolymerSize) {
		my $regexString = '';
		$regexString .= '{A,C,G,T}' foreach (1..$polymerSize);
		$freqHsh_ref->{$polymerSize}{$_} = 0 foreach (glob $regexString);
		foreach my $ID (keys %{$seqHsh_ref}) {
			my $length = length $seqHsh_ref->{$ID};
			foreach my $pos (0..$length-$polymerSize-1) {
				my $polymer = substr $seqHsh_ref->{$ID}, $pos, $polymerSize;
				next if $polymer =~ m/[^ATGC]/i;
				$freqHsh_ref->{$polymerSize}{$polymer}++;
			}
		}
	}
	
	open (FREQ, ">", "$bfilePath");
	foreach my $polymerSize (sort keys %{$freqHsh_ref}) {
		my $sum = sum values %{$freqHsh_ref->{$polymerSize}};
		foreach my $polymer (sort {$freqHsh_ref->{$polymerSize}{$b} <=> $freqHsh_ref->{$polymerSize}{$a}} keys %{$freqHsh_ref->{$polymerSize}}) {
			my $proportion = sprintf "%.10f", $freqHsh_ref->{$polymerSize}{$polymer}/$sum;
			$proportion = 0.0000000001 if $proportion == 0; #--- stupid MAST doesnt like zero!!!
			$proportionHsh_ref->{$polymerSize}{$polymer} = $proportion;
			print FREQ join "", ((join "\t", ($polymer, $proportion)), "\n");
		}
	}
	close FREQ;
	
	return ($freqHsh_ref, $proportionHsh_ref);
}
sub generateShuffleSeq {
#....................................................................................................................................................#
#	subroutineCategory: unassigned
#	dependOnSub: reportStatus|3137
#	appearInSub: >none
#	primaryAppearInSection: 7_retrieveSequenceSurroundingPredefinedTSS|212
#	secondaryAppearInSection: >none
#	input: $resultFastaDir, $resultStorableDir, $seqAroundSiteHsh_ref
#	output: $shuffleSeqAroundSiteHsh_ref, $shuffleSeqAroundSiteInfoHsh_ref
#	toCall: my ($shuffleSeqAroundSiteHsh_ref, $shuffleSeqAroundSiteInfoHsh_ref) = &generateShuffleSeq($seqAroundSiteHsh_ref, $resultFastaDir, $resultStorableDir);
#	calledInLine: 221
#....................................................................................................................................................#
	my ($seqAroundSiteHsh_ref, $resultFastaDir, $resultStorableDir) = @_;
	
	my $shuffleSeqAroundSiteHsh_ref = {};
	my $shuffleSeqAroundSiteInfoHsh_ref = {};
	
	my $shuffleSeqAroundSiteHshPlsPath = "$resultStorableDir/shuffleSeqAroundSiteHsh.pls";
	my $shuffleSeqAroundSiteInfoHshPlsPath = "$resultStorableDir/shuffleSeqAroundSiteInfoHsh.pls";
	
	if (-s $shuffleSeqAroundSiteHshPlsPath and -s $shuffleSeqAroundSiteInfoHshPlsPath) {

		&reportStatus("shuffleSeqAroundSite storable found. Retrieving", 0, "\n");#->3137

		$shuffleSeqAroundSiteHsh_ref = retrieve($shuffleSeqAroundSiteHshPlsPath);
		$shuffleSeqAroundSiteInfoHsh_ref = retrieve($shuffleSeqAroundSiteInfoHshPlsPath);
	
	} else {

		my @itemAry = qw/upStrm dnStrm full/;
		foreach my $siteType (keys %{$seqAroundSiteHsh_ref}) {
			&reportStatus("Shuffling sequences around $siteType", 0, "\n");#->3137
			my %fastaFHHsh = ();
			foreach my $fullOrUpStrmOrDnStrm (@itemAry) {
				my $fastaPath = "$resultFastaDir/shuffle.$siteType.$fullOrUpStrmOrDnStrm.fasta";
				open $fastaFHHsh{$fullOrUpStrmOrDnStrm}, ">", $fastaPath;
				$shuffleSeqAroundSiteInfoHsh_ref->{$siteType}{$fullOrUpStrmOrDnStrm}{'fastaPath'} = $fastaPath;
				$shuffleSeqAroundSiteInfoHsh_ref->{$siteType}{$fullOrUpStrmOrDnStrm}{'totalSeqNum'} = keys %{$seqAroundSiteHsh_ref->{$siteType}};
			}
			foreach my $mRNAID (keys %{$seqAroundSiteHsh_ref->{$siteType}}) {
				my %tmpAryHsh = ();
				@{$tmpAryHsh{'upStrm'}} = split //, $seqAroundSiteHsh_ref->{$siteType}{$mRNAID}{'upStrm'};
				@{$tmpAryHsh{'dnStrm'}} = split //, $seqAroundSiteHsh_ref->{$siteType}{$mRNAID}{'dnStrm'};
				@{$tmpAryHsh{'upStrm'}} = shuffle(@{$tmpAryHsh{'upStrm'}});
				@{$tmpAryHsh{'dnStrm'}} = shuffle(@{$tmpAryHsh{'dnStrm'}});
				@{$tmpAryHsh{'full'}} = (@{$tmpAryHsh{'upStrm'}}, @{$tmpAryHsh{'dnStrm'}});
				foreach my $fullOrUpStrmOrDnStrm (@itemAry) {
					$shuffleSeqAroundSiteHsh_ref->{$siteType}{$fullOrUpStrmOrDnStrm}{$mRNAID} = join '', @{$tmpAryHsh{$fullOrUpStrmOrDnStrm}};
					$shuffleSeqAroundSiteInfoHsh_ref->{$siteType}{$fullOrUpStrmOrDnStrm}{'length'} = length $shuffleSeqAroundSiteHsh_ref->{$siteType}{$mRNAID}{$fullOrUpStrmOrDnStrm} if not $shuffleSeqAroundSiteInfoHsh_ref->{$siteType}{$fullOrUpStrmOrDnStrm}{'length'};
					print {$fastaFHHsh{$fullOrUpStrmOrDnStrm}} ">$mRNAID\n";
					print {$fastaFHHsh{$fullOrUpStrmOrDnStrm}} $shuffleSeqAroundSiteHsh_ref->{$siteType}{$fullOrUpStrmOrDnStrm}{$mRNAID}."\n";
				}
			}
		}

		store($shuffleSeqAroundSiteHsh_ref, $shuffleSeqAroundSiteHshPlsPath);
		store($shuffleSeqAroundSiteInfoHsh_ref, $shuffleSeqAroundSiteInfoHshPlsPath);
		
	}

	return ($shuffleSeqAroundSiteHsh_ref, $shuffleSeqAroundSiteInfoHsh_ref);
}
sub generateThreadHshWithRandomCntg {
#....................................................................................................................................................#
#	subroutineCategory: multithread
#	dependOnSub: >none
#	appearInSub: predictGenomeWideTSS|2591
#	primaryAppearInSection: >none
#	secondaryAppearInSection: 12_predictTSS|282
#	input: $cntgAry_ref, $threadToSpawn
#	output: $randCntgInThreadHsh_ref
#	toCall: my ($randCntgInThreadHsh_ref) = &generateThreadHshWithRandomCntg($threadToSpawn, $cntgAry_ref);
#	calledInLine: 2621
#....................................................................................................................................................#

	my ($threadToSpawn, $cntgAry_ref) = @_;

	my @shuffleCntgAry = shuffle(@{$cntgAry_ref});
	my $threadNum = 1;
	my $randCntgInThreadHsh_ref = {};
	foreach my $cntg (@{$cntgAry_ref}) {
		$threadNum = 1 if $threadNum > $threadToSpawn;
		push @{$randCntgInThreadHsh_ref->{$threadNum}}, $cntg;
		$threadNum++;
	}
	
	return ($randCntgInThreadHsh_ref);

}
sub getAverageTSSScoreForEnd3NATGenes {
#....................................................................................................................................................#
#	subroutineCategory: unassigned
#	dependOnSub: reportStatus|3137
#	appearInSub: >none
#	primaryAppearInSection: 13_analyzeTSSScore|293
#	secondaryAppearInSection: >none
#	input: $TSSScoreIndivmRNAHsh_ref, $mRNANATInfoHsh_ref, $mRNARefPtHsh_ref
#	output: $end3NATGeneAvgTSSScoreHsh_ref
#	toCall: my ($end3NATGeneAvgTSSScoreHsh_ref) = &getAverageTSSScoreForEnd3NATGenes($TSSScoreIndivmRNAHsh_ref, $mRNANATInfoHsh_ref, $mRNARefPtHsh_ref);
#	calledInLine: 306
#....................................................................................................................................................#
	my ($TSSScoreIndivmRNAHsh_ref, $mRNANATInfoHsh_ref, $mRNARefPtHsh_ref) = @_;
	
	my $end3NATGeneAvgTSSScoreHsh_ref = {};
	
	foreach my $refPosType (keys %{$TSSScoreIndivmRNAHsh_ref}) {

		&reportStatus("Getting average TSSScore in $refPosType", 0, "\n");#->3137
		
		foreach my $mRNAID (keys %{$TSSScoreIndivmRNAHsh_ref->{$refPosType}}) {
			my %tmpAryHsh;
			@{$tmpAryHsh{'a'}} = ();
			@{$tmpAryHsh{'s'}} = ();
			my $numPos = keys %{$TSSScoreIndivmRNAHsh_ref->{$refPosType}{$mRNAID}};
			my $end3NAT = $mRNANATInfoHsh_ref->{$mRNAID}{'end3NAT'};
			foreach my $dirtn (keys %tmpAryHsh) {
				foreach my $rltvPos (keys %{$TSSScoreIndivmRNAHsh_ref->{$refPosType}{$mRNAID}}) {
					push @{$tmpAryHsh{$dirtn}}, $TSSScoreIndivmRNAHsh_ref->{$refPosType}{$mRNAID}{$rltvPos}{$dirtn};
				}
				my $sum = my $max = 0;
				if (@{$tmpAryHsh{$dirtn}} > 0) {
					@{$tmpAryHsh{$dirtn}} = sort {$b <=> $a} @{$tmpAryHsh{$dirtn}};
					$sum = sum(@{$tmpAryHsh{$dirtn}});
					$max = ${$tmpAryHsh{$dirtn}}[0];
				}
				
				if (not exists $mRNARefPtHsh_ref->{'NAT_TSS'}{$mRNAID}) {
					push @{$end3NATGeneAvgTSSScoreHsh_ref->{$refPosType}{$end3NAT}{$dirtn}}, $max;
				} else {
					push @{$end3NATGeneAvgTSSScoreHsh_ref->{$refPosType}{'withNATTSS'}{$dirtn}}, $max;
				}
			}
		}
	}
	
	return ($end3NATGeneAvgTSSScoreHsh_ref);
}
sub getBaseAtTSSAndExon {
#....................................................................................................................................................#
#	subroutineCategory: unassigned
#	dependOnSub: reportStatus|3137
#	appearInSub: getmRNAReferencePoints|1943
#	primaryAppearInSection: >none
#	secondaryAppearInSection: 7_retrieveSequenceSurroundingPredefinedTSS|212
#	input: $mRNAInfoHsh_ref, $strndEndValidTSSInfoHsh_ref
#	output: $mRNARefPtHsh_ref
#	toCall: my ($mRNARefPtHsh_ref) = &getBaseAtTSSAndExon($strndEndValidTSSInfoHsh_ref, $mRNAInfoHsh_ref);
#	calledInLine: 1971
#....................................................................................................................................................#

	my ($strndEndValidTSSInfoHsh_ref, $mRNAInfoHsh_ref) = @_;
	
	my $mRNARefPtHsh_ref = {}; #----collect the mRNA_TSS, NAT_TSS, ATG and TAA as reference point, only the genes with NAT_TSS with collect TAA and only gene with mRNA_TSS with collect ATG
	
	my $siteTypeHsh_ref = {};
	$siteTypeHsh_ref->{'ATG'}{'sense'} = "mRNA_TSS";
	$siteTypeHsh_ref->{'TAA'}{'antisense'} = "NAT_TSS";
	
	&reportStatus("Picking valid TSS postion for base scanning", 40, "\n");#->3137
	
	#---get mRNA with TSS at ATG or TAA and the site with highest rd5Ends
	my $siteByCntgBymRNAHsh_ref = {};
	my $allTSSSiteHash_ref = {};
	foreach my $ATGOrTAA (keys %{$strndEndValidTSSInfoHsh_ref}) {
		foreach my $senseOrAntisense (keys %{$strndEndValidTSSInfoHsh_ref->{$ATGOrTAA}}) {
			my $siteType = $siteTypeHsh_ref->{$ATGOrTAA}{$senseOrAntisense};
			foreach my $mRNAID (keys %{$strndEndValidTSSInfoHsh_ref->{$ATGOrTAA}{$senseOrAntisense}}) {
				#----sort the $rltvPos by cov
				foreach my $rltvPos (sort {$strndEndValidTSSInfoHsh_ref->{$ATGOrTAA}{$senseOrAntisense}{$mRNAID}{$b}{'cov'} <=> $strndEndValidTSSInfoHsh_ref->{$ATGOrTAA}{$senseOrAntisense}{$mRNAID}{$a}{'cov'}} keys %{$strndEndValidTSSInfoHsh_ref->{$ATGOrTAA}{$senseOrAntisense}{$mRNAID}}) {
					$mRNARefPtHsh_ref->{$siteType}{$mRNAID} = $strndEndValidTSSInfoHsh_ref->{$ATGOrTAA}{$senseOrAntisense}{$mRNAID}{$rltvPos}{'absPos'};
					last;
				}
				
				if ($siteType eq 'NAT_TSS') {
					if ($mRNAInfoHsh_ref->{$mRNAID}{'strnd'} eq '-') {
						$mRNARefPtHsh_ref->{'NAT_TSS'}{$mRNAID}--; #---get the outside position as the reference
					} else {
						$mRNARefPtHsh_ref->{'NAT_TSS'}{$mRNAID}++; #---get the outside position as the reference
					}
				}
			}
		}
	}

	foreach my $mRNAID (keys %{$mRNAInfoHsh_ref}) {
		my $strnd = $mRNAInfoHsh_ref->{$mRNAID}{'strnd'};
		my @CDSRngAry = sort {$a <=> $b} @{$mRNAInfoHsh_ref->{$mRNAID}{'CDSRng'}};
		my ($CDSStart, $CDSEnd) = ($CDSRngAry[0], $CDSRngAry[-1]);
		my ($ATGPos, $TAAPos) = ($CDSStart, $CDSEnd);
		if ($strnd eq '-') {
			($ATGPos, $TAAPos) = ($TAAPos, $ATGPos);
			$TAAPos--; #---get the outside position as the reference
		} else {
			$TAAPos++; #---get the outside position as the reference
		}
		$mRNARefPtHsh_ref->{'mRNA_ATG'}{$mRNAID} = $ATGPos if $mRNARefPtHsh_ref->{'mRNA_TSS'}{$mRNAID};
		$mRNARefPtHsh_ref->{'mRNA_TAA'}{$mRNAID} = $TAAPos if $mRNARefPtHsh_ref->{'NAT_TSS'}{$mRNAID};
	}
	
	return ($mRNARefPtHsh_ref);
}
sub getCtgryGeneInfo {
#....................................................................................................................................................#
#	subroutineCategory: general, gff
#	dependOnSub: reportStatus|3137
#	appearInSub: >none
#	primaryAppearInSection: 5_processInputData|174
#	secondaryAppearInSection: >none
#	input: $ctgryAry_ref, $geneInfoHsh_ref
#	output: $geneCtgryByCntgHsh_ref, $geneCtgryInfoHsh_ref
#	toCall: my ($geneCtgryInfoHsh_ref, $geneCtgryByCntgHsh_ref) = &getCtgryGeneInfo($geneInfoHsh_ref, $ctgryAry_ref);
#	calledInLine: 191
#....................................................................................................................................................#

	my ($geneInfoHsh_ref, $ctgryAry_ref) = @_;
	
	my $geneCtgryInfoHsh_ref = {};
	my $geneCtgryByCntgHsh_ref = {};
	
	my $ctgryStr = join ",", @{$ctgryAry_ref};

	&reportStatus("Filtering GFF on cgtry $ctgryStr", 0, "\n");#->3137
	
	foreach my $geneID (keys %{$geneInfoHsh_ref}) {
		my $ctgry = $geneInfoHsh_ref->{$geneID}{'ctgry'};
		my $cntg = $geneInfoHsh_ref->{$geneID}{'cntg'};
		if (grep /^$ctgry$/, @{$ctgryAry_ref}) {
			%{$geneCtgryInfoHsh_ref->{$geneID}} = %{$geneInfoHsh_ref->{$geneID}};
			$geneCtgryByCntgHsh_ref->{$cntg}{$geneID}++;
		}
	}
	
	my $numGene = keys %{$geneCtgryInfoHsh_ref};
	
	&reportStatus("$numGene gene filtered on cgtry $ctgryStr", 0, "\n");#->3137
	
	return $geneCtgryInfoHsh_ref, $geneCtgryByCntgHsh_ref;
}
sub getGeneWithValidTSSAtGeneEnd {
#....................................................................................................................................................#
#	subroutineCategory: unassigned
#	dependOnSub: >none
#	appearInSub: getmRNAReferencePoints|1943
#	primaryAppearInSection: >none
#	secondaryAppearInSection: 7_retrieveSequenceSurroundingPredefinedTSS|212
#	input: $TSSmRNAEndFreqHsh_ref, $geneBasedTSSInfoHsh_ref, $validTSSLimitHsh_ref
#	output: $strndEndValidTSSInfoHsh_ref
#	toCall: my ($strndEndValidTSSInfoHsh_ref) = &getGeneWithValidTSSAtGeneEnd($geneBasedTSSInfoHsh_ref, $TSSmRNAEndFreqHsh_ref, $validTSSLimitHsh_ref);
#	calledInLine: 1970
#....................................................................................................................................................#

	my ($geneBasedTSSInfoHsh_ref, $TSSmRNAEndFreqHsh_ref, $validTSSLimitHsh_ref) = @_;
	#----set the Pct limit and count the number of genes with valid TSS
	
	my $strndEndValidTSSInfoHsh_ref = {};
	foreach my $ATGOrTAA (keys %{$validTSSLimitHsh_ref}) {
		foreach my $senseOrAntisense (keys %{$validTSSLimitHsh_ref->{$ATGOrTAA}}) {

			#----add the pct line to define TSS to be stored
			my $valueStatObj = Statistics::Descriptive::Full->new();
			$valueStatObj->add_data(@{$TSSmRNAEndFreqHsh_ref->{$ATGOrTAA}{$senseOrAntisense}});
	
			my $lowerPctLimit = $validTSSLimitHsh_ref->{$ATGOrTAA}{$senseOrAntisense}{'lowerPctLimit'};
			my $upperPctLimit = $validTSSLimitHsh_ref->{$ATGOrTAA}{$senseOrAntisense}{'upperPctLimit'};
	
			my $lowerValLimit = sprintf "%.3f", $valueStatObj->percentile($lowerPctLimit);
			my $upperValLimit = sprintf "%.3f", $valueStatObj->percentile($upperPctLimit);
			
			$validTSSLimitHsh_ref->{$ATGOrTAA}{$senseOrAntisense}{'lowerValLimit'} = $lowerValLimit;
			$validTSSLimitHsh_ref->{$ATGOrTAA}{$senseOrAntisense}{'upperValLimit'} = $upperValLimit;
			
			foreach my $mRNAID (keys %{$geneBasedTSSInfoHsh_ref}) {
				foreach my $i (0..$#{$geneBasedTSSInfoHsh_ref->{$mRNAID}{$ATGOrTAA}{$senseOrAntisense}{'rltvPos'}}) {
					my $rltvPos = ${$geneBasedTSSInfoHsh_ref->{$mRNAID}{$ATGOrTAA}{$senseOrAntisense}{'rltvPos'}}[$i];
					if ($rltvPos <= $upperValLimit and $rltvPos >= $lowerValLimit) {
						$strndEndValidTSSInfoHsh_ref->{$ATGOrTAA}{$senseOrAntisense}{$mRNAID}{$rltvPos}{'cov'} = ${$geneBasedTSSInfoHsh_ref->{$mRNAID}{$ATGOrTAA}{$senseOrAntisense}{'cov'}}[$i];
						$strndEndValidTSSInfoHsh_ref->{$ATGOrTAA}{$senseOrAntisense}{$mRNAID}{$rltvPos}{'absPos'} = ${$geneBasedTSSInfoHsh_ref->{$mRNAID}{$ATGOrTAA}{$senseOrAntisense}{'absPos'}}[$i];
					}
				}
			}
		}
	}
	
	return ($strndEndValidTSSInfoHsh_ref);
	
}
sub getIndivCntgCovPlsPath {
#....................................................................................................................................................#
#	subroutineCategory: storable
#	dependOnSub: reportStatus|3137
#	appearInSub: getmRNAReferencePoints|1943
#	primaryAppearInSection: 5_processInputData|174
#	secondaryAppearInSection: 7_retrieveSequenceSurroundingPredefinedTSS|212
#	input: $cntgCovPlsIndexPath
#	output: $cntgCovInPlsPathHsh_ref
#	toCall: my ($cntgCovInPlsPathHsh_ref) = &getIndivCntgCovPlsPath($cntgCovPlsIndexPath);
#	calledInLine: 194, 1968
#....................................................................................................................................................#
	
	my ($cntgCovPlsIndexPath) = @_;

	$cntgCovPlsIndexPath =~ s/\.gz$//;

	my $cntgCovInPlsPathHsh_ref = {};
	
	system ("gzip -fd $cntgCovPlsIndexPath.gz") if -s "$cntgCovPlsIndexPath.gz";
	my %plsIndexHsh = %{retrieve($cntgCovPlsIndexPath)};
	
	my (undef, $cntgCovStroableDir, undef) = fileparse($cntgCovPlsIndexPath, qr/\.[^.]*/);
	foreach my $cntg (keys %plsIndexHsh) {
		my $cntgCovPlsPath = "$cntgCovStroableDir/$plsIndexHsh{$cntg}";
		die "cntgCovPlsPath $cntgCovPlsPath is invalid\n" if ((not -s $cntgCovPlsPath) and (not -s $cntgCovPlsPath.".gz"));
		$cntgCovInPlsPathHsh_ref->{$cntg} = $cntgCovPlsPath;
	}
	my $numCntg = keys %{$cntgCovInPlsPathHsh_ref};
	&reportStatus("pls path of $numCntg contig stored.", 0, "\n");#->3137
	
	return $cntgCovInPlsPathHsh_ref;
}
sub getInrMotif {
#....................................................................................................................................................#
#	subroutineCategory: unassigned
#	dependOnSub: createWeblogo|667, reportStatus|3137
#	appearInSub: >none
#	primaryAppearInSection: 8_getInitiatorMotif|228
#	secondaryAppearInSection: >none
#	input: $seqAroundmRNATSSHsh_ref, $weblogoDirHsh_ref
#	output: $InrMotifFreqHsh_ref, $InrSeqAlignHsh_ref
#	toCall: my ($InrMotifFreqHsh_ref, $InrSeqAlignHsh_ref) = &getInrMotif($seqAroundmRNATSSHsh_ref, $weblogoDirHsh_ref);
#	calledInLine: 235
#....................................................................................................................................................#
	my ($seqAroundmRNATSSHsh_ref, $weblogoDirHsh_ref) = @_;
	my $startRltvPos = -4;
	my $endRltvPos = 2;
	
	my $InrMotifFreqHsh_ref = {};
	my $InrSeqAlignHsh_ref = {};
	
	&reportStatus("Analyzing Inr element sequence", 0, "\n");#->3137
	my $numGene = 0;
	foreach my $mRNAID (keys %{$seqAroundmRNATSSHsh_ref->{"mRNA_TSS"}}) {
		$numGene++;
		my $InrSeq = (substr $seqAroundmRNATSSHsh_ref->{"mRNA_TSS"}{$mRNAID}{'upStrm'}, $startRltvPos).(substr $seqAroundmRNATSSHsh_ref->{"mRNA_TSS"}{$mRNAID}{'dnStrm'}, 0, $endRltvPos);
		$InrSeqAlignHsh_ref->{$mRNAID} = $InrSeq;
		$InrMotifFreqHsh_ref->{'count'}{$InrSeq}++;
	}
	
	{#---plot weblogo
		my $nameTag = "InrMotif";
		my $seqAlignHsh_ref = $InrSeqAlignHsh_ref;
		my $pdfPath = "$weblogoDirHsh_ref->{pdf}/$nameTag.pdf";
		my $fastaPath = "$weblogoDirHsh_ref->{fasta}/$nameTag.fasta";
		my $cmdPath = "$weblogoDirHsh_ref->{cmd}/$nameTag.cmd";
		my $seqType = 'dna';
		my $title = "InrMotif N=$numGene";
		&createWeblogo($seqAlignHsh_ref, $pdfPath, $fastaPath, $cmdPath, $seqType, $title);#->667
	}
	
	my $totalCount = sum values %{$InrMotifFreqHsh_ref->{'count'}};
	my $cmltveProportion = 0;
	foreach my $InrSeq (sort {$InrMotifFreqHsh_ref->{'count'}{$b} <=> $InrMotifFreqHsh_ref->{'count'}{$a}} keys %{$InrMotifFreqHsh_ref->{'count'}}) {
		$InrMotifFreqHsh_ref->{'proportion'}{$InrSeq} = sprintf "%.10f", $InrMotifFreqHsh_ref->{'count'}{$InrSeq}/$totalCount;
		$cmltveProportion += $InrMotifFreqHsh_ref->{'proportion'}{$InrSeq};
		$InrMotifFreqHsh_ref->{'cmltveProportion'}{$InrSeq} = $cmltveProportion;
		
		#print TMPLOG $InrSeq."\t".$InrMotifFreqHsh_ref->{'proportion'}{$InrSeq}."\t".$InrMotifFreqHsh_ref->{'cmltveProportion'}{$InrSeq}."\n";
	}
	
	return ($InrMotifFreqHsh_ref, $InrSeqAlignHsh_ref);
}
sub getMastGenomeBothStrandHit {
#....................................................................................................................................................#
#	subroutineCategory: unassigned
#	dependOnSub: >none
#	appearInSub: scanMotifWholeGenomeWithMAST|3277
#	primaryAppearInSection: >none
#	secondaryAppearInSection: 9_scanMotifOccurence|246
#	input: $fastaLengthHsh_ref, $mastHitLog
#	output: $allHitBySeqHsh_ref
#	toCall: my ($allHitBySeqHsh_ref) = &getMastGenomeBothStrandHit($mastHitLog, $fastaLengthHsh_ref);
#	calledInLine: 3312
#....................................................................................................................................................#

	my ($mastHitLog, $fastaLengthHsh_ref) = @_;
	
	my %hitCountHsh = ();
	my $allHitBySeqHsh_ref = {};
	
	open MASTLOG, "<", $mastHitLog;
	while (<MASTLOG>) {
		next if $_ =~ m/^#/;
		chomp;
		my ($sequence_name, $motif, $hit_start, $hit_end, $score, $hit_pValue) = split / +/;
		my $strnd = $1 if $sequence_name =~ m/([\+\-])$/;
		$sequence_name =~ s/[\+\-]$//;
		#print TMPLOG "$sequence_name\t$strnd\t$hit_start\t$hit_pValue\n";
		if ($strnd eq '-') {
			$hit_start = $fastaLengthHsh_ref->{$sequence_name} - $hit_start + 1;
		}
		$allHitBySeqHsh_ref->{$sequence_name}{$strnd}{$hit_start} = $hit_pValue;
	}
	close MASTLOG;
	
	system ("pigz $mastHitLog");
	
	return ($allHitBySeqHsh_ref);
}
sub getMotifBoundCutoff {
#....................................................................................................................................................#
#	subroutineCategory: unassigned
#	dependOnSub: >none
#	appearInSub: defineMotifPositionBoundsUsingShuffleBackground|926
#	primaryAppearInSection: >none
#	secondaryAppearInSection: 10_defineMotifBoundary|259
#	input: $abvBkgrdFactor, $extraMargin, $minStreakHsh_ref, $motif, $posPctHsh_ref, $siteSeqRng
#	output: $leftBound, $postionFactorHsh_ref, $rightBound
#	toCall: my ($leftBound, $rightBound, $postionFactorHsh_ref) = &getMotifBoundCutoff($posPctHsh_ref, $minStreakHsh_ref, $abvBkgrdFactor, $extraMargin, $siteSeqRng, $motif);
#	calledInLine: 951
#....................................................................................................................................................#
	my ($posPctHsh_ref, $minStreakHsh_ref, $abvBkgrdFactor, $extraMargin, $siteSeqRng, $motif) = @_;
	
	my %regionPosHsh = ();
	my $regionNum = 0;
	my $postvStreak = 0;
	my $negtvStreak = 0;
	my $withinPostvStreak = 'no';
	my %regionSizeHsh = ();
	
	my $queryPeakPct = 0;
	my $queryPeakPos;
	
	foreach my $pos (sort {$a <=> $b} keys %{$posPctHsh_ref}) {
	
		if ($posPctHsh_ref->{$pos}{'query'} > $queryPeakPct) {
			$queryPeakPos = $pos;
			$queryPeakPct = $posPctHsh_ref->{$pos}{'query'};
		}
		
		if ($posPctHsh_ref->{$pos}{'query'} > $abvBkgrdFactor*$posPctHsh_ref->{$pos}{'shuffle'}) {
			$postvStreak++;
			$negtvStreak = 0;
			$regionNum++ if $postvStreak == $minStreakHsh_ref->{'positive'} and $withinPostvStreak eq "no";
			$withinPostvStreak = 'yes' if ($postvStreak >= $minStreakHsh_ref->{'positive'});
		} else {
			$negtvStreak++;
			$postvStreak = 0;
			$withinPostvStreak = 'no' if ($negtvStreak >= $minStreakHsh_ref->{'negative'});
		}

		if ($withinPostvStreak eq 'yes') {
			push @{$regionPosHsh{$regionNum}}, $pos;
			$regionSizeHsh{$regionNum}++;
		}
	}
	
	my $leftBound = my $rightBound = 0;
	my $postionFactorHsh_ref = {};
	
	if (exists $regionPosHsh{1} and $motif ne 'Inr') {
		my @regionAry = (sort {$regionSizeHsh{$b} <=> $regionSizeHsh{$a}} keys %regionSizeHsh);
		my $largestRegion = $regionAry[0];
		$leftBound = ${$regionPosHsh{$largestRegion}}[0] - $minStreakHsh_ref->{'positive'} - $extraMargin;
		$rightBound = ${$regionPosHsh{$largestRegion}}[-1] - $minStreakHsh_ref->{'negative'} + $extraMargin;
		my $peakValue = 0;
		for my $pos ($leftBound-1..$rightBound+1) {
			my $pct = $posPctHsh_ref->{$pos}{'query'};
			my $posToTSS = $pos - $siteSeqRng - 1;
			$pct = 0 if $pos < $leftBound or $pos > $rightBound; #---add end zero for better plotting
			$postionFactorHsh_ref->{$posToTSS} = $pct;
			$peakValue = $pct if $pct > $peakValue;
		}
		#---scale positional factor using peak pct value
		$postionFactorHsh_ref->{$_} = $postionFactorHsh_ref->{$_}/$peakValue foreach (keys %{$postionFactorHsh_ref});
	}

	if ($motif eq 'Inr') {
		$leftBound = $rightBound = $queryPeakPos;
		my $posToTSS = $queryPeakPos - $siteSeqRng - 1;
		$postionFactorHsh_ref->{$posToTSS} = 1; 
		$postionFactorHsh_ref->{$posToTSS+1} = 0; #---add end zero for better plotting
		$postionFactorHsh_ref->{$posToTSS-1} = 0; #---add end zero for better plotting
	}

	return ($leftBound, $rightBound, $postionFactorHsh_ref);
}
sub getProportionOfGenesWithTSSScoreAboveCutOff {
#....................................................................................................................................................#
#	subroutineCategory: unassigned
#	dependOnSub: reportStatus|3137
#	appearInSub: >none
#	primaryAppearInSection: 13_analyzeTSSScore|293
#	secondaryAppearInSection: >none
#	input: $TSSScoreIndivmRNAHsh_ref, $mRNANATInfoHsh_ref, $mRNARefPtHsh_ref
#	output: $geneAboveTSSScoreHsh_ref
#	toCall: my ($geneAboveTSSScoreHsh_ref) = &getProportionOfGenesWithTSSScoreAboveCutOff($TSSScoreIndivmRNAHsh_ref, $mRNANATInfoHsh_ref, $mRNARefPtHsh_ref);
#	calledInLine: 310
#....................................................................................................................................................#
	my ($TSSScoreIndivmRNAHsh_ref, $mRNANATInfoHsh_ref, $mRNARefPtHsh_ref) = @_;
	
	my $geneAboveTSSScoreHsh_ref = {};
	
	my @cutoffAry = (0..10);
	
	foreach my $cutoff (@cutoffAry) {
		foreach my $refPosType (keys %{$TSSScoreIndivmRNAHsh_ref}) {

			&reportStatus("Getting proportion of genes with TSSScore > $cutoff in $refPosType", 20, "\r");#->3137
		
			foreach my $dirtn (qw/a s/) {
				my %geneAboveCutoffHsh = ();
				my %geneNumHsh = ();
				foreach my $mRNAID (keys %{$TSSScoreIndivmRNAHsh_ref->{$refPosType}}) {
				
					my $end3NAT = $mRNANATInfoHsh_ref->{$mRNAID}{'end3NAT'};
					$geneNumHsh{$end3NAT}++;
					$geneNumHsh{'withNAT_TSS'}++ if exists $mRNARefPtHsh_ref->{'NAT_TSS'}{$mRNAID};
				
					foreach my $rltvPos (keys %{$TSSScoreIndivmRNAHsh_ref->{$refPosType}{$mRNAID}}) {
						if (not exists $mRNARefPtHsh_ref->{'NAT_TSS'}{$mRNAID}) {
							$geneAboveCutoffHsh{$end3NAT}{$mRNAID}++ if ($TSSScoreIndivmRNAHsh_ref->{$refPosType}{$mRNAID}{$rltvPos}{$dirtn} > $cutoff);
						} else {
							$geneAboveCutoffHsh{'withNAT_TSS'}{$mRNAID}++ if ($TSSScoreIndivmRNAHsh_ref->{$refPosType}{$mRNAID}{$rltvPos}{$dirtn} > $cutoff);
						}
					}
				}
			
				foreach my $category (keys %geneNumHsh) {
					my $geneNum = $geneNumHsh{$category};
					my $numAbvCutoff = 0;
					$numAbvCutoff = keys %{$geneAboveCutoffHsh{$category}} if exists $geneAboveCutoffHsh{$category};
					my $prprtnAbvCutoff = $numAbvCutoff/$geneNum;
					$geneAboveTSSScoreHsh_ref->{$refPosType}{$category}{$dirtn}{$cutoff}{'geneNum'} = $geneNum;
					$geneAboveTSSScoreHsh_ref->{$refPosType}{$category}{$dirtn}{$cutoff}{'numAbvCutoff'} = $numAbvCutoff;
					$geneAboveTSSScoreHsh_ref->{$refPosType}{$category}{$dirtn}{$cutoff}{'prprtnAbvCutoff'} = $prprtnAbvCutoff;
					#print TMPLOG join "", ((join "\t", ($refPosType, $category, $dirtn, $cutoff, $geneNum, $numAbvCutoff, $prprtnAbvCutoff)), "\n");
				}
			}
		}
	}
	
	return ($geneAboveTSSScoreHsh_ref);

	return ();
}
sub getSequenceAroundmRNAReferencePoint {
#....................................................................................................................................................#
#	subroutineCategory: unassigned
#	dependOnSub: reportStatus|3137
#	appearInSub: >none
#	primaryAppearInSection: 7_retrieveSequenceSurroundingPredefinedTSS|212
#	secondaryAppearInSection: >none
#	input: $fastaHsh_ref, $mRNAInfoHsh_ref, $mRNARefPtHsh_ref, $resultFastaDir, $resultStorableDir, $siteSeqRng
#	output: $seqAroundSiteHsh_ref, $seqAroundSiteInfoHsh_ref
#	toCall: my ($seqAroundSiteHsh_ref, $seqAroundSiteInfoHsh_ref) = &getSequenceAroundmRNAReferencePoint($mRNAInfoHsh_ref, $mRNARefPtHsh_ref, $fastaHsh_ref, $siteSeqRng, $resultFastaDir, $resultStorableDir);
#	calledInLine: 219
#....................................................................................................................................................#

	my ($mRNAInfoHsh_ref, $mRNARefPtHsh_ref, $fastaHsh_ref, $siteSeqRng, $resultFastaDir, $resultStorableDir) = @_;
	
	my $seqAroundSiteHsh_ref = {};
	my $seqAroundSiteInfoHsh_ref = {};
	
	my $seqAroundSiteHshPlsPath = "$resultStorableDir/seqAroundSiteHsh.pls";
	my $seqAroundSiteInfoHshPlsPath = "$resultStorableDir/seqAroundSiteInfoHsh.pls";
	
	foreach my $siteType (keys %{$mRNARefPtHsh_ref}) {
		my $fastaPath = "$resultFastaDir/$siteType.full.fasta";
		$seqAroundSiteInfoHsh_ref->{$siteType}{'fastaPath'} = $fastaPath;
		$seqAroundSiteInfoHsh_ref->{$siteType}{'length'} = $siteSeqRng*2;
		$seqAroundSiteInfoHsh_ref->{$siteType}{'totalSeqNum'} = keys %{$mRNARefPtHsh_ref->{$siteType}};
	
		open (FASTA, ">", $fastaPath);
		&reportStatus("Getting sequences around $siteType", 0, "\n");#->3137
		foreach my $mRNAID (keys %{$mRNARefPtHsh_ref->{$siteType}}) {
			if ($mRNAInfoHsh_ref->{$mRNAID}{'strnd'} eq '+') {

				$seqAroundSiteHsh_ref->{$siteType}{$mRNAID}{'upStrm'} = substr $fastaHsh_ref->{$mRNAInfoHsh_ref->{$mRNAID}{'cntg'}}, $mRNARefPtHsh_ref->{$siteType}{$mRNAID}-$siteSeqRng-1, $siteSeqRng;
				$seqAroundSiteHsh_ref->{$siteType}{$mRNAID}{'dnStrm'} = substr $fastaHsh_ref->{$mRNAInfoHsh_ref->{$mRNAID}{'cntg'}}, $mRNARefPtHsh_ref->{$siteType}{$mRNAID}-1, $siteSeqRng;

			} else {
			
				$seqAroundSiteHsh_ref->{$siteType}{$mRNAID}{'upStrm'} = substr $fastaHsh_ref->{$mRNAInfoHsh_ref->{$mRNAID}{'cntg'}}, $mRNARefPtHsh_ref->{$siteType}{$mRNAID}-1+1, $siteSeqRng;
				$seqAroundSiteHsh_ref->{$siteType}{$mRNAID}{'dnStrm'} = substr $fastaHsh_ref->{$mRNAInfoHsh_ref->{$mRNAID}{'cntg'}}, $mRNARefPtHsh_ref->{$siteType}{$mRNAID}-$siteSeqRng-1+1, $siteSeqRng;
			
				foreach my $upStrmOrDnStrm (keys %{$seqAroundSiteHsh_ref->{$siteType}{$mRNAID}}) {
					$seqAroundSiteHsh_ref->{$siteType}{$mRNAID}{$upStrmOrDnStrm} = reverse $seqAroundSiteHsh_ref->{$siteType}{$mRNAID}{$upStrmOrDnStrm};
					$seqAroundSiteHsh_ref->{$siteType}{$mRNAID}{$upStrmOrDnStrm} =~ tr/ACGTacgt/TGCAtgca/;
				}
			}
		
			#---reverse complement the seq if the purpose is to investigate NAT
			if ($siteType eq 'mRNA_TAA' or $siteType eq 'NAT_TSS') {
				($seqAroundSiteHsh_ref->{$siteType}{$mRNAID}{'upStrm'}, $seqAroundSiteHsh_ref->{$siteType}{$mRNAID}{'dnStrm'}) = ($seqAroundSiteHsh_ref->{$siteType}{$mRNAID}{'dnStrm'}, $seqAroundSiteHsh_ref->{$siteType}{$mRNAID}{'upStrm'});
				foreach my $upStrmOrDnStrm (keys %{$seqAroundSiteHsh_ref->{$siteType}{$mRNAID}}) {
					$seqAroundSiteHsh_ref->{$siteType}{$mRNAID}{$upStrmOrDnStrm} = reverse $seqAroundSiteHsh_ref->{$siteType}{$mRNAID}{$upStrmOrDnStrm};
					$seqAroundSiteHsh_ref->{$siteType}{$mRNAID}{$upStrmOrDnStrm} =~ tr/ACGTacgt/TGCAtgca/;
				}
			}
		
			#print TMPLOG join "\t", ($siteType, $mRNAInfoHsh_ref->{$mRNAID}{'strnd'}, $mRNARefPtHsh_ref->{$siteType}{$mRNAID}, $mRNAID, $seqAroundSiteHsh_ref->{$siteType}{$mRNAID}{'upStrm'}, $seqAroundSiteHsh_ref->{$siteType}{$mRNAID}{'dnStrm'}."\n");
		
			print FASTA ">$mRNAID\n";
			print FASTA $seqAroundSiteHsh_ref->{$siteType}{$mRNAID}{'upStrm'}.$seqAroundSiteHsh_ref->{$siteType}{$mRNAID}{'dnStrm'}."\n";
		}
		close FASTA
	}
	
	store($seqAroundSiteHsh_ref, $seqAroundSiteHshPlsPath);
	store($seqAroundSiteInfoHsh_ref, $seqAroundSiteInfoHshPlsPath);
	
	return ($seqAroundSiteHsh_ref, $seqAroundSiteInfoHsh_ref);
}
sub getSingleMotifMASTLogPostionalData {
#....................................................................................................................................................#
#	subroutineCategory: unassigned
#	dependOnSub: >none
#	appearInSub: scanMotifAroundSiteWithMAST|3210
#	primaryAppearInSection: >none
#	secondaryAppearInSection: 9_scanMotifOccurence|246
#	input: $mastHitLog, $maxHitPVal, $maxPos, $totalSeqNum
#	output: $allHitBySeqHsh_ref, $motifPctHsh_ref
#	toCall: my ($motifPctHsh_ref, $allHitBySeqHsh_ref) = &getSingleMotifMASTLogPostionalData($mastHitLog, $maxHitPVal, $maxPos, $totalSeqNum);
#	calledInLine: 3253, 3264
#....................................................................................................................................................#
	my ($mastHitLog, $maxHitPVal, $maxPos, $totalSeqNum) = @_;
	
	my $motifPctHsh_ref = {};
	my %hitCountHsh = ();
	my $allHitBySeqHsh_ref = {};
	
	open MASTLOG, "<", $mastHitLog;
	while (<MASTLOG>) {
		next if $_ =~ m/^#/;
		chomp;
		my ($sequence_name, $motif, $hit_start, $hit_end, $score, $hit_pValue) = split / +/;
		$hitCountHsh{$hit_start}++ if ($hit_pValue <= $maxHitPVal);
		$allHitBySeqHsh_ref->{$sequence_name}{$hit_start} = $hit_pValue;
	}
	close MASTLOG;
	
	foreach my $pos (1..$maxPos) {
		$motifPctHsh_ref->{$pos} = 0;
		$motifPctHsh_ref->{$pos} = 100*$hitCountHsh{$pos}/$totalSeqNum if $hitCountHsh{$pos};
	}

	return ($motifPctHsh_ref, $allHitBySeqHsh_ref);
}
sub getTSSScoreAroundmRNAReferencePoint {
#....................................................................................................................................................#
#	subroutineCategory: unassigned
#	dependOnSub: reportStatus|3137
#	appearInSub: >none
#	primaryAppearInSection: 13_analyzeTSSScore|293
#	secondaryAppearInSection: >none
#	input: $genomeWideTSSPlsPathHsh_ref, $mRNAByCntgHsh_ref, $mRNAInfoHsh_ref, $mRNARefPtHsh_ref, $siteSeqRng
#	output: $TSSScoreRefPtBymRNAHsh_ref, $TSSScoreRefPtPlotHsh_ref
#	toCall: my ($TSSScoreRefPtPlotHsh_ref, $TSSScoreRefPtBymRNAHsh_ref) = &getTSSScoreAroundmRNAReferencePoint($mRNAInfoHsh_ref, $mRNARefPtHsh_ref, $genomeWideTSSPlsPathHsh_ref, $siteSeqRng, $mRNAByCntgHsh_ref);
#	calledInLine: 301
#....................................................................................................................................................#
	my ($mRNAInfoHsh_ref, $mRNARefPtHsh_ref, $genomeWideTSSPlsPathHsh_ref, $siteSeqRng, $mRNAByCntgHsh_ref) = @_;
	
	my $TSSScoreRefPtPlotHsh_ref = {};
	my $TSSScoreRefPtBymRNAHsh_ref = {};
	
	my %tmpStrndDirtnHsh = ();
	$tmpStrndDirtnHsh{'+'}{'s'} = '+';
	$tmpStrndDirtnHsh{'+'}{'a'} = '-';
	$tmpStrndDirtnHsh{'-'}{'s'} = '-';
	$tmpStrndDirtnHsh{'-'}{'a'} = '+';
	
	foreach my $sumOrCount (qw/sum count/) {
		foreach my $siteType (keys %{$mRNARefPtHsh_ref}) {
			foreach my $mode (qw/allPos mostUpStrmPos/) {
				foreach my $rltvPos (-1*$siteSeqRng..$siteSeqRng-1) {
					foreach my $dirtn (qw/a s/) {
						$TSSScoreRefPtPlotHsh_ref->{$siteType}{$mode}{$rltvPos}{$dirtn}{$sumOrCount} = 0;
					}
				}
			}
		}
	}

	foreach my $cntg (keys %{$mRNAByCntgHsh_ref}) {
		&reportStatus("Getting TSS scores around site in $cntg", 20, "\r");#->3137

		my $cntgTSSScoreAry_ref = retrieve($genomeWideTSSPlsPathHsh_ref->{$cntg});
		
		foreach my $mRNAID (keys %{$mRNAByCntgHsh_ref->{$cntg}}) {
			foreach my $siteType (keys %{$mRNARefPtHsh_ref}) {
				my @upStrmPosRng = my @dnStrmPosRng = ();
				if (exists $mRNARefPtHsh_ref->{$siteType}{$mRNAID}) {
					my $geneStrnd = $mRNAInfoHsh_ref->{$mRNAID}{'strnd'};
					my $refPos = $mRNARefPtHsh_ref->{$siteType}{$mRNAID};
					if ($geneStrnd eq '+') {
						@upStrmPosRng = ($refPos-$siteSeqRng-1..$refPos-1);
						@dnStrmPosRng = ($refPos..$refPos+$siteSeqRng);
					} else {
						@upStrmPosRng = reverse ($refPos+1..$refPos+$siteSeqRng+1);
						@dnStrmPosRng = reverse ($refPos-$siteSeqRng..$refPos);
					}
					
					my @searchRngAry = (@upStrmPosRng, @dnStrmPosRng);
					
					my $rltvPos = -1*$siteSeqRng;
					foreach my $pos (@searchRngAry) {
						my $i = $pos - 1;
						if ($cntgTSSScoreAry_ref->[$i]) {
							my %tmpScoreHsh = ();
							($tmpScoreHsh{'+'}, $tmpScoreHsh{'-'}) = split ",", $cntgTSSScoreAry_ref->[$i];
							foreach my $dirtn (qw/a s/) {
								my $dirtnStrnd = $tmpStrndDirtnHsh{$geneStrnd}{$dirtn};
								$TSSScoreRefPtBymRNAHsh_ref->{$mRNAID}{$siteType}{$dirtn}{$rltvPos} = $tmpScoreHsh{$dirtnStrnd} if $tmpScoreHsh{$dirtnStrnd} > 0;
							}
						}
						$rltvPos++;
					}
					foreach my $dirtn (qw/a s/) {
						if (exists $TSSScoreRefPtBymRNAHsh_ref->{$mRNAID}{$siteType}{$dirtn}) {
							#---[12/10/2013 13:46] get sorted rltvPos for the $TSSScoreRefPtMode eq 'upStrmMost'
							my @rltvPosAry = sort {$a <=> $b} (keys %{$TSSScoreRefPtBymRNAHsh_ref->{$mRNAID}{$siteType}{$dirtn}});
							@rltvPosAry = reverse @rltvPosAry if $dirtn eq 'a'; #---[12/10/2013 13:46] reverse for antisense direction, max rltv pos is the upstream
							
							#---[12/10/2013 13:58] go thr the rltv pos from up to down stream
							foreach my $rltvPos (@rltvPosAry) {
								$TSSScoreRefPtPlotHsh_ref->{$siteType}{'allPos'}{$rltvPos}{$dirtn}{'sum'} += $TSSScoreRefPtBymRNAHsh_ref->{$mRNAID}{$siteType}{$dirtn}{$rltvPos};
								$TSSScoreRefPtPlotHsh_ref->{$siteType}{'allPos'}{$rltvPos}{$dirtn}{'count'}++;
							}
							
							#---[12/10/2013 14:05] get only the most upstream position
							my $mostUpStrmPos = $rltvPosAry[0];
							$TSSScoreRefPtPlotHsh_ref->{$siteType}{'mostUpStrmPos'}{$mostUpStrmPos}{$dirtn}{'sum'} += $TSSScoreRefPtBymRNAHsh_ref->{$mRNAID}{$siteType}{$dirtn}{$mostUpStrmPos};
							$TSSScoreRefPtPlotHsh_ref->{$siteType}{'mostUpStrmPos'}{$mostUpStrmPos}{$dirtn}{'count'}++;
						}
					}
				}
			}
		}
	}
	
	return ($TSSScoreRefPtPlotHsh_ref, $TSSScoreRefPtBymRNAHsh_ref);
}
sub getTSSScoreForAllGenes {
#....................................................................................................................................................#
#	subroutineCategory: unassigned
#	dependOnSub: reportStatus|3137
#	appearInSub: >none
#	primaryAppearInSection: 13_analyzeTSSScore|293
#	secondaryAppearInSection: >none
#	input: $TSSScoreStorableDir, $TSScoreRegSizeHsh_ref, $fastaLengthHsh_ref, $forceCalculateTSSScoreIndivmRNA, $genomeWideTSSPlsPathHsh_ref, $mRNAByCntgHsh_ref, $mRNAInfoHsh_ref, $skipEndNt
#	output: $TSSScoreIndivmRNAHsh_ref
#	toCall: my ($TSSScoreIndivmRNAHsh_ref) = &getTSSScoreForAllGenes($mRNAInfoHsh_ref, $genomeWideTSSPlsPathHsh_ref, $mRNAByCntgHsh_ref, $skipEndNt, $TSScoreRegSizeHsh_ref, $fastaLengthHsh_ref, $TSSScoreStorableDir, $forceCalculateTSSScoreIndivmRNA);
#	calledInLine: 304
#....................................................................................................................................................#
	my ($mRNAInfoHsh_ref, $genomeWideTSSPlsPathHsh_ref, $mRNAByCntgHsh_ref, $skipEndNt, $TSScoreRegSizeHsh_ref, $fastaLengthHsh_ref, $TSSScoreStorableDir, $forceCalculateTSSScoreIndivmRNA) = @_;
	
	my $TSSScoreIndivmRNAHsh_ref = {};
	
	my $TSSScoreIndivmRNAHshPlsPath = "$TSSScoreStorableDir/TSSScoreIndivmRNAHsh.pls";
	
	my $upStrmSize = $TSScoreRegSizeHsh_ref->{'upStrm'};
	my $dnStrmSize = $TSScoreRegSizeHsh_ref->{'dnStrm'};
	my $safetyBoundSize = $upStrmSize;
	$safetyBoundSize = $dnStrmSize if $dnStrmSize > $upStrmSize;
	
	if (-s $TSSScoreIndivmRNAHshPlsPath and $forceCalculateTSSScoreIndivmRNA eq 'no') {
	
		&reportStatus("Retrieving TSSScoreIndivmRNAHshPlsPath", 0, "\n");#->3137
		$TSSScoreIndivmRNAHsh_ref = retrieve($TSSScoreIndivmRNAHshPlsPath);
	
	} else {
		
		my %tmpStrndDirtnHsh = ();
		$tmpStrndDirtnHsh{'+'}{'s'} = '+';
		$tmpStrndDirtnHsh{'+'}{'a'} = '-';
		$tmpStrndDirtnHsh{'-'}{'s'} = '-';
		$tmpStrndDirtnHsh{'-'}{'a'} = '+';

		foreach my $cntg (keys %{$mRNAByCntgHsh_ref}) {
		
			&reportStatus("Getting TSS scores all mRNAs in $cntg", 20, "\r");#->3137
		
			my $cntgTSSScoreAry_ref = retrieve($genomeWideTSSPlsPathHsh_ref->{$cntg});

			foreach my $mRNAID (keys %{$mRNAByCntgHsh_ref->{$cntg}}) {

				my $geneStrnd = $mRNAInfoHsh_ref->{$mRNAID}{'strnd'};

				#----get the gene ends as reference points
				my @geneRngAry = sort {$a <=> $b} @{$mRNAInfoHsh_ref->{$mRNAID}{'geneRng'}};
				my ($geneStart, $geneEnd) = ($geneRngAry[0], $geneRngAry[-1]);

				#---skip the gene ta the edge and the small genes
				next if $geneStart+2*$safetyBoundSize > $geneEnd-2*$safetyBoundSize or $geneStart - $safetyBoundSize <= 0 or $geneEnd + $safetyBoundSize > $fastaLengthHsh_ref->{$cntg};

				#---get a random within gene position
				my @withinGeneRngAry = shuffle($geneStart+2*$safetyBoundSize..$geneEnd-2*$safetyBoundSize);
				my $randWithinGenePos = $withinGeneRngAry[0];

				my ($ATGPos, $TAAPos) = ($geneStart, $geneEnd);
				($ATGPos, $TAAPos) = ($TAAPos, $ATGPos) if $geneStrnd eq '-';
			
				my %refPosHsh = ('ATG'=>$ATGPos, 'TAA'=>$TAAPos, 'randWithin'=>$randWithinGenePos);
			
				foreach my $refPosType (keys %refPosHsh) {
					my $refPos = $refPosHsh{$refPosType};
					my @upStrmPosRng = my @dnStrmPosRng = ();
					if ($geneStrnd eq '+') {
						@upStrmPosRng = ($refPos-$upStrmSize-1-$skipEndNt..$refPos-1-$skipEndNt);
						@dnStrmPosRng = ($refPos+$skipEndNt..$refPos+$dnStrmSize+$skipEndNt);
					} else {
						@upStrmPosRng = reverse ($refPos+1+$skipEndNt..$refPos+$upStrmSize+1+$skipEndNt);
						@dnStrmPosRng = reverse ($refPos-$dnStrmSize-$skipEndNt..$refPos-$skipEndNt);
					}
					
					my @searchRngAry = (@upStrmPosRng, @dnStrmPosRng);
					
					my $rltvPos = -1*$upStrmSize;
					foreach my $pos (@searchRngAry) {
						my $i = $pos - 1;
						if ($cntgTSSScoreAry_ref->[$i]) {
							my %tmpScoreHsh = ();
							($tmpScoreHsh{'+'}, $tmpScoreHsh{'-'}) = split ",", $cntgTSSScoreAry_ref->[$i];
							foreach my $dirtn (keys %{$tmpStrndDirtnHsh{$geneStrnd}}) {
								$TSSScoreIndivmRNAHsh_ref->{$refPosType}{$mRNAID}{$rltvPos}{$dirtn} = $tmpScoreHsh{$tmpStrndDirtnHsh{$geneStrnd}{$dirtn}};
							}
						} else {
							foreach my $dirtn (keys %{$tmpStrndDirtnHsh{$geneStrnd}}) {
								$TSSScoreIndivmRNAHsh_ref->{$refPosType}{$mRNAID}{$rltvPos}{$dirtn} = 0;
							}
						}
						$rltvPos++;
					}
				}
			}
		}
		
		store($TSSScoreIndivmRNAHsh_ref, $TSSScoreIndivmRNAHshPlsPath);
	}

	return ($TSSScoreIndivmRNAHsh_ref);
}
sub getTSSScoreInExonAndTSS {
#....................................................................................................................................................#
#	subroutineCategory: unassigned
#	dependOnSub: reportStatus|3137
#	appearInSub: >none
#	primaryAppearInSection: 13_analyzeTSSScore|293
#	secondaryAppearInSection: >none
#	input: $genomeWideTSSPlsPathHsh_ref, $mRNAByCntgHsh_ref, $mRNAInfoHsh_ref, $mRNARefPtHsh_ref, $siteSeqRng
#	output: $TSSScoreExonTSSHshNonZero_ref, $TSSScoreExonTSSHshWithZero_ref
#	toCall: my ($TSSScoreExonTSSHshNonZero_ref, $TSSScoreExonTSSHshWithZero_ref) = &getTSSScoreInExonAndTSS($genomeWideTSSPlsPathHsh_ref, $mRNAByCntgHsh_ref, $mRNARefPtHsh_ref, $siteSeqRng, $mRNAInfoHsh_ref);
#	calledInLine: 298
#....................................................................................................................................................#
	my ($genomeWideTSSPlsPathHsh_ref, $mRNAByCntgHsh_ref, $mRNARefPtHsh_ref, $siteSeqRng, $mRNAInfoHsh_ref) = @_;
	
	
	my $TSSScoreExonTSSHshNonZero_ref = {};
	my $TSSScoreExonTSSHshWithZero_ref = {};
	
	foreach my $cntg (keys %{$mRNAByCntgHsh_ref}) {
		&reportStatus("Getting TSS scores at TSS and exon of genes on $cntg", 20, "\r");#->3137

		my $cntgTSSScoreAry_ref = retrieve($genomeWideTSSPlsPathHsh_ref->{$cntg});
		
		foreach my $mRNAID (keys %{$mRNAByCntgHsh_ref->{$cntg}}) {
			next if not exists $mRNARefPtHsh_ref->{'mRNA_TSS'}{$mRNAID};

			my @geneRngAry = sort {$a <=> $b} @{$mRNAInfoHsh_ref->{$mRNAID}{'geneRng'}};
			my ($geneStart, $geneEnd) = ($geneRngAry[0], $geneRngAry[-1]);
			next if $geneStart+$siteSeqRng > $geneEnd-$siteSeqRng;

			my %exonTSSPosHsh = ();
			@{$exonTSSPosHsh{'TSS'}} = ($mRNARefPtHsh_ref->{'mRNA_TSS'}{$mRNAID});
			
			@{$exonTSSPosHsh{'exon'}} = shuffle($geneStart+$siteSeqRng..$geneEnd-$siteSeqRng);

			my $geneStrnd = $mRNAInfoHsh_ref->{$mRNAID}{'strnd'};
			foreach my $exonOrTSS (keys %exonTSSPosHsh) {
				#----will sample non-zero numbers on exon and TSS, for histogram
				foreach my $pos (@{$exonTSSPosHsh{$exonOrTSS}}) {
					my $i = $pos-1;
					my %tmpScoreHsh = ();
					$tmpScoreHsh{$exonOrTSS}{$_} = 0 foreach (qw/+ -/);
					($tmpScoreHsh{$exonOrTSS}{'+'}, $tmpScoreHsh{$exonOrTSS}{'-'}) = split ",", $cntgTSSScoreAry_ref->[$i] if $cntgTSSScoreAry_ref->[$i];
					push @{$TSSScoreExonTSSHshNonZero_ref->{$exonOrTSS}}, $tmpScoreHsh{$exonOrTSS}{$geneStrnd} if $tmpScoreHsh{$exonOrTSS}{$geneStrnd} > 0;
				}

				#----will sample a single position on exon and TSS no matter zero or not for histogram
				my $singlePos = ${$exonTSSPosHsh{$exonOrTSS}}[0];
				my $i = $singlePos-1;
				my %tmpScoreHsh = ();
				$tmpScoreHsh{$exonOrTSS}{$_} = 0 foreach (qw/+ -/);
				($tmpScoreHsh{$exonOrTSS}{'+'}, $tmpScoreHsh{$exonOrTSS}{'-'}) = split ",", $cntgTSSScoreAry_ref->[$i] if $cntgTSSScoreAry_ref->[$i];
				push @{$TSSScoreExonTSSHshWithZero_ref->{$exonOrTSS}}, $tmpScoreHsh{$exonOrTSS}{$geneStrnd};
			}
		}
	}

	return ($TSSScoreExonTSSHshNonZero_ref, $TSSScoreExonTSSHshWithZero_ref);
}
sub getmRNAReferencePoints {
#....................................................................................................................................................#
#	subroutineCategory: specific
#	dependOnSub: getBaseAtTSSAndExon|1182, getGeneWithValidTSSAtGeneEnd|1283, getIndivCntgCovPlsPath|1331, investigateTSSRelativeToATGAndTAA|2169
#	appearInSub: >none
#	primaryAppearInSection: 7_retrieveSequenceSurroundingPredefinedTSS|212
#	secondaryAppearInSection: >none
#	input: $mRNAInfoHsh_ref, $validTSSCntgPlsPath
#	output: $mRNARefPtHsh_ref
#	toCall: my ($mRNARefPtHsh_ref) = &getmRNAReferencePoints($mRNAInfoHsh_ref, $validTSSCntgPlsPath);
#	calledInLine: 217
#....................................................................................................................................................#
	my ($mRNAInfoHsh_ref, $validTSSCntgPlsPath) = @_;
	
	my $TSSSearchRngHsh_ref = {};
	$TSSSearchRngHsh_ref->{'ATG'}{'up'} = 150;
	$TSSSearchRngHsh_ref->{'ATG'}{'dn'} = 150;
	$TSSSearchRngHsh_ref->{'TAA'}{'up'} = 150;
	$TSSSearchRngHsh_ref->{'TAA'}{'dn'} = 150;

	my $validTSSLimitHsh_ref = {};
	$validTSSLimitHsh_ref->{'ATG'}{'sense'}{'lowerPctLimit'} = 5;
	$validTSSLimitHsh_ref->{'ATG'}{'sense'}{'upperPctLimit'} = 95;
	$validTSSLimitHsh_ref->{'TAA'}{'antisense'}{'lowerPctLimit'} = 5;
	$validTSSLimitHsh_ref->{'TAA'}{'antisense'}{'upperPctLimit'} = 95;
	
	my $validTSSCovPlsPathHsh_ref = &getIndivCntgCovPlsPath($validTSSCntgPlsPath);#->1331
	my ($geneBasedTSSInfoHsh_ref, $TSSmRNAEndFreqHsh_ref) = &investigateTSSRelativeToATGAndTAA($validTSSCovPlsPathHsh_ref, $mRNAInfoHsh_ref, $TSSSearchRngHsh_ref);#->2169
	my ($strndEndValidTSSInfoHsh_ref) = &getGeneWithValidTSSAtGeneEnd($geneBasedTSSInfoHsh_ref, $TSSmRNAEndFreqHsh_ref, $validTSSLimitHsh_ref);#->1283
	my ($mRNARefPtHsh_ref) = &getBaseAtTSSAndExon($strndEndValidTSSInfoHsh_ref, $mRNAInfoHsh_ref);#->1182

	return ($mRNARefPtHsh_ref);
}
sub ggplotCulmulativeFrequency {
#....................................................................................................................................................#
#	subroutineCategory: unassigned
#	dependOnSub: >none
#	appearInSub: plotmRNANATInfo|2552
#	primaryAppearInSection: >none
#	secondaryAppearInSection: 6_defineGeneWithEnd3NAT|200
#	input: $RScriptPath, $dataPath, $extraArg, $height, $log2OrLinear, $logPath, $pdfPath, $plotAry_ref, $width, $xAxis
#	output: 
#	toCall: &ggplotCulmulativeFrequency($plotAry_ref, $dataPath, $pdfPath, $RScriptPath, $logPath, $xAxis, $extraArg, $log2OrLinear, $height, $width);
#	calledInLine: 2584
#....................................................................................................................................................#
	
	my ($plotAry_ref, $dataPath, $pdfPath, $RScriptPath, $logPath, $xAxis, $extraArg, $log2OrLinear, $height, $width) = @_;
	
	open (PLOTDATA, ">", $dataPath);
	print PLOTDATA $xAxis."\n";
	foreach my $value (@{$plotAry_ref}) {
		print PLOTDATA $value."\n";
	}
	close PLOTDATA;
	
	my $scale = '';
	$scale = ' + scale_x_continuous(trans=log2_trans())' if $log2OrLinear eq 'log2';
	
	my $numValue = @{$plotAry_ref};
	open (R, ">", $RScriptPath);
	print R "library(ggplot2)"."\n";
	print R "library(scales)"."\n";
	print R "dataFrame = read.table(file=\"$dataPath\", header=TRUE)"."\n";
	print R "ggplot(dataFrame, aes(x=$xAxis)) + ggtitle(\"Culmulative frequency of $xAxis $log2OrLinear scale [n=$numValue]\") + stat_ecdf() $scale $extraArg"."\n";
	print R "ggsave(file=\"$pdfPath\", height=$height, width=$width)\n";
	close R;
	
	system ("R --slave --vanilla --file=$RScriptPath &>$logPath ");

	return ();
	
}
sub ggplotMultiSampleBoxWhisker {
#....................................................................................................................................................#
#	subroutineCategory: unassigned
#	dependOnSub: >none
#	appearInSub: plotAverageTSSScoreForEnd3NATGenes|2269, plotTSSScoreInExonAndTSS|2435
#	primaryAppearInSection: >none
#	secondaryAppearInSection: 13_analyzeTSSScore|293
#	input: $RScriptPath, $dataPath, $dataPtMax, $extraArg, $height, $log2OrLinear, $logPath, $pdfPath, $plotAryHsh_ref, $width, $yAxis
#	output: 
#	toCall: &ggplotMultiSampleBoxWhisker($plotAryHsh_ref, $dataPath, $pdfPath, $RScriptPath, $logPath, $yAxis, $extraArg, $log2OrLinear, $dataPtMax, $height, $width);
#	calledInLine: 2307, 2487
#....................................................................................................................................................#
	my ($plotAryHsh_ref, $dataPath, $pdfPath, $RScriptPath, $logPath, $yAxis, $extraArg, $log2OrLinear, $dataPtMax, $height, $width) = @_;
	
	open (PLOTDATA, ">", $dataPath);
	print PLOTDATA "sample\t$yAxis\n";

	foreach my $sample (keys %{$plotAryHsh_ref}) {
		my $plotAry_ref = $plotAryHsh_ref->{$sample};

		my $indivDataPtMax = $dataPtMax;
		$indivDataPtMax = @{$plotAry_ref} if $indivDataPtMax > @{$plotAry_ref};

		#---down sample the data point number
		my @shuffleIndexAry = shuffle(0..$#{$plotAry_ref});
		foreach my $i (0..$indivDataPtMax-1) {
			my $value = $plotAry_ref->[$shuffleIndexAry[$i]];
			print PLOTDATA "$sample\t$value\n";
		}
		print PLOTDATA "\n";
	}
	
	my $scale = '';
	$scale = ' + scale_y_continuous(trans=log2_trans())' if $log2OrLinear eq 'log2';
	open (R, ">", $RScriptPath);
	print R "library(ggplot2)"."\n";
	print R "library(scales)"."\n";
	print R "dataFrame = read.table(file=\"$dataPath\", header=TRUE)"."\n";
	print R "ggplot(dataFrame, aes(x=sample, y=$yAxis)) + ggtitle(\"Boxplot of $yAxis in $log2OrLinear scale\") + geom_boxplot(aes(fill = sample)) $scale $extraArg"."\n";
	print R "ggsave(file=\"$pdfPath\", height=$height, width=$width)\n";
	
	system ("R --slave --vanilla --file=$RScriptPath &>$logPath ");

	return ();
}
sub ggplotTwoSampleHistogram {
#....................................................................................................................................................#
#	subroutineCategory: ggplot, plotInR
#	dependOnSub: >none
#	appearInSub: plotTSSScoreInExonAndTSS|2435
#	primaryAppearInSection: >none
#	secondaryAppearInSection: 13_analyzeTSSScore|293
#	input: $RScriptPath, $binWidth, $dataPath, $dataPtMax, $densityOrFrequency, $extraArg, $leftxAxisPercentileLimit, $log2OrLinear, $logPath, $pdfPath, $plotAryHsh_ref, $rightxAxisPercentileLimit, $xAxis
#	output: none
#	toCall: &ggplotTwoSampleHistogram($plotAryHsh_ref, $dataPath, $pdfPath, $RScriptPath, $logPath, $leftxAxisPercentileLimit, $rightxAxisPercentileLimit, $xAxis, $binWidth, $dataPtMax, $extraArg, $log2OrLinear, $densityOrFrequency);
#	calledInLine: 2473
#....................................................................................................................................................#

	my ($plotAryHsh_ref, $dataPath, $pdfPath, $RScriptPath, $logPath, $leftxAxisPercentileLimit, $rightxAxisPercentileLimit, $xAxis, $binWidth, $dataPtMax, $extraArg, $log2OrLinear, $densityOrFrequency) = @_;
	
	open (PLOTDATA, ">", $dataPath);
	print PLOTDATA "sample\t$xAxis\n";

	foreach my $sample (keys %{$plotAryHsh_ref}) {
		
		my $plotAry_ref = $plotAryHsh_ref->{$sample};
	
		my $valueStatObj = Statistics::Descriptive::Full->new();
		$valueStatObj->add_data(@{$plotAry_ref});

		my $leftxAxisLimitValue;
		if ($leftxAxisPercentileLimit eq 'min') {
			$leftxAxisLimitValue = $valueStatObj->min();
		} else {
			$leftxAxisLimitValue = $valueStatObj->percentile($leftxAxisPercentileLimit);
		}

		my $rightxAxisLimitValue;
		if ($rightxAxisPercentileLimit eq 'max') {
			$rightxAxisLimitValue = $valueStatObj->max();
		} else {
			$rightxAxisLimitValue = $valueStatObj->percentile($rightxAxisPercentileLimit);
		}
	
		#---trim the end values
		my @trimmedAry = ();
		foreach my $value (@{$plotAry_ref}) {
			push @trimmedAry, $value if $value <= $rightxAxisLimitValue and $value >= $leftxAxisLimitValue;
		}
		
		my $indivDataPtMax = $dataPtMax;
		$indivDataPtMax = @trimmedAry if $indivDataPtMax > @trimmedAry;

		#---down sample the data point number
		my @shuffleIndexAry = shuffle(0..$#trimmedAry);
		foreach my $i (0..$indivDataPtMax-1) {
			my $shuffleValue = $trimmedAry[$shuffleIndexAry[$i]];
			print PLOTDATA "$sample\t$shuffleValue\n";
		}
	}
	close PLOTDATA;

	my $scale = '';
	$scale = ' + scale_x_continuous(trans=log2_trans())' if $log2OrLinear eq 'log2';

	my $plotGeom = "ggplot(dataFrame, aes($xAxis, fill = sample)) + geom_density(alpha = 0.2)";#---default = density
	$plotGeom = "ggplot(dataFrame, aes($xAxis, color = sample)) + geom_freqpoly(binwidth = $binWidth)" if $densityOrFrequency eq 'frequency';

	open (R, ">", $RScriptPath);
	print R "library(ggplot2)"."\n";
	print R "library(scales)"."\n";
	print R "dataFrame = read.table(file=\"$dataPath\", header=TRUE)"."\n";
	print R "$plotGeom + ggtitle(\"Distribution of $xAxis $log2OrLinear scale\") $scale $extraArg"."\n";
	print R "ggsave(file=\"$pdfPath\")\n";
	close R;
	
	system ("R --slave --vanilla --file=$RScriptPath &>$logPath ");

}
sub ggplotXYLinesMultipleSamples {
#....................................................................................................................................................#
#	subroutineCategory: ggplot, plotInR
#	dependOnSub: >none
#	appearInSub: defineMotifPositionBoundsUsingShuffleBackground|926, plotInrMotifProportion|2312, plotMotifPostionFactor|2354, plotTSSScoreAroundmRNAReferencePoint|2391, plotTSSScoreInGeneWithAndWithoutEnd3NAT|2493
#	primaryAppearInSection: >none
#	secondaryAppearInSection: 10_defineMotifBoundary|259, 13_analyzeTSSScore|293, 8_getInitiatorMotif|228
#	input: $RScriptPath, $YAxis, $YVariable, $dataPath, $extraArg, $height, $logPath, $pdfPath, $plotDataHsh_ref, $width, $xAxis
#	output: none
#	toCall: &ggplotXYLinesMultipleSamples($plotDataHsh_ref, $dataPath, $pdfPath, $RScriptPath, $logPath, $xAxis, $YAxis, $YVariable, $extraArg, $height, $width);
#	calledInLine: 974, 2348, 2386, 2427, 2546
#....................................................................................................................................................#

	my ($plotDataHsh_ref, $dataPath, $pdfPath, $RScriptPath, $logPath, $xAxis, $YAxis, $YVariable, $extraArg, $height, $width) = @_;
	
	open (PLOTDATA, ">", $dataPath);
	print PLOTDATA join "", (join "\t", ($YVariable, $YAxis, $xAxis)), "\n";
	foreach my $YCategory (sort keys %{$plotDataHsh_ref}) {
		foreach my $XVal (sort {$a <=> $b} keys %{$plotDataHsh_ref->{$YCategory}}) {
			my $YVal = $plotDataHsh_ref->{$YCategory}{$XVal};
			print PLOTDATA join "", (join "\t", ($YCategory, $YVal, $XVal)), "\n";
		}
	}
	close PLOTDATA;

	open (R, ">", $RScriptPath);
	print R "library(ggplot2)"."\n";
	print R "dataFrame = read.table(file=\"$dataPath\", header=TRUE)"."\n";
	print R "ggplot(dataFrame, aes(x=$xAxis, y=$YAxis, colour=$YVariable)) + geom_line() $extraArg"."\n";
	print R "ggsave(file=\"$pdfPath\", height=$height, width=$width)\n";
	close R;
	
	system ("R --slave --vanilla --file=$RScriptPath &>$logPath");

}
sub investigateTSSRelativeToATGAndTAA {
#....................................................................................................................................................#
#	subroutineCategory: unassigned
#	dependOnSub: reportStatus|3137
#	appearInSub: getmRNAReferencePoints|1943
#	primaryAppearInSection: >none
#	secondaryAppearInSection: 7_retrieveSequenceSurroundingPredefinedTSS|212
#	input: $TSSSearchRngHsh_ref, $mRNAInfoHsh_ref, $validTSSCovPlsPathHsh_ref
#	output: $TSSmRNAEndFreqHsh_ref, $geneBasedTSSInfoHsh_ref
#	toCall: my ($geneBasedTSSInfoHsh_ref, $TSSmRNAEndFreqHsh_ref) = &investigateTSSRelativeToATGAndTAA($validTSSCovPlsPathHsh_ref, $mRNAInfoHsh_ref, $TSSSearchRngHsh_ref);
#	calledInLine: 1969
#....................................................................................................................................................#

	my ($validTSSCovPlsPathHsh_ref, $mRNAInfoHsh_ref, $TSSSearchRngHsh_ref) = @_;
	
	#---generate mRNAID by cntg hash
	my $mRNAIDByCntgHsh_ref = {};
	foreach my $mRNAID (keys %{$mRNAInfoHsh_ref}) {
		my $cntg = $mRNAInfoHsh_ref->{$mRNAID}{'cntg'};
		push @{$mRNAIDByCntgHsh_ref->{$cntg}}, $mRNAID;
	}

	my $TSSmRNAEndFreqHsh_ref = {};
	my $geneBasedTSSInfoHsh_ref = {};

	#---get the exon positions of the goldenmRNA set
	foreach my $cntg (keys %{$mRNAIDByCntgHsh_ref}) {
	
		&reportStatus("Getting TSS frequency around ATG and TAA", 40, "\r");#->3137

		my $filterTEXCntgCovAry_ref = retrieve($validTSSCovPlsPathHsh_ref->{$cntg});

		foreach my $mRNAID (@{$mRNAIDByCntgHsh_ref->{$cntg}}) {
			my $strnd = $mRNAInfoHsh_ref->{$mRNAID}{'strnd'};
			
			#---assume $RNAStart = ATG and RNAEnd = TAA
			my ($RNAStart, $RNAEnd) = @{$mRNAInfoHsh_ref->{$mRNAID}{'RNARng'}};
			my $tmpPosSrchAryHsh_ref = {};
			
			foreach my $rltvPos (-1*$TSSSearchRngHsh_ref->{'ATG'}{'up'}..$TSSSearchRngHsh_ref->{'ATG'}{'dn'}) {
				if ($strnd eq '+') {
					$tmpPosSrchAryHsh_ref->{'ATG'}{$rltvPos} = $RNAStart+$rltvPos;
				} else {
					$tmpPosSrchAryHsh_ref->{'ATG'}{$rltvPos} = $RNAEnd-$rltvPos;
				}
			}

			foreach my $rltvPos (-1*$TSSSearchRngHsh_ref->{'TAA'}{'up'}..$TSSSearchRngHsh_ref->{'TAA'}{'dn'}) {
				if ($strnd eq '+') {
					$tmpPosSrchAryHsh_ref->{'TAA'}{$rltvPos} = $RNAEnd+$rltvPos;
				} else {
					$tmpPosSrchAryHsh_ref->{'TAA'}{$rltvPos} = $RNAStart-$rltvPos;
				}
			}
			
			foreach my $ATGOrTAA (keys %{$tmpPosSrchAryHsh_ref}) {
				
				foreach my $senseOrAntisense ('sense', 'antisense') {
					foreach my $rltvPosOrCov ('rltvPos', 'cov') {
						@{$geneBasedTSSInfoHsh_ref->{$mRNAID}{$ATGOrTAA}{$senseOrAntisense}{$rltvPosOrCov}} = ();
					}
				}
	
				foreach my $rltvPos (keys %{$tmpPosSrchAryHsh_ref->{$ATGOrTAA}}) {
					my $absPos = $tmpPosSrchAryHsh_ref->{$ATGOrTAA}{$rltvPos};

					#---pos > 0 to ensure the search is not out of bound
					if ($absPos > 0 and $filterTEXCntgCovAry_ref->[$absPos-1]) {
						my $i = $absPos-1;
						my ($senseCov, $antisenseCov) = split /,/, $filterTEXCntgCovAry_ref->[$i];
						($antisenseCov, $senseCov) = ($senseCov, $antisenseCov) if $strnd eq '-';

						if ($senseCov > 0) {
							push @{$geneBasedTSSInfoHsh_ref->{$mRNAID}{$ATGOrTAA}{'sense'}{'rltvPos'}}, $rltvPos;
							push @{$geneBasedTSSInfoHsh_ref->{$mRNAID}{$ATGOrTAA}{'sense'}{'cov'}}, $senseCov;
							push @{$geneBasedTSSInfoHsh_ref->{$mRNAID}{$ATGOrTAA}{'sense'}{'absPos'}}, $absPos;
							push @{$TSSmRNAEndFreqHsh_ref->{$ATGOrTAA}{'sense'}}, $rltvPos;
						}

						if ($antisenseCov > 0) {
							push @{$geneBasedTSSInfoHsh_ref->{$mRNAID}{$ATGOrTAA}{'antisense'}{'rltvPos'}}, $rltvPos;
							push @{$geneBasedTSSInfoHsh_ref->{$mRNAID}{$ATGOrTAA}{'antisense'}{'cov'}}, $antisenseCov;
							push @{$geneBasedTSSInfoHsh_ref->{$mRNAID}{$ATGOrTAA}{'antisense'}{'absPos'}}, $absPos;
							push @{$TSSmRNAEndFreqHsh_ref->{$ATGOrTAA}{'antisense'}}, $rltvPos;
						}
					}
				}
				
				foreach my $senseOrAntisense ('sense', 'antisense') {
					if (@{$geneBasedTSSInfoHsh_ref->{$mRNAID}{$ATGOrTAA}{$senseOrAntisense}{'rltvPos'}} == 0) {
						@{$geneBasedTSSInfoHsh_ref->{$mRNAID}{$ATGOrTAA}{$senseOrAntisense}{'rltvPos'}} = (-9999);
						@{$geneBasedTSSInfoHsh_ref->{$mRNAID}{$ATGOrTAA}{$senseOrAntisense}{'cov'}} = (-9999);
					}
				}
			}
		}
	}
	
	return ($geneBasedTSSInfoHsh_ref, $TSSmRNAEndFreqHsh_ref);
}
sub plotAverageTSSScoreForEnd3NATGenes {
#....................................................................................................................................................#
#	subroutineCategory: unassigned
#	dependOnSub: ggplotMultiSampleBoxWhisker|2016
#	appearInSub: >none
#	primaryAppearInSection: 13_analyzeTSSScore|293
#	secondaryAppearInSection: >none
#	input: $TSSScoreggplotDirHsh_ref, $end3NATGeneAvgTSSScoreHsh_ref
#	output: 
#	toCall: &plotAverageTSSScoreForEnd3NATGenes($end3NATGeneAvgTSSScoreHsh_ref, $TSSScoreggplotDirHsh_ref);
#	calledInLine: 308
#....................................................................................................................................................#
	my ($end3NATGeneAvgTSSScoreHsh_ref, $TSSScoreggplotDirHsh_ref) = @_;

	my $plotAryHsh_ref = {};
	
	#@{$plotAryHsh_ref->{'withinWithNATTSS_a'}} = @{$end3NATGeneAvgTSSScoreHsh_ref->{'randWithin'}{'withNATTSS'}{'a'}};
	#@{$plotAryHsh_ref->{'withinWithNATTSS_s'}} = @{$end3NATGeneAvgTSSScoreHsh_ref->{'randWithin'}{'withNATTSS'}{'s'}};
	@{$plotAryHsh_ref->{'TAAWithNATTSS_a'}} = @{$end3NATGeneAvgTSSScoreHsh_ref->{'TAA'}{'withNATTSS'}{'a'}};
	@{$plotAryHsh_ref->{'TAAWithNATTSS_s'}} = @{$end3NATGeneAvgTSSScoreHsh_ref->{'TAA'}{'withNATTSS'}{'s'}};
	@{$plotAryHsh_ref->{'TAANoNAT_a'}} = @{$end3NATGeneAvgTSSScoreHsh_ref->{'TAA'}{'no'}{'a'}};
	@{$plotAryHsh_ref->{'TAANoNAT_s'}} = @{$end3NATGeneAvgTSSScoreHsh_ref->{'TAA'}{'no'}{'s'}};
	#@{$plotAryHsh_ref->{'TAAYesNAT_a'}} = @{$end3NATGeneAvgTSSScoreHsh_ref->{'TAA'}{'yes'}{'a'}};
	#@{$plotAryHsh_ref->{'TAAYesNAT_s'}} = @{$end3NATGeneAvgTSSScoreHsh_ref->{'TAA'}{'yes'}{'s'}};
	@{$plotAryHsh_ref->{'ATGWithNATTSS_a'}} = @{$end3NATGeneAvgTSSScoreHsh_ref->{'ATG'}{'withNATTSS'}{'a'}};
	@{$plotAryHsh_ref->{'ATGWithNATTSS_s'}} = @{$end3NATGeneAvgTSSScoreHsh_ref->{'ATG'}{'withNATTSS'}{'s'}};

	{
		my $nameTag = "end3NATGeneAvgTSSScore.box";
		my $dataPath = "$TSSScoreggplotDirHsh_ref->{dat}/$nameTag.dat";
		my $pdfPath = "$TSSScoreggplotDirHsh_ref->{pdf}/$nameTag.pdf";
		my $RScriptPath = "$TSSScoreggplotDirHsh_ref->{R}/$nameTag.R";
		my $logPath = "$TSSScoreggplotDirHsh_ref->{log}/$nameTag.log";
		my $yAxis = "TSSScore";
		my $log2OrLinear = 'linear';
		my $extraArg = '+  stat_summary(fun.y=mean, geom="line", aes(group=1)) + stat_summary(fun.y=mean, geom="point")';
		my $height = 8;
		my $width = 20;
		my $dataPtMax = 5000;
		&ggplotMultiSampleBoxWhisker($plotAryHsh_ref, $dataPath, $pdfPath, $RScriptPath, $logPath, $yAxis, $extraArg, $log2OrLinear, $dataPtMax, $height, $width);#->2016
	}
	return ();
}
sub plotInrMotifProportion {
#....................................................................................................................................................#
#	subroutineCategory: unassigned
#	dependOnSub: ggplotXYLinesMultipleSamples|2134
#	appearInSub: >none
#	primaryAppearInSection: 8_getInitiatorMotif|228
#	secondaryAppearInSection: >none
#	input: $InrMotifFreqHsh_ref, $cmltvePprtnLimit, $generalggplotDirHsh_ref
#	output: 
#	toCall: &plotInrMotifProportion($InrMotifFreqHsh_ref, $generalggplotDirHsh_ref, $cmltvePprtnLimit);
#	calledInLine: 237
#....................................................................................................................................................#
	my ($InrMotifFreqHsh_ref, $generalggplotDirHsh_ref, $cmltvePprtnLimit) = @_;
	
	my $InrCmltvePprtnPlotDataHsh_ref = {};
	
	#---plot the culmulative proportion
	my $motifNum = 0;
	foreach my $InrSeq (sort {$InrMotifFreqHsh_ref->{'cmltveProportion'}{$a} <=> $InrMotifFreqHsh_ref->{'cmltveProportion'}{$b}} keys %{$InrMotifFreqHsh_ref->{'cmltveProportion'}}) {
		$motifNum++;
		my $cmltveProportion = $InrMotifFreqHsh_ref->{'cmltveProportion'}{$InrSeq};
		$InrCmltvePprtnPlotDataHsh_ref->{'cmltveProportion'}{$motifNum} = $cmltveProportion;
	}
	
	{#--ggplot culmulative pct
		my $nameTag = "InrMotifCumulativeProportion";
		my $plotDataHsh_ref = $InrCmltvePprtnPlotDataHsh_ref;
		my $dataPath = "$generalggplotDirHsh_ref->{dat}/$nameTag.dat";
		my $pdfPath = "$generalggplotDirHsh_ref->{pdf}/$nameTag.pdf";;
		my $RScriptPath = "$generalggplotDirHsh_ref->{R}/$nameTag.R";;
		my $logPath = "$generalggplotDirHsh_ref->{log}/$nameTag.log";;
		my $xAxis = "motifNum";
		my $YAxis = "culmulativeProportion";
		my $YVariable = "category";
		my $extraArg = " + scale_y_continuous(breaks=seq(0, 1, by=0.1))";
		my $height = 6;
		my $width = 12;
		&ggplotXYLinesMultipleSamples($plotDataHsh_ref, $dataPath, $pdfPath, $RScriptPath, $logPath, $xAxis, $YAxis, $YVariable, $extraArg, $height, $width);#->2134
	}

	return ();
}
sub plotMotifPostionFactor {
#....................................................................................................................................................#
#	subroutineCategory: unassigned
#	dependOnSub: ggplotXYLinesMultipleSamples|2134
#	appearInSub: >none
#	primaryAppearInSection: 10_defineMotifBoundary|259
#	secondaryAppearInSection: >none
#	input: $generalggplotDirHsh_ref, $motifPostionFactorHsh_ref
#	output: 
#	toCall: &plotMotifPostionFactor($motifPostionFactorHsh_ref, $generalggplotDirHsh_ref);
#	calledInLine: 266
#....................................................................................................................................................#
	my ($motifPostionFactorHsh_ref, $generalggplotDirHsh_ref) = @_;
	
	my $plotDataHsh_ref = {};
	foreach my $motif (keys %{$motifPostionFactorHsh_ref}) {
		foreach my $pos (keys %{$motifPostionFactorHsh_ref->{$motif}}) {
			$plotDataHsh_ref->{$motif}{$pos} = $motifPostionFactorHsh_ref->{$motif}{$pos};
		}
	}
	
	my $nameTag = "motifPostionFactor";
	my $pdfPath = $generalggplotDirHsh_ref->{'pdf'}."/$nameTag.pdf";
	my $dataPath = $generalggplotDirHsh_ref->{'dat'}."/$nameTag.dat";
	my $RScriptPath = $generalggplotDirHsh_ref->{'R'}."/$nameTag.R";
	my $logPath = $generalggplotDirHsh_ref->{'log'}."/$nameTag.log";
	my $xAxis = 'relativePositon';
	my $YAxis = 'factor';
	my $YVariable = 'motif';
	#my $extraArg = '+ ylim(0, 100)';
	my $extraArg = '';
	my $height = 6;
	my $width = 14;
	&ggplotXYLinesMultipleSamples($plotDataHsh_ref, $dataPath, $pdfPath, $RScriptPath, $logPath, $xAxis, $YAxis, $YVariable, $extraArg, $height, $width);#->2134

	return ();
}
sub plotTSSScoreAroundmRNAReferencePoint {
#....................................................................................................................................................#
#	subroutineCategory: unassigned
#	dependOnSub: ggplotXYLinesMultipleSamples|2134
#	appearInSub: >none
#	primaryAppearInSection: 13_analyzeTSSScore|293
#	secondaryAppearInSection: >none
#	input: $TSSScoreRefPtPlotHsh_ref, $TSSScoreggplotDirHsh_ref, $mRNARefPtHsh_ref
#	output: 
#	toCall: &plotTSSScoreAroundmRNAReferencePoint($TSSScoreRefPtPlotHsh_ref, $TSSScoreggplotDirHsh_ref, $mRNARefPtHsh_ref);
#	calledInLine: 302
#....................................................................................................................................................#
	my ($TSSScoreRefPtPlotHsh_ref, $TSSScoreggplotDirHsh_ref, $mRNARefPtHsh_ref) = @_;
	
	#---[12/10/2013 16:42] choose to plot the s=um of the TSSSocre or just the frequency count
	#foreach my $sumOrCount (qw/sum count/) {
	foreach my $sumOrCount (qw/sum/) {
		foreach my $siteType (keys %{$TSSScoreRefPtPlotHsh_ref}) {
			my $totalNum = keys %{$mRNARefPtHsh_ref->{$siteType}};
			foreach my $mode (keys %{$TSSScoreRefPtPlotHsh_ref->{$siteType}}) {
				my $plotDataHsh_ref = {};
				foreach my $rltvPos (sort keys %{$TSSScoreRefPtPlotHsh_ref->{$siteType}{$mode}}) {
					foreach my $dirtn (sort keys %{$TSSScoreRefPtPlotHsh_ref->{$siteType}{$mode}{$rltvPos}}) {
						$plotDataHsh_ref->{$dirtn}{$rltvPos} = $TSSScoreRefPtPlotHsh_ref->{$siteType}{$mode}{$rltvPos}{$dirtn}{$sumOrCount};
					}
				}
				my $nameTag = "$sumOrCount.$mode.$siteType.TSS.Score";
				my $pdfPath = $TSSScoreggplotDirHsh_ref->{'pdf'}."/$nameTag.pdf";
				my $dataPath = $TSSScoreggplotDirHsh_ref->{'dat'}."/$nameTag.dat";
				my $RScriptPath = $TSSScoreggplotDirHsh_ref->{'R'}."/$nameTag.R";
				my $logPath = $TSSScoreggplotDirHsh_ref->{'log'}."/$nameTag.log";
				my $xAxis = 'relative_Positon';
				my $YAxis = 'average_TSS_Score';
				my $YVariable = 'direction';
				my $extraArg = "+ ggtitle(\"N=$totalNum\")";
				my $height = 6;
				my $width = 14;
				&ggplotXYLinesMultipleSamples($plotDataHsh_ref, $dataPath, $pdfPath, $RScriptPath, $logPath, $xAxis, $YAxis, $YVariable, $extraArg, $height, $width);#->2134
			}
		}
	}
	
	return ();
}
sub plotTSSScoreInExonAndTSS {
#....................................................................................................................................................#
#	subroutineCategory: unassigned
#	dependOnSub: ggplotMultiSampleBoxWhisker|2016, ggplotTwoSampleHistogram|2060
#	appearInSub: >none
#	primaryAppearInSection: 13_analyzeTSSScore|293
#	secondaryAppearInSection: >none
#	input: $TSSScoreExonTSSHshNonZero_ref, $TSSScoreExonTSSHshWithZero_ref, $TSSScoreggplotDirHsh_ref
#	output: 
#	toCall: &plotTSSScoreInExonAndTSS($TSSScoreExonTSSHshNonZero_ref, $TSSScoreExonTSSHshWithZero_ref, $TSSScoreggplotDirHsh_ref);
#	calledInLine: 299
#....................................................................................................................................................#
	my ($TSSScoreExonTSSHshNonZero_ref, $TSSScoreExonTSSHshWithZero_ref, $TSSScoreggplotDirHsh_ref) = @_;
	
	my %tmpHsh = ();
	$tmpHsh{'withZero'}{'plotAryHsh_ref'} = $TSSScoreExonTSSHshWithZero_ref;
	$tmpHsh{'nonZero'}{'plotAryHsh_ref'} = $TSSScoreExonTSSHshNonZero_ref;
	$tmpHsh{'withZero'}{'densityOrFrequency'} = 'density';
	$tmpHsh{'nonZero'}{'densityOrFrequency'} = 'density';
	$tmpHsh{'withZero'}{'log2OrLinear'} = 'linear';
	$tmpHsh{'nonZero'}{'log2OrLinear'} = 'linear';
	
	foreach my $withZeroOrNonZero (keys %tmpHsh) {
		my $plotAryHsh_ref = $tmpHsh{$withZeroOrNonZero}{'plotAryHsh_ref'};
		my $dataPtMax = 3000;

		{
			my $nameTag = "TSSScore.TSS.vs.Exon.histogram.$withZeroOrNonZero";
			my $dataPath = "$TSSScoreggplotDirHsh_ref->{dat}/$nameTag.dat";
			my $pdfPath = "$TSSScoreggplotDirHsh_ref->{pdf}/$nameTag.pdf";
			my $RScriptPath = "$TSSScoreggplotDirHsh_ref->{R}/$nameTag.R";
			my $logPath = "$TSSScoreggplotDirHsh_ref->{log}/$nameTag.log";
			my $leftxAxisPercentileLimit = 'min';
			my $rightxAxisPercentileLimit = 'max';
			my $xAxis = "TSSScore";
			my $binWidth = 0.1;
			my $log2OrLinear = $tmpHsh{$withZeroOrNonZero}{'log2OrLinear'};
			my $extraArg = '';
			my $densityOrFrequency = $tmpHsh{$withZeroOrNonZero}{'densityOrFrequency'};
			&ggplotTwoSampleHistogram($plotAryHsh_ref, $dataPath, $pdfPath, $RScriptPath, $logPath, $leftxAxisPercentileLimit, $rightxAxisPercentileLimit, $xAxis, $binWidth, $dataPtMax, $extraArg, $log2OrLinear, $densityOrFrequency);#->2060
		}
	
		{
			my $nameTag = "TSSScore.TSS.vs.Exon.box.$withZeroOrNonZero";
			my $dataPath = "$TSSScoreggplotDirHsh_ref->{dat}/$nameTag.dat";
			my $pdfPath = "$TSSScoreggplotDirHsh_ref->{pdf}/$nameTag.pdf";
			my $RScriptPath = "$TSSScoreggplotDirHsh_ref->{R}/$nameTag.R";
			my $logPath = "$TSSScoreggplotDirHsh_ref->{log}/$nameTag.log";
			my $yAxis = "TSSScore";
			my $log2OrLinear = $tmpHsh{$withZeroOrNonZero}{'log2OrLinear'};
			my $extraArg = '';
			my $height = 12;
			my $width = 8;
			&ggplotMultiSampleBoxWhisker($plotAryHsh_ref, $dataPath, $pdfPath, $RScriptPath, $logPath, $yAxis, $extraArg, $log2OrLinear, $dataPtMax, $height, $width);#->2016
		}
	}
	return ();
}
sub plotTSSScoreInGeneWithAndWithoutEnd3NAT {
#....................................................................................................................................................#
#	subroutineCategory: unassigned
#	dependOnSub: ggplotXYLinesMultipleSamples|2134
#	appearInSub: >none
#	primaryAppearInSection: 13_analyzeTSSScore|293
#	secondaryAppearInSection: >none
#	input: $TSSScoreIndivmRNAHsh_ref, $TSSScoreggplotDirHsh_ref, $mRNANATInfoHsh_ref, $mRNARefPtHsh_ref
#	output: 
#	toCall: &plotTSSScoreInGeneWithAndWithoutEnd3NAT($TSSScoreIndivmRNAHsh_ref, $TSSScoreggplotDirHsh_ref, $mRNARefPtHsh_ref, $mRNANATInfoHsh_ref);
#	calledInLine: 314
#....................................................................................................................................................#
	my ($TSSScoreIndivmRNAHsh_ref, $TSSScoreggplotDirHsh_ref, $mRNARefPtHsh_ref, $mRNANATInfoHsh_ref) = @_;
	
	my %geneCountHsh = ();
	my %tmpSumHsh = ();
	my $plotDataHsh_ref = {};
	foreach my $mRNAID (keys %{$TSSScoreIndivmRNAHsh_ref->{'TAA'}}) {
		my $withOrWithoutNAT = '';
		
		if (exists $mRNARefPtHsh_ref->{'NAT_TSS'}{$mRNAID}) {
			$withOrWithoutNAT = "withNAT_TSS";
		} elsif ($mRNANATInfoHsh_ref->{$mRNAID}{'end3NAT'} eq 'no') {
			$withOrWithoutNAT = "no";
		} else {
			next;
		}
		
		$geneCountHsh{$withOrWithoutNAT}++;
		
		foreach my $rltvPos (keys %{$TSSScoreIndivmRNAHsh_ref->{'TAA'}{$mRNAID}}) {
			$tmpSumHsh{$withOrWithoutNAT}{$rltvPos} = 0 if not exists $tmpSumHsh{$withOrWithoutNAT}{$rltvPos};
			$tmpSumHsh{$withOrWithoutNAT}{$rltvPos} += $TSSScoreIndivmRNAHsh_ref->{'TAA'}{$mRNAID}{$rltvPos}{'a'};
		}
	}
	
	foreach my $withOrWithoutNAT (keys %tmpSumHsh) {
		foreach my $rltvPos (keys %{$tmpSumHsh{$withOrWithoutNAT}}) {
			$plotDataHsh_ref->{$withOrWithoutNAT}{$rltvPos} = $tmpSumHsh{$withOrWithoutNAT}{$rltvPos}/$geneCountHsh{$withOrWithoutNAT};
		}
	}
	
	{
		my $nameTag = "TSSScore.With.or.Without.End3NAT";
		my $pdfPath = $TSSScoreggplotDirHsh_ref->{'pdf'}."/$nameTag.pdf";
		my $dataPath = $TSSScoreggplotDirHsh_ref->{'dat'}."/$nameTag.dat";
		my $RScriptPath = $TSSScoreggplotDirHsh_ref->{'R'}."/$nameTag.R";
		my $logPath = $TSSScoreggplotDirHsh_ref->{'log'}."/$nameTag.log";
		my $xAxis = 'relative_Positon';
		my $YAxis = 'average_TSS_Score';
		my $YVariable = 'withOrWithout';
		my $extraArg = '';
		my $height = 6;
		my $width = 14;
		&ggplotXYLinesMultipleSamples($plotDataHsh_ref, $dataPath, $pdfPath, $RScriptPath, $logPath, $xAxis, $YAxis, $YVariable, $extraArg, $height, $width);#->2134
	}
	
	return ();
}
sub plotmRNANATInfo {
#....................................................................................................................................................#
#	subroutineCategory: unassigned
#	dependOnSub: ggplotCulmulativeFrequency|1977
#	appearInSub: >none
#	primaryAppearInSection: 6_defineGeneWithEnd3NAT|200
#	secondaryAppearInSection: >none
#	input: $generalggplotDirHsh_ref, $mRNANATInfoHsh_ref
#	output: 
#	toCall: &plotmRNANATInfo($mRNANATInfoHsh_ref, $generalggplotDirHsh_ref);
#	calledInLine: 206
#....................................................................................................................................................#
	my ($mRNANATInfoHsh_ref, $generalggplotDirHsh_ref) = @_;
	
	my %tmpScaleHsh = ('posHitPct'=>'linear', 'covPerNt'=>'log2');
	
	foreach my $posHitPctOrCovPerNt (qw/posHitPct covPerNt/) {
		my $plotAry_ref = ();
		foreach my $mRNAID (keys %{$mRNANATInfoHsh_ref}) {
			push @{$plotAry_ref}, $mRNANATInfoHsh_ref->{$mRNAID}{$posHitPctOrCovPerNt};
		}
		
		{
			my $nameTag = "NATAtmRNA3End.$posHitPctOrCovPerNt.culmulative";
			my $dataPath = "$generalggplotDirHsh_ref->{dat}/$nameTag.dat";
			my $pdfPath = "$generalggplotDirHsh_ref->{pdf}/$nameTag.pdf";;
			my $RScriptPath = "$generalggplotDirHsh_ref->{R}/$nameTag.R";;
			my $logPath = "$generalggplotDirHsh_ref->{log}/$nameTag.log";;
			my $xAxis = "$posHitPctOrCovPerNt";
			my $extraArg = '';
			my $log2OrLinear = $tmpScaleHsh{$posHitPctOrCovPerNt};
			my $height = 6;
			my $width = 12;
			&ggplotCulmulativeFrequency($plotAry_ref, $dataPath, $pdfPath, $RScriptPath, $logPath, $xAxis, $extraArg, $log2OrLinear, $height, $width);#->1977
		}
	}

	return ();
}
sub predictGenomeWideTSS {
#....................................................................................................................................................#
#	subroutineCategory: unassigned
#	dependOnSub: checkRunningThreadAndWaitToJoin|522, createEmptyStorableForGenowideTSSPredictionData|550, generateThreadHshWithRandomCntg|1107, reportStatus|3137
#	appearInSub: >none
#	primaryAppearInSection: 12_predictTSS|282
#	secondaryAppearInSection: >none
#	input: $TSSScoreStorableDir, $cntgMotifHitPlsPathHsh_ref, $fastaHsh_ref, $logTransformPVal, $motifPostionFactorHsh_ref, $predictTSSMotifInfoHsh_ref
#	output: $genomeWideTSSPlsPathHsh_ref
#	toCall: my ($genomeWideTSSPlsPathHsh_ref) = &predictGenomeWideTSS($motifPostionFactorHsh_ref, $cntgMotifHitPlsPathHsh_ref, $predictTSSMotifInfoHsh_ref, $fastaHsh_ref, $TSSScoreStorableDir, $logTransformPVal);
#	calledInLine: 287
#....................................................................................................................................................#
	my ($motifPostionFactorHsh_ref, $cntgMotifHitPlsPathHsh_ref, $predictTSSMotifInfoHsh_ref, $fastaHsh_ref, $TSSScoreStorableDir, $logTransformPVal) = @_;
	
	my ($genomeWideTSSPlsPathHsh_ref, $allGenomeWideTSSPlsPathExist) = &createEmptyStorableForGenowideTSSPredictionData($fastaHsh_ref, $TSSScoreStorableDir);#->550
	
	#----calculate the TSS Score only if not all GenomeWideTSSPlsPath exists
	if ($allGenomeWideTSSPlsPathExist eq 'no') {
		&reportStatus("Start calculating the TSS Score in all cntg", 0, "\n");#->3137

		my $motifInfoHsh_ref = {};
		foreach my $motif (keys %{$motifPostionFactorHsh_ref}) {
			$motifInfoHsh_ref->{$motif}{'postionFactorHsh_ref'} = $motifPostionFactorHsh_ref->{$motif};
			@{$motifInfoHsh_ref->{$motif}{'rltvPosAry'}} = (sort {$a <=> $b} keys %{$motifPostionFactorHsh_ref->{$motif}});
			$motifInfoHsh_ref->{$motif}{'mustValid'} = $predictTSSMotifInfoHsh_ref->{$motif}{'mustValid'};
			$motifInfoHsh_ref->{$motif}{'maxPValGenomePredict'} = $predictTSSMotifInfoHsh_ref->{$motif}{'maxPValGenomePredict'};
			$motifInfoHsh_ref->{$motif}{'scoreFactor'} = $predictTSSMotifInfoHsh_ref->{$motif}{'scoreFactor'};
		}
	
		my $threadToSpawn = 10;
		my @cntgAry = (keys %{$genomeWideTSSPlsPathHsh_ref});
		my $randCntgInThreadHsh_ref = &generateThreadHshWithRandomCntg($threadToSpawn, \@cntgAry);#->1107
		my $cntgProc :shared = 0;
		foreach my $threadNum (sort {$a <=> $b} keys %{$randCntgInThreadHsh_ref}) {
			my $cntgAry_ref = \@{$randCntgInThreadHsh_ref->{$threadNum}};
			my $cntgNum = @{$randCntgInThreadHsh_ref->{$threadNum}};

			&reportStatus("$cntgNum cntgs spawned to thread $threadNum", 0, "\n");#->3137

			#---spawn a new thread
			threads->create(
		
				sub {
					my ($cntgAry_ref) = @_;

						foreach my $cntg (@{$cntgAry_ref}) {
							$cntgProc++;
	
							&reportStatus("$cntgProc cntgs processed", 20, "\r");#->3137

							my $cntgTSSScoreAry_ref = retrieve($genomeWideTSSPlsPathHsh_ref->{$cntg});
							my $cntgMotifHitHsh_ref = retrieve($cntgMotifHitPlsPathHsh_ref->{$cntg});
						
							foreach my $i (0..$#{$cntgTSSScoreAry_ref}) {
								my $pos = $i + 1;
								my %strndScoreHsh = ();
								foreach my $strnd (qw/+ -/) {
									$strndScoreHsh{$strnd} = 0;
									foreach my $motif (keys %{$motifInfoHsh_ref}) {
										my %motifScoreHsh = ();
										$motifScoreHsh{$motif} = 0;
										my $postionFactorHsh_ref = $motifInfoHsh_ref->{$motif}{'postionFactorHsh_ref'};
					
										foreach my $rltvPos (@{$motifInfoHsh_ref->{$motif}{'rltvPosAry'}}) {
											my $srchPos;
											if ($strnd eq '+') {
												$srchPos = $pos + $rltvPos;
											} else {
												$srchPos = $pos - $rltvPos;
											}
						
											if ($cntgMotifHitHsh_ref->{$strnd}{$srchPos}{$motif}) {
												my $pval = $cntgMotifHitHsh_ref->{$strnd}{$srchPos}{$motif};
												if ($pval <= $motifInfoHsh_ref->{$motif}{'maxPValGenomePredict'}) {
													my $postionFactor = $postionFactorHsh_ref->{$rltvPos};
													my $scoreFactor = $motifInfoHsh_ref->{$motif}{'scoreFactor'};
													my $transformPVal = 1/$pval;
													$transformPVal = log($transformPVal) if $logTransformPVal eq 'yes';
													my $score = $scoreFactor*$postionFactor*$transformPVal;
													$motifScoreHsh{$motif} += $score;
												}
											}
										}

										#---collect the score
										$strndScoreHsh{$strnd} += $motifScoreHsh{$motif};
					
										#---if mustValid but not invalid
										if ($motifScoreHsh{$motif} == 0 and $motifInfoHsh_ref->{$motif}{'mustValid'} eq 'yes') {
											$strndScoreHsh{$strnd} = 0;
											last;
										}
									}
								}
			
								if ($strndScoreHsh{'+'} > 0 or $strndScoreHsh{'-'} > 0) {
									$strndScoreHsh{$_} = sprintf "%.5f", $strndScoreHsh{$_} foreach (qw/+ -/);
									$cntgTSSScoreAry_ref->[$i] = join ',', ($strndScoreHsh{'+'}, $strndScoreHsh{'-'});
									#print TMPLOG join "\t", ($strndScoreHsh{'+'}, $strndScoreHsh{'-'}."\n");
								}
							}
						
							store($cntgTSSScoreAry_ref, $genomeWideTSSPlsPathHsh_ref->{$cntg});
						
						}
					}
				,($cntgAry_ref)
			);
		}

		#---wait until all threads are finished
		&checkRunningThreadAndWaitToJoin('yes', 1);#->522

	} else {
		
		&reportStatus("TSS Score storable found. Skip calculating", 0, "\n");#->3137
	
	}
	
	return ($genomeWideTSSPlsPathHsh_ref);
}
sub printAllMotifWiggle {
#....................................................................................................................................................#
#	subroutineCategory: unassigned
#	dependOnSub: printMotifWiggle|2802
#	appearInSub: >none
#	primaryAppearInSection: >none
#	secondaryAppearInSection: >none
#	input: $fullGenomeMotifMastHsh_ref, $resultWigDir
#	output: 
#	toCall: &printAllMotifWiggle($fullGenomeMotifMastHsh_ref, $resultWigDir);
#	calledInLine: none
#....................................................................................................................................................#
	my ($fullGenomeMotifMastHsh_ref, $resultWigDir) = @_;
	
	my %tmpMotifHitRefHsh = ();

	$tmpMotifHitRefHsh{$_} = $fullGenomeMotifMastHsh_ref->{$_} foreach (keys %{$fullGenomeMotifMastHsh_ref});
	
	foreach my $nameTag (keys %tmpMotifHitRefHsh) {
		my $motifPosHsh_ref = $tmpMotifHitRefHsh{$nameTag};
		&printMotifWiggle($resultWigDir, $motifPosHsh_ref, $nameTag);#->2802
	}
	
	return ();
}
sub printCMDLogOrFinishMessage {
#....................................................................................................................................................#
#	subroutineCategory: general
#	dependOnSub: currentTime|694
#	appearInSub: >none
#	primaryAppearInSection: 0_startingTasks|78, 14_finishingTasks|320
#	secondaryAppearInSection: >none
#	input: $CMDLogOrFinishMessage
#	output: none
#	toCall: &printCMDLogOrFinishMessage($CMDLogOrFinishMessage);
#	calledInLine: 84, 329
#....................................................................................................................................................#

	my ($CMDLogOrFinishMessage) = @_;
	
	if ($CMDLogOrFinishMessage eq "CMDLog") {
		#---open a log file if it doesnt exists
		my $absoluteScriptPath = abs_path($0);
		my $dirPath = dirname(rel2abs($absoluteScriptPath));
		my ($scriptName, $callScriptPath, $scriptSuffix) = fileparse($absoluteScriptPath, qr/\.[^.]*/);
		open (CMDLOG, ">>$dirPath/$scriptName.cmd.log.txt"); #---append the CMD log file
		print CMDLOG "[".&currentTime()."]\t"."$dirPath/$scriptName$scriptSuffix ".(join " ", @ARGV)."\n";#->694
		close CMDLOG;
		print "\n=========================================================================\n";
		print "[".&currentTime()."] starts running ...... \n";#->694
		print "=========================================================================\n\n";

	} elsif ($CMDLogOrFinishMessage eq "finishMessage") {
		print "\n=========================================================================\n";
		print "[".&currentTime()."] finished running .......\n";#->694
		print "=========================================================================\n\n";
	}
}
sub printGeneAboveTSSScoreLog {
#....................................................................................................................................................#
#	subroutineCategory: unassigned
#	dependOnSub: >none
#	appearInSub: >none
#	primaryAppearInSection: 13_analyzeTSSScore|293
#	secondaryAppearInSection: >none
#	input: $TSSScoreLogDir, $geneAboveTSSScoreHsh_ref
#	output: 
#	toCall: &printGeneAboveTSSScoreLog($geneAboveTSSScoreHsh_ref, $TSSScoreLogDir);
#	calledInLine: 312
#....................................................................................................................................................#
	my ($geneAboveTSSScoreHsh_ref, $TSSScoreLogDir) = @_;
	
	my $geneAboveTSSScoreLogPath = "$TSSScoreLogDir/geneAboveTSSScore.xls";
	open LOG, ">", $geneAboveTSSScoreLogPath;
	foreach my $refPosType (sort keys %{$geneAboveTSSScoreHsh_ref}) {
		foreach my $category (sort qw/withNAT_TSS no/) {
			foreach my $dirtn (sort keys %{$geneAboveTSSScoreHsh_ref->{$refPosType}{$category}}) {
				foreach my $cutoff (sort {$a <=> $b} keys %{$geneAboveTSSScoreHsh_ref->{$refPosType}{$category}{$dirtn}}) {
					my $geneNum = $geneAboveTSSScoreHsh_ref->{$refPosType}{$category}{$dirtn}{$cutoff}{'geneNum'};
					my $numAbvCutoff = $geneAboveTSSScoreHsh_ref->{$refPosType}{$category}{$dirtn}{$cutoff}{'numAbvCutoff'};
					my $prprtnAbvCutoff = $geneAboveTSSScoreHsh_ref->{$refPosType}{$category}{$dirtn}{$cutoff}{'prprtnAbvCutoff'};
					print LOG join "", ((join "\t", ($refPosType, $category, $dirtn, $cutoff, $geneNum, $numAbvCutoff, $prprtnAbvCutoff)), "\n");
				}
			}
		}
	}
	close LOG;
	
	return ();
}
sub printMotifWiggle {
#....................................................................................................................................................#
#	subroutineCategory: unassigned
#	dependOnSub: >none
#	appearInSub: printAllMotifWiggle|2712
#	primaryAppearInSection: >none
#	secondaryAppearInSection: >none
#	input: $motifPosHsh_ref, $nameTag, $resultWigDir
#	output: 
#	toCall: &printMotifWiggle($resultWigDir, $motifPosHsh_ref, $nameTag);
#	calledInLine: 2731
#....................................................................................................................................................#
	my ($resultWigDir, $motifPosHsh_ref, $nameTag) = @_;
	
	my %tmpFHHsh = ();
	my %strndWordHsh = ('+'=>'plus', '-'=>'minus');
	foreach my $strnd (qw/+ -/) {
		my $path = "$resultWigDir/$nameTag.$strndWordHsh{$strnd}.wig";
		open $tmpFHHsh{$strnd}, ">", $path;
	}
	
	foreach my $cntg (sort keys %{$motifPosHsh_ref}) {
		foreach my $strnd (sort keys %{$motifPosHsh_ref->{$cntg}}) {
			print {$tmpFHHsh{$strnd}} "variableStep chrom=$cntg span=1\n";
			foreach my $pos (sort keys %{$motifPosHsh_ref->{$cntg}{$strnd}}) {
				print {$tmpFHHsh{$strnd}} join '', ((join "\t", ($pos, 1)), "\n");
			}
		}
	}
	
	return ();
}
sub printTSSScoreWiggle {
#....................................................................................................................................................#
#	subroutineCategory: unassigned
#	dependOnSub: checkRunningThreadAndWaitToJoin|522, printWiggleSingleTrackFromCntgCovPlsPathHsh|2870, reportStatus|3137
#	appearInSub: >none
#	primaryAppearInSection: 12_predictTSS|282
#	secondaryAppearInSection: >none
#	input: $TSSScoreWigDir, $genomeWideTSSPlsPathHsh_ref
#	output: 
#	toCall: &printTSSScoreWiggle($TSSScoreWigDir, $genomeWideTSSPlsPathHsh_ref);
#	calledInLine: 288
#....................................................................................................................................................#
	my ($TSSScoreWigDir, $genomeWideTSSPlsPathHsh_ref) = @_;
	
	my %tmpStrndInfoHsh = ();
	$tmpStrndInfoHsh{'+'}{'wigPath'} = "$TSSScoreWigDir/genomeWideTSSScore.plus.wig";
	$tmpStrndInfoHsh{'+'}{'aryIndex'} = 0;
	$tmpStrndInfoHsh{'-'}{'wigPath'} = "$TSSScoreWigDir/genomeWideTSSScore.minus.wig";
	$tmpStrndInfoHsh{'-'}{'aryIndex'} = 1;
	
	&reportStatus("Printing TSS Score Wiggle", 0, "\n");#->3137

	foreach my $strnd (keys %tmpStrndInfoHsh) {
		my $wigPath = $tmpStrndInfoHsh{$strnd}{'wigPath'};
		my $aryIndex = $tmpStrndInfoHsh{$strnd}{'aryIndex'};
		my $gzip = 'no';
		my $cntgCovPlsPathHsh_ref = $genomeWideTSSPlsPathHsh_ref;
		if (not -s $wigPath) {
			threads->create(\&printWiggleSingleTrackFromCntgCovPlsPathHsh, ($cntgCovPlsPathHsh_ref, $aryIndex, $wigPath, $gzip));#->2870
		}
	}

	&checkRunningThreadAndWaitToJoin('yes', 1);#->522

	return ();
}
sub printWiggleSingleTrackFromCntgCovPlsPathHsh {
#....................................................................................................................................................#
#	subroutineCategory: wiggle
#	dependOnSub: >none
#	appearInSub: printTSSScoreWiggle|2834
#	primaryAppearInSection: >none
#	secondaryAppearInSection: 12_predictTSS|282
#	input: $aryIndex, $cntgCovPlsPathHsh_ref, $gzip, $wigPath
#	output: none
#	toCall: &printWiggleSingleTrackFromCntgCovPlsPathHsh($cntgCovPlsPathHsh_ref, $aryIndex, $wigPath, $gzip);
#	calledInLine: 2861
#....................................................................................................................................................#
	
	my ($cntgCovPlsPathHsh_ref, $aryIndex, $wigPath, $gzip) = @_;
	
	open (WIGGLE, ">", $wigPath);
	
	foreach my $cntg (sort keys %{$cntgCovPlsPathHsh_ref}) {

		print WIGGLE "variableStep chrom=$cntg span=1\n";

		my $cntgCovPlsPath = "$cntgCovPlsPathHsh_ref->{$cntg}";
 		system ("gzip -df $cntgCovPlsPath.gz") if (-s "$cntgCovPlsPath.gz" and $gzip eq 'yes');
		my $cntgCovAry_ref = retrieve($cntgCovPlsPath);
		system ("gzip -f $cntgCovPlsPath") if (-s $cntgCovPlsPath and $gzip eq 'yes');
		for my $i (0..$#{$cntgCovAry_ref}) {
			if ($cntgCovAry_ref->[$i]) {
				my @tmpCovAry = split /,/, $cntgCovAry_ref->[$i];
				my $cov = $tmpCovAry[$aryIndex];
				if ($cov > 0) {
					my $pos = $i + 1;
					print WIGGLE join '', ((join "\t", ($pos, $cov)), "\n");
				}
			}
		}
	}
	close WIGGLE;
}
sub printmRNANATInfoLog {
#....................................................................................................................................................#
#	subroutineCategory: unassigned
#	dependOnSub: >none
#	appearInSub: >none
#	primaryAppearInSection: 6_defineGeneWithEnd3NAT|200
#	secondaryAppearInSection: >none
#	input: $mRNANATInfoHsh_ref, $mRNANATInfoLogPath
#	output: none
#	toCall: &printmRNANATInfoLog($mRNANATInfoHsh_ref, $mRNANATInfoLogPath);
#	calledInLine: 207
#....................................................................................................................................................#

	my ($mRNANATInfoHsh_ref, $mRNANATInfoLogPath) = @_;
	
	open (LOG, ">", $mRNANATInfoLogPath);
	foreach my $mRNAID (sort keys %{$mRNANATInfoHsh_ref}) {
		my @outputAry = ();
		push @outputAry, 'mRNAID';
		foreach my $item (sort keys %{$mRNANATInfoHsh_ref->{$mRNAID}}) {
			push @outputAry, $item;
		}
		print LOG join "", ((join "\t", @outputAry), "\n");
		last;
	}
	
	foreach my $mRNAID (sort keys %{$mRNANATInfoHsh_ref}) {
		my @outputAry = ();
		push @outputAry, $mRNAID;
		foreach my $item (sort keys %{$mRNANATInfoHsh_ref->{$mRNAID}}) {
			push @outputAry, $mRNANATInfoHsh_ref->{$mRNAID}{$item};
		}
		print LOG join "", ((join "\t", @outputAry), "\n");
	}
	
	close LOG;
}
sub readGFF_oneRNAPerGene {
#....................................................................................................................................................#
#	subroutineCategory: general, gff
#	dependOnSub: currentTime|694
#	appearInSub: >none
#	primaryAppearInSection: 5_processInputData|174
#	secondaryAppearInSection: >none
#	input: $gffPath
#	output: $geneInfoHsh_ref
#	toCall: my ($geneInfoHsh_ref) = &readGFF_oneRNAPerGene($gffPath);
#	calledInLine: 186
#....................................................................................................................................................#

	my ($gffPath) = @_;

	my $geneInfoHsh_ref = {};
	
	#---read the gff
	my $geneByRNAHsh_ref = {};

	open (GFF, $gffPath);
	print "[".&currentTime()."] Reading: $gffPath\n";#->694
	while (my $theLine = <GFF>) {

		chomp $theLine;
		
		last if $theLine =~ m/^##FASTA/;
		
		if ($theLine !~ m/^\#|^\@/ and $theLine !~ m/\tsupercontig\t/) {

			my ($seq, undef, $geneCategory, $featureStart, $featureEnd, undef, $geneStrd, undef, $dscrptns) = split (/\t/, $theLine);
			
			#----assigne all non -/+ will be treated as plus
			$geneStrd = "+" if (($geneStrd ne "-") and ($geneStrd ne "+"));
			
			my @dscrptnsSplt = split /;/, $dscrptns;
			my ($unqID, $parent);
			my $geneName = "unknown";
			foreach my $theDscptn (@dscrptnsSplt) {
				if ($theDscptn =~ m/^ID=/) {$unqID = substr ($theDscptn, index ($theDscptn, "=")+1);}
				if ($theDscptn =~ m/^Parent=/) {$parent = substr ($theDscptn, index ($theDscptn, "=")+1);}
				if ($theDscptn =~ m/^description=/) {$geneName = substr ($theDscptn, index ($theDscptn, "=")+1);}
			}

			if ($geneCategory eq "gene") {#---gene
				
				my $geneID = $unqID;
				
				$geneInfoHsh_ref->{$geneID}{'strnd'} = $geneStrd;
				$geneInfoHsh_ref->{$geneID}{'cntg'} = $seq;
				$geneInfoHsh_ref->{$geneID}{'description'} = uri_unescape($geneName);
				$geneInfoHsh_ref->{$geneID}{'description'} =~ s/\+/ /g;
				@{$geneInfoHsh_ref->{$geneID}{'geneRng'}} = ($featureStart, $featureEnd);

			} elsif ($geneCategory eq "CDS") {#---Only for coding genes
				
				my $RNAID = $parent;
				my $geneID = $geneByRNAHsh_ref->{$RNAID};
				push @{$geneInfoHsh_ref->{$geneID}{'CDSRng'}}, ($featureStart, $featureEnd);
				
			} elsif ($geneCategory eq "exon") {#---exon, may be exons of alternative transcripts, wiull sort out later
				my $RNAID = $parent;
				my $geneID = $geneByRNAHsh_ref->{$RNAID};
				push @{$geneInfoHsh_ref->{$geneID}{'exonRng'}}, ($featureStart, $featureEnd);
				
			} else {#---can be tRNA, rRNA, mRNA, repRNA, ncRNA
				my $RNAID = $unqID;
				my $geneID = $parent;
				$geneByRNAHsh_ref->{$RNAID} = $geneID;
				$geneInfoHsh_ref->{$geneID}{'ctgry'} = $geneCategory;
				@{$geneInfoHsh_ref->{$geneID}{'RNARng'}} = ($featureStart, $featureEnd);
				$geneInfoHsh_ref->{$geneID}{'RNAID'} = $RNAID;
			}
		}#---end of if (($theLine !~ m/^\#|^\@/) and ($theLine !~ m/\tsupercontig\t/)) {
	}#---end of while (my $theLine = <INFILE>)
	close GFF;
	
	#---get the UTR if any
	my $minUTRLength = 10;
	foreach my $geneID (keys %{$geneInfoHsh_ref}) {

		if (exists $geneInfoHsh_ref->{$geneID}{'CDSRng'}) {
			my $exonMin = min(@{$geneInfoHsh_ref->{$geneID}{'exonRng'}});
			my $exonMax = max(@{$geneInfoHsh_ref->{$geneID}{'exonRng'}});
			my $CDSMin = min(@{$geneInfoHsh_ref->{$geneID}{'CDSRng'}});
			my $CDSMax = max(@{$geneInfoHsh_ref->{$geneID}{'CDSRng'}});

			if ($geneInfoHsh_ref->{$geneID}{'strnd'} eq '+') {
				@{$geneInfoHsh_ref->{$geneID}{'UTR5Rng'}} = ($exonMin, $CDSMin-1) if ($CDSMin-$exonMin > $minUTRLength);
				@{$geneInfoHsh_ref->{$geneID}{'UTR3Rng'}} = ($CDSMax+1, $exonMax) if ($exonMax-$CDSMax > $minUTRLength);
			} else {
				@{$geneInfoHsh_ref->{$geneID}{'UTR3Rng'}} = ($exonMin, $CDSMin-1) if ($CDSMin-$exonMin > $minUTRLength);
				@{$geneInfoHsh_ref->{$geneID}{'UTR5Rng'}} = ($CDSMax+1, $exonMax) if ($exonMax-$CDSMax > $minUTRLength);
			}
		}
	}
	
	return ($geneInfoHsh_ref);
}
sub readMultiFasta {
#....................................................................................................................................................#
#	subroutineCategory: fasta, general
#	dependOnSub: reportStatus|3137
#	appearInSub: >none
#	primaryAppearInSection: 5_processInputData|174
#	secondaryAppearInSection: >none
#	input: $fastaPath
#	output: $fastaHsh_ref
#	toCall: my ($fastaHsh_ref) = &readMultiFasta($fastaPath);
#	calledInLine: 181
#....................................................................................................................................................#

	my ($fastaPath) = @_;

	my ($seq, $seqName);
	my $fastaHsh_ref = {};
	my $i = 0;

	&reportStatus("Reading: $fastaPath", 0, "\n");#->3137
	
	open (INFILE, $fastaPath);
	chomp (my $curntLine = <INFILE>); #get the first line
	while (my $nextLine = <INFILE>) {
		chomp $nextLine;
		
		#---Only two types of line in current line, the header or seq
		if ($curntLine =~ m/^>/) {#-- header line
			my @theLineSplt = split (/\|/, $curntLine);
			$seqName = $theLineSplt[0]; #---get the first tag
			$seqName =~ s/ //g; #---remove space
			$seqName =~ s/>//g; #---remove space
		} else {#--seq line
			$seq = $seq.$curntLine;
		}
		
		#---check if next line has a > or that's the end of file
		if ($nextLine =~ m/^>/) {
			$seq =~ tr/a-z/A-Z/;
			$fastaHsh_ref->{$seqName} = $seq;
			$seq = "";

			#---ad hoc limit
			#my $cntgNum = keys %{$fastaHsh_ref}; last if $cntgNum > 100;

		} elsif (eof(INFILE)) {#---this is the last line
			$seq =~ tr/a-z/A-Z/;
			$seq = $seq.$nextLine;
			$fastaHsh_ref->{$seqName} = $seq;
		}
		
		#---next line becomes current line
		$curntLine = $nextLine;
	}

	close INFILE;
	return ($fastaHsh_ref);
}
sub readParameters {
#....................................................................................................................................................#
#	subroutineCategory: general
#	dependOnSub: >none
#	appearInSub: >none
#	primaryAppearInSection: 0_startingTasks|78
#	secondaryAppearInSection: >none
#	input: none
#	output: $fastaPath, $fullReadDefineNATPileupIndxPath, $gffPath, $outDir, $validTSSCntgPlsPath
#	toCall: my ($fullReadDefineNATPileupIndxPath, $validTSSCntgPlsPath, $fastaPath, $gffPath, $outDir) = &readParameters();
#	calledInLine: 87
#....................................................................................................................................................#
	
	my ($fullReadDefineNATPileupIndxPath, $validTSSCntgPlsPath, $fastaPath, $gffPath, $outDir);

	my $dirPath = dirname(rel2abs($0));
	$outDir = "$dirPath/TSSMotifAnalyzer/";

	GetOptions 	("fullReadDefineNATPileupIndxPath=s" => \$fullReadDefineNATPileupIndxPath,
				 "validTSSCntgPlsPath=s"  => \$validTSSCntgPlsPath,
				 "fastaPath=s"  => \$fastaPath,
				 "gffPath=s"  => \$gffPath,
				 "outDir:s"  => \$outDir)

	or die	("Error in command line arguments\n");
	
	#---check file
	foreach my $fileToCheck ($fullReadDefineNATPileupIndxPath, $validTSSCntgPlsPath, $fastaPath, $gffPath) {
		die "Can't read $fileToCheck" if (not -s $fileToCheck and not -s "$fileToCheck.gz");
	}

	system "mkdir -p -m 777 $outDir/";

	return($fullReadDefineNATPileupIndxPath, $validTSSCntgPlsPath, $fastaPath, $gffPath, $outDir);
}
sub reportStatus {
#....................................................................................................................................................#
#	subroutineCategory: general
#	dependOnSub: currentTime|694
#	appearInSub: calculateBackgroundNucleotideFrequency|402, checkGeneInfo|496, checkRunningThreadAndWaitToJoin|522, createEmptyStorableForGenowideTSSPredictionData|550, createMotifHitHshStorable|628, defineEnd3NATInmRNA|746, generateShuffleSeq|1045, getAverageTSSScoreForEnd3NATGenes|1134, getBaseAtTSSAndExon|1182, getCtgryGeneInfo|1246, getIndivCntgCovPlsPath|1331, getInrMotif|1364, getProportionOfGenesWithTSSScoreAboveCutOff|1530, getSequenceAroundmRNAReferencePoint|1589, getTSSScoreAroundmRNAReferencePoint|1692, getTSSScoreForAllGenes|1784, getTSSScoreInExonAndTSS|1884, investigateTSSRelativeToATGAndTAA|2169, predictGenomeWideTSS|2591, printTSSScoreWiggle|2834, readMultiFasta|3044, scanMotifAroundSiteWithMAST|3210, scanMotifWholeGenomeWithMAST|3277, zipUnzipCntgCovInPlsPathHsh|3349
#	primaryAppearInSection: >none
#	secondaryAppearInSection: 12_predictTSS|282, 13_analyzeTSSScore|293, 14_finishingTasks|320, 5_processInputData|174, 6_defineGeneWithEnd3NAT|200, 7_retrieveSequenceSurroundingPredefinedTSS|212, 8_getInitiatorMotif|228, 9_scanMotifOccurence|246
#	input: $lineEnd, $message, $numTrailingSpace
#	output: 
#	toCall: &reportStatus($message, $numTrailingSpace, $lineEnd);
#	calledInLine: 420, 425, 510, 518, 545, 577, 645, 764, 786, 836, 1066, 1075, 1151, 1202, 1265, 1278, 1359, 1382, 1550, 1616, 1727, 1808, 1821, 1902, 2196, 2608, 2627, 2638, 2705, 2853, 3062, 3228, 3242, 3296, 3299, 3308, 3311, 3314, 3364
#....................................................................................................................................................#
	my ($message, $numTrailingSpace, $lineEnd) = @_;

	my $trailingSpaces = '';
	$trailingSpaces .= " " for (1..$numTrailingSpace);
	
	print "[".&currentTime()."] ".$message.$trailingSpaces.$lineEnd;#->694

	return ();
}
sub reverseComplementRefFasta {
#....................................................................................................................................................#
#	subroutineCategory: unassigned
#	dependOnSub: >none
#	appearInSub: >none
#	primaryAppearInSection: 5_processInputData|174
#	secondaryAppearInSection: >none
#	input: $fastaHsh_ref, $resultFastaDir
#	output: $fastaWithRevComHsh_ref, $revComFastaPath
#	toCall: my ($fastaWithRevComHsh_ref, $revComFastaPath) = &reverseComplementRefFasta($fastaHsh_ref, $resultFastaDir);
#	calledInLine: 182
#....................................................................................................................................................#
	my ($fastaHsh_ref, $resultFastaDir) = @_;
	
	my $fastaWithRevComHsh_ref = {};
	my $revComFastaPath = "$resultFastaDir/ref.with.rev.com.fasta";
	open FASTA, ">", $revComFastaPath;
	foreach my $seqName (sort keys %{$fastaHsh_ref}) {
		my $revComSeq = reverse $fastaHsh_ref->{$seqName};
		$revComSeq =~ tr/ACGTacgt/TGCAtgca/;
		print FASTA ">$seqName\+\n";
		print FASTA "$fastaHsh_ref->{$seqName}\n";
		print FASTA ">$seqName\-\n";
		print FASTA "$revComSeq\n";
		$fastaWithRevComHsh_ref->{$seqName}{'+'} = $fastaHsh_ref->{$seqName};
		$fastaWithRevComHsh_ref->{$seqName}{'-'} = $revComSeq;
	}
	close FASTA;

	return ($fastaWithRevComHsh_ref, $revComFastaPath);
}
sub runMAST {
#....................................................................................................................................................#
#	subroutineCategory: thridPartyApp
#	dependOnSub: >none
#	appearInSub: scanMotifAroundSiteWithMAST|3210, scanMotifWholeGenomeWithMAST|3277
#	primaryAppearInSection: >none
#	secondaryAppearInSection: 9_scanMotifOccurence|246
#	input: $bfilePath, $extraOption, $fastaPath, $mastOutDir, $motifFilePath
#	output: $mastHitLog
#	toCall: my ($mastHitLog) = &runMAST($fastaPath, $motifFilePath, $bfilePath, $mastOutDir, $extraOption);
#	calledInLine: 3252, 3263, 3309
#....................................................................................................................................................#
	my ($fastaPath, $motifFilePath, $bfilePath, $mastOutDir, $extraOption) = @_;
	
	my $mastHitLog = "$mastOutDir/mast.hit.txt";
	my $mastErrorLog = "$mastOutDir/mast.error.txt";
	my $mastCmd = "mast $motifFilePath $fastaPath -oc $mastOutDir -ev 1e+10 -mt 10 -hit_list $extraOption >$mastHitLog 2>$mastErrorLog";
	system ("$mastCmd");

	return ($mastHitLog);
}
sub scanMotifAroundSiteWithMAST {
#....................................................................................................................................................#
#	subroutineCategory: unassigned
#	dependOnSub: getSingleMotifMASTLogPostionalData|1657, reportStatus|3137, runMAST|3189
#	appearInSub: >none
#	primaryAppearInSection: 9_scanMotifOccurence|246
#	secondaryAppearInSection: >none
#	input: $bkgdNtFreqHsh_ref, $mastRunDir, $motifFilePathHsh_ref, $predictTSSMotifInfoHsh_ref, $resultStorableDir, $seqAroundSiteInfoHsh_ref, $shuffleSeqAroundSiteInfoHsh_ref
#	output: $mastAroundSiteResultHsh_ref
#	toCall: my ($mastAroundSiteResultHsh_ref) = &scanMotifAroundSiteWithMAST($seqAroundSiteInfoHsh_ref, $shuffleSeqAroundSiteInfoHsh_ref, $motifFilePathHsh_ref, $bkgdNtFreqHsh_ref, $mastRunDir, $predictTSSMotifInfoHsh_ref, $resultStorableDir);
#	calledInLine: 253
#....................................................................................................................................................#
	my ($seqAroundSiteInfoHsh_ref, $shuffleSeqAroundSiteInfoHsh_ref, $motifFilePathHsh_ref, $bkgdNtFreqHsh_ref, $mastRunDir, $predictTSSMotifInfoHsh_ref, $resultStorableDir) = @_;

	my $mastAroundSiteResultHshPlsPath = "$resultStorableDir/mastAroundSiteResultHsh.pls";
	my $mastAroundSiteResultHsh_ref = {};
	
	if (-s $mastAroundSiteResultHshPlsPath) {
		
		&reportStatus('Retrieving mastAroundSiteResultHshPlsPath. Skip mast around mRNA_TSS', 0, "\n");#->3137
		
		$mastAroundSiteResultHsh_ref = retrieve($mastAroundSiteResultHshPlsPath);

	} else {
		#my $extraOption = ' -best '; #---extrat option in MAST
		#my $extraOption = ' -comp '; #---extrat option in MAST
		my $extraOption = ' -norc ';#---extrat option in MAST
		#foreach my $siteType (keys %{$seqAroundSiteInfoHsh_ref}) {
		foreach my $siteType ('mRNA_TSS') {#----do the mRNA_TSS onl,
			foreach my $motif (keys %{$predictTSSMotifInfoHsh_ref}) {
				my $motifFilePath = $motifFilePathHsh_ref->{$motif};
				my $maxHitPVal = $predictTSSMotifInfoHsh_ref->{$motif}{'maxPValDefineBound'};
	
				&reportStatus("Running mast for $motif around $siteType", 0, "\n");#->3137
				
				my $maxPos = $seqAroundSiteInfoHsh_ref->{$siteType}{'length'};
				my $totalSeqNum = $seqAroundSiteInfoHsh_ref->{$siteType}{'totalSeqNum'};
	
				{#---on query seq
					my $fastaPath = $seqAroundSiteInfoHsh_ref->{$siteType}{'fastaPath'};
					my $subMastOutDir = "$mastRunDir/$motif/query/$siteType/";
					system ("mkdir -pm 777 $subMastOutDir");
					my $bfilePath = $bkgdNtFreqHsh_ref->{$siteType}{'full'}{'bfilePath'};
					my $mastHitLog = &runMAST($fastaPath, $motifFilePath, $bfilePath, $subMastOutDir, $extraOption);#->3189
					my ($motifPctHsh_ref, $tmpHitBySeqHsh_ref) = &getSingleMotifMASTLogPostionalData($mastHitLog, $maxHitPVal, $maxPos, $totalSeqNum);#->1657
					$mastAroundSiteResultHsh_ref->{$siteType}{$motif}{'query'}{'tmpHitBySeqHsh_ref'} = $tmpHitBySeqHsh_ref;#---take the has out of the lexcial scope
					$mastAroundSiteResultHsh_ref->{$siteType}{$motif}{'query'}{'motifPctHsh_ref'} = $motifPctHsh_ref;#---take the has out of the lexcial scope
				}
	
				{#---on shuffle seq
					my $fastaPath = $shuffleSeqAroundSiteInfoHsh_ref->{$siteType}{'full'}{'fastaPath'};
					my $subMastOutDir = "$mastRunDir/$motif/shuffle/$siteType/";
					system ("mkdir -pm 777 $subMastOutDir");
					my $bfilePath = $bkgdNtFreqHsh_ref->{$siteType}{'full'}{'bfilePath'};
					my $mastHitLog = &runMAST($fastaPath, $motifFilePath, $bfilePath, $subMastOutDir, $extraOption);#->3189
					my ($motifPctHsh_ref, $tmpHitBySeqHsh_ref) = &getSingleMotifMASTLogPostionalData($mastHitLog, $maxHitPVal, $maxPos, $totalSeqNum);#->1657
					$mastAroundSiteResultHsh_ref->{$siteType}{$motif}{'shuffle'}{'tmpHitBySeqHsh_ref'} = $tmpHitBySeqHsh_ref;#---take the has out of the lexcial scope
					$mastAroundSiteResultHsh_ref->{$siteType}{$motif}{'shuffle'}{'motifPctHsh_ref'} = $motifPctHsh_ref;#---take the has out of the lexcial scope
				}
			}
		}
		
		store($mastAroundSiteResultHsh_ref, $mastAroundSiteResultHshPlsPath);
	}
	
	return ($mastAroundSiteResultHsh_ref);
}
sub scanMotifWholeGenomeWithMAST {
#....................................................................................................................................................#
#	subroutineCategory: unassigned
#	dependOnSub: createMotifHitHshStorable|628, generateMASTBackgroundFile|1000, getMastGenomeBothStrandHit|1415, reportStatus|3137, runMAST|3189, storeMotifHitToHshStorable|3324
#	appearInSub: >none
#	primaryAppearInSection: 9_scanMotifOccurence|246
#	secondaryAppearInSection: >none
#	input: $fastaHsh_ref, $fastaLengthHsh_ref, $mastBackgroundDir, $mastRunDir, $maxPolymerSize, $motifFilePathHsh_ref, $resultStorableDir, $revComFastaPath
#	output: $cntgMotifHitPlsPathHsh_ref
#	toCall: my ($cntgMotifHitPlsPathHsh_ref) = &scanMotifWholeGenomeWithMAST($revComFastaPath, $fastaHsh_ref, $motifFilePathHsh_ref, $mastRunDir, $maxPolymerSize, $mastBackgroundDir, $fastaLengthHsh_ref, $resultStorableDir);
#	calledInLine: 251
#....................................................................................................................................................#
	my ($revComFastaPath, $fastaHsh_ref, $motifFilePathHsh_ref, $mastRunDir, $maxPolymerSize, $mastBackgroundDir, $fastaLengthHsh_ref, $resultStorableDir) = @_;
	
	my ($cntgMotifHitIdxPlsPath, $cntgMotifHitPlsPathHsh_ref, $allCntgMotifHitPlsExist) = &createMotifHitHshStorable($fastaHsh_ref, $resultStorableDir);#->628

	#---will do the prediction only when cntgMotifHitIdxPls doesnt exist, checked in &createMotifHitHshStorable
	if ($allCntgMotifHitPlsExist eq 'no') {
		my $bfilePath = "$mastBackgroundDir/full.genome.freq.txt";
		if (not -s $bfilePath) {
			&reportStatus("Generating full genome background file", 0, "\n");#->3137
			my (undef, undef) = &generateMASTBackgroundFile($fastaHsh_ref, $bfilePath, $maxPolymerSize);#->1000
		} else {
			&reportStatus("Full genome background file found", 0, "\n");#->443	#->3137
		}

		foreach my $motif (keys %{$motifFilePathHsh_ref}) {
			my $motifFilePath = $motifFilePathHsh_ref->{$motif};
			my $subMastOutDir = "$mastRunDir/$motif/query/fullGenome/";
			system ("mkdir -pm 777 $subMastOutDir");
			my $extraOption = ' -norc ';

			&reportStatus("Running mast for $motif on full genome", 0, "\n");#->3137
			my ($mastHitLog) = &runMAST($revComFastaPath, $motifFilePath, $bfilePath, $subMastOutDir, $extraOption);#->3189

			&reportStatus("Getting mast results for $motif", 0, "\n");#->3137
			my ($allHitBySeqHsh_ref) = &getMastGenomeBothStrandHit($mastHitLog, $fastaLengthHsh_ref);#->1415

			&reportStatus("Storing mast results for $motif", 0, "\n");#->3137
			&storeMotifHitToHshStorable($cntgMotifHitPlsPathHsh_ref, $allHitBySeqHsh_ref, $motif);#->3324

			$allHitBySeqHsh_ref = {};
		}
	}

	return ($cntgMotifHitPlsPathHsh_ref);
}
sub storeMotifHitToHshStorable {
#....................................................................................................................................................#
#	subroutineCategory: unassigned
#	dependOnSub: >none
#	appearInSub: scanMotifWholeGenomeWithMAST|3277
#	primaryAppearInSection: >none
#	secondaryAppearInSection: 9_scanMotifOccurence|246
#	input: $allHitBySeqHsh_ref, $cntgMotifHitPlsPathHsh_ref, $motif
#	output: 
#	toCall: &storeMotifHitToHshStorable($cntgMotifHitPlsPathHsh_ref, $allHitBySeqHsh_ref, $motif);
#	calledInLine: 3315
#....................................................................................................................................................#
	my ($cntgMotifHitPlsPathHsh_ref, $allHitBySeqHsh_ref, $motif) = @_;
	
	foreach my $cntg (keys %{$allHitBySeqHsh_ref}) {
		my $cntgMotifHitHsh_ref = retrieve($cntgMotifHitPlsPathHsh_ref->{$cntg});
		foreach my $strnd (keys %{$allHitBySeqHsh_ref->{$cntg}}) {
			foreach my $pos (keys %{$allHitBySeqHsh_ref->{$cntg}{$strnd}}) {
				$cntgMotifHitHsh_ref->{$strnd}{$pos}{$motif} = $allHitBySeqHsh_ref->{$cntg}{$strnd}{$pos};
			}
		}
		store($cntgMotifHitHsh_ref, $cntgMotifHitPlsPathHsh_ref->{$cntg});
	}
	return ();
}
sub zipUnzipCntgCovInPlsPathHsh {
#....................................................................................................................................................#
#	subroutineCategory: storable
#	dependOnSub: reportStatus|3137
#	appearInSub: >none
#	primaryAppearInSection: 14_finishingTasks|320, 5_processInputData|174
#	secondaryAppearInSection: >none
#	input: $cntgCovInPlsPathHsh_ref, $zipUnzip
#	output: none
#	toCall: &zipUnzipCntgCovInPlsPathHsh($zipUnzip, $cntgCovInPlsPathHsh_ref);
#	calledInLine: 195
#....................................................................................................................................................#

	my ($zipUnzip, $cntgCovInPlsPathHsh_ref) = @_;
	
	foreach my $cntg (sort keys %{$cntgCovInPlsPathHsh_ref}) {
		&reportStatus("$zipUnzip cntg ary", 20, "\r");#->3137
		my $cntgCovPlsPath = "$cntgCovInPlsPathHsh_ref->{$cntg}";
		if ($zipUnzip eq 'unzip') {
			system ("gzip -df $cntgCovPlsPath.gz") if (-s "$cntgCovPlsPath.gz");
		} else {
			system ("gzip -f $cntgCovPlsPath") if (-s "$cntgCovPlsPath");
		}
	}
	print "\n";
}

exit;
