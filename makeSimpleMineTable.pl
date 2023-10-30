#!/usr/bin/perl -w
use strict;

print "Create AGR SimpleMine source data table ...\n";

my @spe = ("HUMAN", "MGI", "RGD", "ZFIN", "FB", "WB", "SGD", "XBXL", "XBXT");
my %speMod = ("HUMAN" => "Homo sapiens", 
	       "MGI" => "Mus musculus", 
	       "RGD" => "Rattus norvegicus", 
	       "ZFIN" => "Danio rerio", 
	       "FB" => "Drosophila melanogaster", 
	       "WB" => "Caenorhabditis elegans", 
	       "SGD" => "Saccharomyces cerevisiae",
	       "XBXL" => "Xenopus laevis",
	       "XBXT" => "Xenopus tropicalis");
my ($line, $columns);
my @tmp;
my @stuff;
my %duplicateInfo; #check if the information was aready taken;
my $checkstring;


#------------ Get Gene Identifers and Interaction -----------------
open (INT, "Interaction.csv") || die "can't open Interaction.csv!";
my ($geneid, $modID, $genename, $longname, $s, $ncbiID, $ensemblID, $uniproID, $pantherID, $refseqID, $synonym, $gIntDes, $mIntDes);
my $i = 0; #this index geneList
my $totalGenes;
my @geneList;
my %gName; #public name for each gene id, called Gene Symbol in AGR
my %gSpe; #species for each gene id
my %gModID; #Gene ID in its own Model Organism Database
my %gLName; #Long names for each gene id, called Gene Name in AGR
my %gSYNO; #synonym ids
my %gUniPro; #UniProtKB ID for gene ID
my %gNCBI; #NCBI ID for gene ID
my %gPANTHER; #PANTHER ID for gene ID
my %gENSEM; #ENSEMBL ID for gene ID
my %gRefSeq; #RefSeq ID for gene ID

my %totalSpeGene; #total genes in each species
foreach $s (@spe) {
    $totalSpeGene{$s} = 0;
}
my %gGeneticInt; #genetic interaction for each gene
my %gMolInt; #molecular interaction for each gene

$line = <INT>;
while ($line = <INT>) {
    chomp ($line);
    @tmp = ();
    ($geneid, $modID, $genename, $longname, $s, $ncbiID, $ensemblID, $uniproID, $pantherID, $refseqID, $synonym, $gIntDes, $mIntDes) = split /\t/, $line;
    $geneList[$i] = $geneid; #add this gene to gene list.
    $i++;    
    $gModID{$geneid} = $modID;
    $gName{$geneid} = $genename; #short name, AGR Gene Symbol
    $gLName{$geneid} = $longname; #long name, AGR Gene Name
    $gSpe{$geneid} = $s;
    $gSYNO{$geneid} = $synonym;
    $gNCBI{$geneid} = $ncbiID;
    $gENSEM{$geneid} = $ensemblID;
    $gUniPro{$geneid} = $uniproID;
    $gPANTHER{$geneid} = $pantherID;
    $gRefSeq{$geneid} = $refseqID;
    $gGeneticInt{$geneid} = $gIntDes;
    $gMolInt{$geneid} = $mIntDes;   

    foreach $s (@spe) {
	if ($s eq "$gSpe{$geneid}") {#this gene matches this species
	    $totalSpeGene{$s}++;
	}
    }
}
close (INT);
$totalGenes = $i;
foreach $s (@spe) {
    print "$totalSpeGene{$s} genes found in $speMod{$s}.\n";
}
print "$totalGenes genes found for all species.\n";
#-------- Done getting gene identifiers, genetic and molecular interaction ----


#---------------- Get Gene Description ----------------
open (DES, "/home/wen/agrSimpleMine/src/desFile") || die "can't open desFile!";
my $des;
my %gDes; #Description for each gene id
$i = 0;
while ($line = <DES>) {
        next unless !($line =~ /^\#/);
	chomp ($line);
	if (($line =~ /^HGNC:/)||($line =~ /^MGI:/)||($line =~ /^RGD:/)||($line =~ /^ZFIN:/) ||($line =~ /^FB:/)||($line =~ /^WB:/)||($line =~ /^SGD:/)||($line =~ /Xenbase/)) {
	    ($geneid, $genename) = split /\s+/, $line;
	} elsif ($line ne "") {
	    $des = $line;
	    $gDes{$geneid} = $des;
	    $i++;
	} 
}    
$totalGenes = $i;
close (DES);
print "$totalGenes genes found with description.\n";
#---------------------- Done with Gene Description --------------


#------------ Get Orthology -----------------
my %humanHomo;
my %ratHomo;
my %mouseHomo;
my %fishHomo;
my %flyHomo;
my %wormHomo;
my %yeastHomo;
my %xbxlHomo;
my %xbxtHomo;
my ($h, $hName, $hDes);
my ($date, $version);

open (HOM, "/home/wen/agrSimpleMine/src/orthoFile") || die "can't open orthoFile!";
while ($line = <HOM>) {
    chomp ($line);
    if ($line =~ /Alliance Database Version/) {
	$version = $line;
    } elsif  ($line =~ /Date file generated/) {
	$date = $line;
    }
       next unless (($line =~ /^HGNC:/)||($line =~ /^MGI:/)||($line =~ /^RGD:/)||($line =~ /^ZFIN:/) ||($line =~ /^FB:/)||($line =~ /^WB:/)||($line =~ /^SGD:/)||($line =~ /^Xenbase:/));

	@tmp = split /\t/, $line;
	$geneid = $tmp[0];
	$h = $tmp[4]; #homolog gene id
	$hName = $tmp[5]; #homolog gene name
	$hDes = join ':', $h, $hName;
  
	if ($h =~ /^RGD/) {
	    if ($ratHomo{$geneid}) {
		$ratHomo{$geneid} =  join " \| ",  $ratHomo{$geneid}, $hDes; 
	    } else {
		$ratHomo{$geneid} = $hDes;
	    }
	} elsif  ($h =~ /^MGI/) {
	    if ($mouseHomo{$geneid}) {
		$mouseHomo{$geneid} =  join " \| ",  $mouseHomo{$geneid}, $hDes; 
	    } else {
		$mouseHomo{$geneid} = $hDes;
	    }
	} elsif  ($h =~ /^ZFIN/) {
	    if ($fishHomo{$geneid}) {
		$fishHomo{$geneid} =  join " \| ",  $fishHomo{$geneid}, $hDes; 
	    } else {
		$fishHomo{$geneid} = $hDes;
	    }
	}  elsif  ($h =~ /^FB/) {
	    if ($flyHomo{$geneid}) {
		$flyHomo{$geneid} =  join " \| ",  $flyHomo{$geneid}, $hDes; 
	    } else {
		$flyHomo{$geneid} = $hDes;
	    }
	}  elsif  ($h =~ /^WB/) {
	    if ($wormHomo{$geneid}) {
		$wormHomo{$geneid} =  join " \| ",  $wormHomo{$geneid}, $hDes; 
	    } else {
		$wormHomo{$geneid} = $hDes;
	    }
	}  elsif  ($h =~ /^SGD/) {
	    if ($yeastHomo{$geneid}) {
		$yeastHomo{$geneid} =  join " \| ",  $yeastHomo{$geneid}, $hDes; 
	    } else {
		$yeastHomo{$geneid} = $hDes;
	    }
	} elsif  ($h =~ /^HGNC/) {
	    if ($humanHomo{$geneid}) {
		$humanHomo{$geneid} =  join " \| ",  $humanHomo{$geneid}, $hDes; 
	    } else {
		$humanHomo{$geneid} = $hDes;
	    }
	}  elsif  ($tmp[7] eq "Xenopus laevis") {
	    if ($xbxlHomo{$geneid}) {
		$xbxlHomo{$geneid} =  join " \| ",  $xbxlHomo{$geneid}, $hDes; 
	    } else {
		$xbxlHomo{$geneid} = $hDes;
	    }
	}  elsif  ($tmp[7] eq "Xenopus tropicalis") {
	    if ($xbxtHomo{$geneid}) {
		$xbxtHomo{$geneid} =  join " \| ",  $xbxtHomo{$geneid}, $hDes; 
	    } else {
		$xbxtHomo{$geneid} = $hDes;
	    }
	} 
}    
close (HOM);
#-------- Done getting orthology info ----

#------------ Get Disease Association -----------------
my %gDisease; #diseases associated with each gene
my %existDisease; #diseases that were already annotated to each gene

open (DIS, "/home/wen/agrSimpleMine/src/diseaseFile") || die "can't open diseaseFile!";
my ($diseasename, $doid, $evidence, $diseaseDes);
while ($line = <DIS>) {
       next unless ($line =~ /^NCBITaxon/);
	chomp ($line);
	@tmp = split /\t/, $line;
	$geneid = $tmp[3];
	$evidence = $tmp[5];
	$doid = $tmp[6];
        $diseasename = $tmp[7];
       $checkstring = join "", $geneid, $doid, $evidence;
       if ($existDisease{$checkstring}) {
	   #do nothing since this disease is already annotated to this gene
       } else {
	   $diseaseDes = join ':', $evidence, $doid, $diseasename;
	   if ($gDisease{$geneid}) {
		$gDisease{$geneid} =  join " \| ",  $gDisease{$geneid}, $diseaseDes; 
	   } else {
		$gDisease{$geneid} = $diseaseDes;
	   }
	   $existDisease{$checkstring} = 1;
       }
}    
close (DIS);
#-------- Done getting disease association ----


#------------ Get Variants -----------------
my %gVar; #variants associated with each gene
open (VAR, "/home/wen/agrSimpleMine/src/varFile") || die "can't open varFile!";
my ($varname, $varDes);
my $totalGeneVar = 0;
while ($line = <VAR>) {
    next unless ($line =~ /^NCBITaxon/);    
    chomp ($line);
    @tmp = ();
    @tmp = split /\t/, $line;
    next unless ($tmp[3] ne "");
    next unless ($tmp[9] ne "");
    $varname = $tmp[3];
    $geneid = $tmp[9];
    if ($gVar{$geneid}) {
       $gVar{$geneid} =  join " \| ",  $gVar{$geneid}, $varname; 
    } else {
	$gVar{$geneid} = $varname;
	$totalGeneVar++;
    }
}    
close (VAR);
print "Found Variants info for $totalGeneVar genes.\n";
#-------- Done getting variants ----

#------------ Get Expression -----------------
my %gExprLocation; #expression location for each gene
my %gExprStage; #expression stage for each gene
open (EXP, "/home/wen/agrSimpleMine/src/exprFile") || die "can't open exprFile!";
my ($exprLocation, $exprStage);
while ($line = <EXP>) {
    chomp ($line);
    @tmp = ();
    @tmp = split /\t/, $line;
    $columns = @tmp;
    next unless ($columns > 6);
	
    $geneid = $tmp[2];	

        next unless (($tmp[4] ne "")&&($tmp[4] ne "C. elegans Cell and Anatomy")) ;
	$checkstring = join "", $geneid, $tmp[4];
	next unless !($duplicateInfo{$checkstring});
        $exprLocation = $tmp[4];

	$duplicateInfo{$checkstring} = 1; #this this information as already recorded
	if ($gExprLocation{$geneid}) {
		$gExprLocation{$geneid} =  join " \| ",  $gExprLocation{$geneid}, $exprLocation; 
	} else {
		$gExprLocation{$geneid} = $exprLocation;
	}
	
	next unless (($tmp[5] ne "")&&($tmp[5] ne "Nematoda Life Stage"));	
	$checkstring = join "", $geneid, $tmp[5];
	next unless !($duplicateInfo{$checkstring});
	$exprStage = $tmp[5];
	$duplicateInfo{$checkstring} = 1; #this this information as already recorded
	if ($gExprStage{$geneid}) {
		$gExprStage{$geneid} =  join " \| ",  $gExprStage{$geneid}, $exprStage; 
	} else {
		$gExprStage{$geneid} = $exprStage;
	}	
}    
close (EXP);
#-------- Done getting expression ----



#----- print the source data table -------------

#---- print headers ------------
#my $header = "Gene ID\tGene Symbol\tGene Name\tDescription\tSpecies\tNCBI ID\tENSEMBL ID\tUniProtKB ID\tPANTHER ID\tSynonym\tDisease Association\tExpression Location\tExpression Stage\tVariants\tGenetic Interaction\tMolecular\/Physical Interaction\tHuman Ortholog\tMouse Ortholog\tRat Ortholog\tFish Ortholog\tFly Ortholog\tWorm Ortholog\tYeast Ortholog";
my $header = "Gene ID\tGene Symbol\tGene Name\tDescription\tSpecies\tNCBI ID\tENSEMBL ID\tUniProtKB ID\tPANTHER ID\tRefSeq ID\tSynonym\tDisease Association\tExpression Location\tExpression Stage\tVariants\tGenetic Interaction\tMolecular\/Physical Interaction\tHomo sapiens Ortholog\tMus musculus Ortholog\tRattus norvegicus Ortholog\tDanio rerio Ortholog\tDrosophila melanogaster Ortholog\tCaenorhabditis elegans Ortholog\tSaccharomyces cerevisiae Ortholog\tXenopus laevis Ortholog\tXenopus tropicalis Ortholog";

my $nameHeader = "Gene ID\tGene Symbol\tGene Name\tMOD ID\tNCBI ID\tENSEMBL ID\tUniProtKB ID\tPANTHER ID\tRefSeq ID\tSynonym";

#open (OUT, ">/home/wen/agrSimpleMine/sourceFile/headers") || die "cannot open headers!\n";
#print OUT "$header\n";
#close (OUT);

#---- print version/date------
open (VER, ">/home/wen/agrSimpleMine/sourceFile/version") || die "cannot open version!\n";
print VER "# SimpleMine Results Based on Alliance of Genome Resources\n";
print VER "$version\n";
print VER "$date\n";
close (VER);

open (FBOUT, ">/home/wen/agrSimpleMine/sourceFile/FB/SimpleMineSourceData.csv") || die "cannot open FB/SimpleMineSourceData.csv!\n";
print FBOUT "$header\n";

open (FBNAM, ">/home/wen/agrSimpleMine/sourceFile/FB/GeneName.csv") || die "cannot open FB/GeneName.csv!\n";
print FBNAM "$nameHeader\n";

open (HUMANOUT, ">/home/wen/agrSimpleMine/sourceFile/HUMAN/SimpleMineSourceData.csv") || die "cannot open HUMAN/SimpleMineSourceData.csv!\n";
print HUMANOUT "$header\n";

open (HUMANNAM, ">/home/wen/agrSimpleMine/sourceFile/HUMAN/GeneName.csv") || die "cannot open HUMAN/GeneName.csv!\n";
print HUMANNAM "$nameHeader\n";

open (MGIOUT, ">/home/wen/agrSimpleMine/sourceFile/MGI/SimpleMineSourceData.csv") || die "cannot open MGI/SimpleMineSourceData.csv!\n";
print MGIOUT "$header\n";

open (MGINAM, ">/home/wen/agrSimpleMine/sourceFile/MGI/GeneName.csv") || die "cannot open MGI/GeneName.csv!\n";
print MGINAM "$nameHeader\n";

open (RGDOUT, ">/home/wen/agrSimpleMine/sourceFile/RGD/SimpleMineSourceData.csv") || die "cannot open RGD/SimpleMineSourceData.csv!\n";
print RGDOUT "$header\n";

open (RGDNAM, ">/home/wen/agrSimpleMine/sourceFile/RGD/GeneName.csv") || die "cannot open RGD/GeneName.csv!\n";
print RGDNAM "$nameHeader\n";

open (SGDOUT, ">/home/wen/agrSimpleMine/sourceFile/SGD/SimpleMineSourceData.csv") || die "cannot open SGD/SimpleMineSourceData.csv!\n";
print SGDOUT "$header\n";

open (SGDNAM, ">/home/wen/agrSimpleMine/sourceFile/SGD/GeneName.csv") || die "cannot open SGD/GeneName.csv!\n";
print SGDNAM "$nameHeader\n";

open (WBOUT, ">/home/wen/agrSimpleMine/sourceFile/WB/SimpleMineSourceData.csv") || die "cannot open WB/SimpleMineSourceData.csv!\n";
print WBOUT "$header\n";

open (WBNAM, ">/home/wen/agrSimpleMine/sourceFile/WB/GeneName.csv") || die "cannot open WB/GeneName.csv!\n";
print WBNAM "$nameHeader\n";

open (ZFINOUT, ">/home/wen/agrSimpleMine/sourceFile/ZFIN/SimpleMineSourceData.csv") || die "cannot open ZFIN/SimpleMineSourceData.csv!\n";
print ZFINOUT "$header\n";

open (ZFINNAM, ">/home/wen/agrSimpleMine/sourceFile/ZFIN/GeneName.csv") || die "cannot open ZFIN/GeneName.csv!\n";
print ZFINNAM "$nameHeader\n";


open (XBXLOUT, ">/home/wen/agrSimpleMine/sourceFile/XBXL/SimpleMineSourceData.csv") || die "cannot open XBXL/SimpleMineSourceData.csv!\n";
print XBXLOUT "$header\n";

open (XBXLNAM, ">/home/wen/agrSimpleMine/sourceFile/XBXL/GeneName.csv") || die "cannot open XBXL/GeneName.csv!\n";
print XBXLNAM "$nameHeader\n";


open (XBXTOUT, ">/home/wen/agrSimpleMine/sourceFile/XBXT/SimpleMineSourceData.csv") || die "cannot open XBXT/SimpleMineSourceData.csv!\n";
print XBXTOUT "$header\n";

open (XBXTNAM, ">/home/wen/agrSimpleMine/sourceFile/XBXT/GeneName.csv") || die "cannot open XBXT/GeneName.csv!\n";
print XBXTNAM "$nameHeader\n";


open (MIXOUT, ">/home/wen/agrSimpleMine/sourceFile/MIX/SimpleMineSourceData.csv") || die "cannot open ZFIN/SimpleMineSourceData.csv!\n";
print MIXOUT "$header\n";

open (MIXNAM, ">/home/wen/agrSimpleMine/sourceFile/MIX/GeneName.csv") || die "cannot open ZFIN/GeneName.csv!\n";
print MIXNAM "$nameHeader\n";

my ($humanOrtho, $ratOrtho, $mouseOrtho, $fishOrtho, $flyOrtho, $wormOrtho, $yeastOrtho, $xbxlOrtho, $xbxtOrtho, $exprLocDes, $exprStaDes);

foreach $geneid (@geneList) {
    next unless ($gName{$geneid});

    if ($gName{$geneid}) {
	$genename = $gName{$geneid};
    } else {
	$genename = "N.A.";
    }
    if ($gLName{$geneid}) {
	$longname = $gLName{$geneid};
    } else {
	$longname = "N.A.";
    }    
    if ($gModID{$geneid}) {
	$modID = $gModID{$geneid};
    } else {
	$modID = "N.A.";
    }
    if  ($gUniPro{$geneid}) {
	$uniproID = $gUniPro{$geneid};
    } else {
	$uniproID = "";
    }
    if  ($gNCBI{$geneid}) {
	$ncbiID = $gNCBI{$geneid};
    } else {
	$ncbiID = "";
    }
    if  ($gPANTHER{$geneid}) {
	$pantherID = $gPANTHER{$geneid};
    } else {
	$pantherID = "";
    }
    if  ($gRefSeq{$geneid}) {
	$refseqID = $gRefSeq{$geneid};
    } else {
	$refseqID = "";
    }
    if  ($gENSEM{$geneid}) {
	$ensemblID = $gENSEM{$geneid};
    } else {
	$ensemblID = "";
    }
    if ($gSYNO{$geneid}) {
       $synonym = $gSYNO{$geneid};
    } else {
       $synonym = "";
    }  
    if ($gDes{$geneid}) {
	$des = $gDes{$geneid};
    } else {
	$des = "N.A.";
    }
    if ($gSpe{$geneid}) {
	$s = $gSpe{$geneid};
    } else {
	$s = "N.A.";
	print "ERROR! No species information for $geneid!\n";
    }
    if ($humanHomo{$geneid}) {
	$humanOrtho = $humanHomo{$geneid};
    } else {
	$humanOrtho = "N.A.";
    }
    if ($mouseHomo{$geneid}) {
	$mouseOrtho = $mouseHomo{$geneid};
    } else {
	$mouseOrtho = "N.A.";
    }
    if ($ratHomo{$geneid}) {
	$ratOrtho = $ratHomo{$geneid};
    } else {
	$ratOrtho = "N.A.";
    }
    if ($fishHomo{$geneid}) {
	$fishOrtho = $fishHomo{$geneid};
    } else {
	$fishOrtho = "N.A.";
    }
    if ($flyHomo{$geneid}) {
	$flyOrtho = $flyHomo{$geneid};
    } else {
	$flyOrtho = "N.A.";
    }
    if ($wormHomo{$geneid}) {
	$wormOrtho = $wormHomo{$geneid};
    } else {
	$wormOrtho = "N.A.";
    }
    if ($yeastHomo{$geneid}) {
	$yeastOrtho = $yeastHomo{$geneid};
    } else {
	$yeastOrtho = "N.A.";
    }
    if ($xbxlHomo{$geneid}) {
	$xbxlOrtho = $xbxlHomo{$geneid};
    } else {
	$xbxlOrtho = "N.A.";
    }
    if ($xbxtHomo{$geneid}) {
	$xbxtOrtho = $xbxtHomo{$geneid};
    } else {
	$xbxtOrtho = "N.A.";
    }    
    if ($gDisease{$geneid}) {
	$diseaseDes = $gDisease{$geneid};
    } else {
        $diseaseDes = "N.A.";
    }
    if ($gExprLocation{$geneid}) {
	$exprLocDes = $gExprLocation{$geneid};	
    } else {
        $exprLocDes = "N.A.";
    }
    if ($gExprStage{$geneid}) {
	$exprStaDes = $gExprStage{$geneid};
    } else {
        $exprStaDes = "N.A.";
    }
     if ($gVar{$geneid}) {
	$varDes = $gVar{$geneid};
    } else {
        $varDes = "N.A.";
    }    
    if ($gGeneticInt{$geneid}) {
    	$gIntDes = $gGeneticInt{$geneid};
    } else {
    	$gIntDes = "N.A.";
    }
    if ($gMolInt{$geneid}) {
    	$mIntDes = $gMolInt{$geneid};
    } else {
    	$mIntDes = "N.A.";
    }
   
    #print species specific source data
#my $header = "AGR Gene ID\tAGR Gene Name\tDescription\tSpecies\t\tNCBI ID\tENSEMBL ID\tUniProtKB ID\tPANTHER ID\tDisease Association\tExpression Location\tExpression Stage\tVariants\tGenetic Interaction\tMolecular\/Physical Interaction\tHuman Ortholog\tMouse Ortholog\tRat Ortholog\tFish Ortholog\tFly Ortholog\tWorm Ortholog\tYeast Ortholog";   
    my $dataline = join "\t", $geneid, $genename, $longname, $des, $speMod{$s}, $ncbiID, $ensemblID, $uniproID, $pantherID, $refseqID, $synonym, $diseaseDes, $exprLocDes, $exprStaDes, $varDes, $gIntDes, $mIntDes, $humanOrtho, $mouseOrtho, $ratOrtho, $fishOrtho, $flyOrtho, $wormOrtho, $yeastOrtho, $xbxlOrtho, $xbxtOrtho;

#my $nameHeader = "AGR Gene ID\tAGR Gene Name\tMOD ID\tNCBI ID\tENSEMBL ID\tUniProtKB ID\tPANTHER ID";
    my $nameline = join "\t", $geneid, $genename, $longname, $modID, $ncbiID, $ensemblID, $uniproID, $pantherID, $refseqID, $synonym;
    
    if ($s eq "HUMAN") {
	print HUMANOUT "$dataline\n";
        print HUMANNAM "$nameline\n";
    } elsif ($s eq "FB") {
	print FBOUT "$dataline\n";
        print FBNAM "$nameline\n";
    } elsif ($s eq "MGI") {
	print MGIOUT "$dataline\n";
        print MGINAM "$nameline\n";
    } elsif ($s eq "RGD") {
	print RGDOUT "$dataline\n";
        print RGDNAM "$nameline\n";
    } elsif ($s eq "SGD") {
	print SGDOUT "$dataline\n";
        print SGDNAM "$nameline\n";
    } elsif ($s eq "WB") {
	print WBOUT "$dataline\n";
        print WBNAM "$nameline\n";
    } elsif ($s eq "ZFIN") {
	print ZFINOUT "$dataline\n";
        print ZFINNAM "$nameline\n";
    } elsif ($s eq "XBXL") {
	print XBXLOUT "$dataline\n";
        print XBXLNAM "$nameline\n";
    } elsif ($s eq "XBXT") {
	print XBXTOUT "$dataline\n";
        print XBXTNAM "$nameline\n";
    } else {
	print "ERROR! Cannot find species for $s!\n";
    }

    print MIXOUT "$dataline\n";
    print MIXNAM "$nameline\n";

    
}

close (FBOUT);
close (FBNAM);
close (HUMANOUT);
close (HUMANNAM);
close (MGIOUT);
close (MGINAM);
close (RGDOUT);
close (RGDNAM);
close (SGDOUT);
close (SGDNAM);
close (WBOUT);
close (WBNAM);
close (ZFINOUT);
close (ZFINNAM);
close (XBXLOUT);
close (XBXLNAM);
close (XBXTOUT);
close (XBXTNAM);
close (MIXOUT);
close (MIXNAM);
