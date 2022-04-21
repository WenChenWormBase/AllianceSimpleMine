#!/usr/bin/perl -w
use strict;

print "Create AGR SimpleMine source data table ...\n";

#my @spe = ("HUMAN", "MGI", "RGD", "ZFIN", "FB", "WB", "SGD");
my %speMod = ("HGNC" => "Homo sapiens", 
	       "MGI" => "Mus musculus", 
	       "RGD" => "Rattus norvegicus", 
	       "ZFIN" => "Danio rerio", 
	       "FB" => "Drosophila melanogaster", 
	       "WB" => "Caenorhabditis elegans", 
	       "SGD" => "Saccharomyces cerevisiae");
my $s; #species name
my $modID; #gene id in its own Model Organism Database

my ($line, $columns);
my @tmp;
my @stuff;
my $i; #this index geneList
my $totalGenes;
my @geneList;
my ($geneid, $genename, $des);
my %duplicateInfo; #check if the information was aready taken;
my $checkstring;


#---------------- Get Gene Description -----------------------------------
open (DES, "desFile") || die "can't open desFile!";
$i = 0;

my %gName; #public name for each gene id
my %gDes; #Description for each gene id
my %gSpe; #species for each gene id
my %gModID; #Gene ID in its own Model Organism Database
my $totalHumanGene = 0;
my $totalMGIGene = 0;
my $totalRGDGene = 0;
my $totalZFINGene = 0;
my $totalFBGene = 0;
my $totalWBGene = 0;
my $totalSGDGene = 0;

while ($line = <DES>) {
        next unless !($line =~ /^\#/);
	chomp ($line);
	if (($line =~ /^HGNC:/)||($line =~ /^MGI:/)||($line =~ /^RGD:/)||($line =~ /^ZFIN:/) ||($line =~ /^FB:/)||($line =~ /^WB:/)||($line =~ /^SGD:/)) {
	    ($geneid, $genename) = split /\s+/, $line;
	    $gName{$geneid} = $genename;
	    $geneList[$i] = $geneid;
	    $i++;
	    ($s, $modID) = split ":", $geneid;
	    $gSpe{$geneid} = $s;
	    $gModID{$geneid} = $modID;
	    if ($s eq "HGNC") { 
		$totalHumanGene++;
	    } elsif ($s eq "MGI") { 
		$totalMGIGene++; 
	    } elsif ($s eq "RGD") { 
		$totalRGDGene++;
	    } elsif ($s eq "ZFIN") { 
		$totalZFINGene++; 
	    } elsif ($s eq "FB") { 
		$totalFBGene++; 
	    } elsif ($s eq "WB") { 
		$totalWBGene++; 
	    } elsif ($s eq "SGD") { 
		$totalSGDGene++; 
	    }		
	} elsif ($line ne "") {
	    $gDes{$geneid} = $line;
	} 
}    
$totalGenes = $i;
close (DES);
print "$totalGenes genes found with description\n";
#---------------------- Done with Gene Description --------------



#---------------- Get Genomic Resource IDs -----------------------------------
open (XREF, "xrefFile") || die "can't open xrefFile!";
$i = 0;

my %gUniPro; #UniProtKB IDs
my %gNCBI; #NCBI IDs
my %gPANTHER; #PANTHER IDs
my %gENSEM; #ENSEMBL IDs
my ($resource, $rID);


while ($line = <XREF>) {
    chomp ($line);
    @tmp = ();
    @tmp = split /\t/, $line;
    $columns = @tmp;
    next unless ($columns > 2);
    $geneid = $tmp[0];
	if ($tmp[1] =~ /^UniProtKB:/) {
	    ($resource, $rID) = split ":", $tmp[1];
	    if ($gUniPro{$geneid}) {
		$gUniPro{$geneid} = join " | ",  $gUniPro{$geneid}, $rID;
	    } else {
		$gUniPro{$geneid} = $rID;
	    }
	} elsif ($tmp[1] =~ /^NCBI_Gene:/) {
	    ($resource, $rID) = split ":", $tmp[1];
	    $gNCBI{$geneid} = $rID;
	} elsif ($tmp[1] =~ /ENSEMBL:/) {
	    ($resource, $rID) = split ":", $tmp[1];
	    $gENSEM{$geneid} = $rID;
	} elsif ($tmp[1] =~ /PANTHER:/) {
	    ($resource, $rID) = split ":", $tmp[1];
	    $gPANTHER{$geneid} = $rID;
	} 
}    
close (XREF);
#---------------------- Done with Gene Description --------------


#--------------Get Synonyms for Gene Names --------------
open (SYNO, "synonymFile") || die "can't open synonymFile!";
my %gSYNO; #synonym ids
my $totalSynGene = 0;
my $totalColumns = 0;
my $synonym;
my @findSyn;

$line = <SYNO>;

while ($line = <SYNO>) {
    #next unless !($line =~ /^\#/);
    chomp ($line);
    @tmp = ();
    if ($line =~ /"/) {
       @findSyn = ();
       @findSyn = split '"', $line;
       $synonym = $findSyn[1]; 
    } else {
       $synonym = "N.A."; 
    }
    @tmp = split ",", $line; 
    $geneid = $tmp[0];	
    if ($synonym eq "N.A.") {
       $synonym = $tmp[1];
    } 
    if ($gSYNO{$geneid}) {
       $gSYNO{$geneid} = join " | ", $gSYNO{$geneid}, $synonym;
    } else {
       $gSYNO{$geneid} = $synonym;
       $totalSynGene++;
    }
}
close (SYNO);
print "$totalSynGene genes have synonyms.\n";
#-------------Done with Gene Name Synonyms -------------


#-----------Get Name into Synonym list -------------

open (LNAM, "nameFile") || die "can't open nameFile!";
my $longname;
my %gLName;
    
$line = <LNAM>;

while ($line = <LNAM>) {
    #next unless !($line =~ /^\#/);
    chomp ($line);
    @tmp = ();
    if ($line =~ /"/) {
       @findSyn = ();
       @findSyn = split '"', $line;
       $synonym = $findSyn[1]; 
    } else {
       $synonym = "N.A."; 
    }
    @tmp = split /\t/, $line; 
    $geneid = $tmp[0];
    $longname = $tmp[2];
    $gLName{$geneid} = $longname;
	
    #if ($synonym eq "N.A.") {
    #   $synonym = $tmp[2];
    #} 
    #if ($gSYNO{$geneid}) {
    #   $gSYNO{$geneid} = join " | ", $gSYNO{$geneid}, $synonym;
    #} else {
    #   $gSYNO{$geneid} = $synonym;
    #   $totalSynGene++;
    #}
}
close (LNAM);
#print "$totalSynGene genes have synonyms after including their official names.\n";
#----------Done getting Name --------------



#------------ Get Orthology -----------------
my %humanHomo;
my %ratHomo;
my %mouseHomo;
my %fishHomo;
my %flyHomo;
my %wormHomo;
my %yeastHomo;
my ($h, $hName, $hDes);
my ($date, $version);

open (HOM, "orthoFile") || die "can't open orthoFile!";
while ($line = <HOM>) {
    chomp ($line);
    if ($line =~ /Alliance Database Version/) {
	$version = $line;
    } elsif  ($line =~ /Date file generated/) {
	$date = $line;
    }

       next unless (($line =~ /^HGNC:/)||($line =~ /^MGI:/)||($line =~ /^RGD:/)||($line =~ /^ZFIN:/) ||($line =~ /^FB:/)||($line =~ /^WB:/)||($line =~ /^SGD:/));

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
	}
}    
close (HOM);
#-------- Done getting orthology info ----

#------------ Get Disease Association -----------------
my %gDisease; #diseases associated with each gene
my %existDisease; #diseases that were already annotated to each gene

open (DIS, "diseaseFile") || die "can't open diseaseFile!";
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
open (VAR, "varFile") || die "can't open diseaseFile!";
my ($varname, $so, $varDes);

my $totalGeneVar = 0;
while ($line = <VAR>) {
    next unless ($line =~ /soTerm/);
    next unless ($line =~ /allele_of_gene_ids/);
    next unless ($line =~ /allele_symbols/);

    chomp ($line);
    @stuff = ();
    ($stuff[0], $stuff[1]) = split "allele_ids", $line;
    @tmp = split '"', $stuff[1];
    $so = $tmp[7];
    $geneid = $tmp[9];
    $varname = $tmp[5];
    $varDes = join ':', $so, $varname;
    if ($gVar{$geneid}) {
       $gVar{$geneid} =  join " \| ",  $gVar{$geneid}, $varDes; 
    } else {
	$gVar{$geneid} = $varDes;
	$totalGeneVar++;
    }
}    
close (VAR);
print "Found Variants info for $totalGeneVar genes.\n";
#-------- Done getting variants ----

#------------ Get Expression -----------------
my %gExprLocation; #expression location for each gene
my %gExprStage; #expression stage for each gene
open (EXP, "exprFile") || die "can't open exprFile!";
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


#------------ Get Interaction -----------------
open (INT, "Interaction.csv") || die "can't open Interaction.csv!";
my %gGeneticInt; #genetic interaction for each gene
my %gMolInt; #molecular interaction for each gene
my ($gIntDes, $mIntDes);
$line = <INT>;
while ($line = <INT>) {
    chomp ($line);
    @tmp = ();
    ($geneid, $gIntDes, $mIntDes) = split /\t/, $line;
    $gGeneticInt{$geneid} = $gIntDes;
    $gMolInt{$geneid} = $mIntDes;    
}
close (INT);
#-------- Done getting both genetic and molecular interaction ----


#----- print the source data table -------------
#open (NAM, ">GeneName.csv") || die "cannot open GeneName.csv!\n";
#print NAM "AGR Gene ID\tAGR Gene Name\tSpecies\t\n";
#open (OUT, ">SimpleMineSourceData.csv") || die "cannot open SimpleMineSourceData.csv!\n";
#print OUT "AGR Gene ID\tAGR Gene Name\tDescription\tSpecies\tHuman Ortholog\tMouse Ortholog\tRat Ortholog\tFish Ortholog\tFly Ortholog\tWorm Ortholog\tYeast Ortholog\tDisease Association\tExpression Location\tExpression Stage\tVariants\tGenetic Interaction\tMolecular\/Physical Interaction\n";

#---- print headers ------------
my $header = "Gene ID\tGene Symbol\tGene Name\tDescription\tSpecies\tNCBI ID\tENSEMBL ID\tUniProtKB ID\tPANTHER ID\tSynonym\tDisease Association\tExpression Location\tExpression Stage\tVariants\tGenetic Interaction\tMolecular\/Physical Interaction\tHuman Ortholog\tMouse Ortholog\tRat Ortholog\tFish Ortholog\tFly Ortholog\tWorm Ortholog\tYeast Ortholog";
#my $nameHeader = "AGR Gene ID\tAGR Gene Name\tMOD ID\tNCBI ID\tENSEMBL ID\tUniProtKB ID\tPANTHER ID";
my $nameHeader = "Gene ID\tGene Symbol\tGene Name\tMOD ID\tNCBI ID\tENSEMBL ID\tUniProtKB ID\tPANTHER ID\tSynonym";

open (OUT, ">/home/wen/agrSimpleMine/sourceFile/headers") || die "cannot open headers!\n";
print OUT "$header\n";
close (OUT);

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


open (MIXOUT, ">/home/wen/agrSimpleMine/sourceFile/MIX/SimpleMineSourceData.csv") || die "cannot open ZFIN/SimpleMineSourceData.csv!\n";
print MIXOUT "$header\n";

open (MIXNAM, ">/home/wen/agrSimpleMine/sourceFile/MIX/GeneName.csv") || die "cannot open ZFIN/GeneName.csv!\n";
print MIXNAM "$nameHeader\n";

my ($uniproID, $ncbiID, $pantherID, $ensemblID, $humanOrtho, $ratOrtho, $mouseOrtho, $fishOrtho, $flyOrtho, $wormOrtho, $yeastOrtho, $exprLocDes, $exprStaDes);

foreach $geneid (@geneList) {
    next unless ($gName{$geneid});


#my %gUniPro; #UniProtKB IDs
#my %gNCBI; #NCBI IDs
#my %gPANTHER; #PANTHER IDs
#my %gENSEM; #ENSEMBL IDs
    
    #--- fill in the blank fields ----------
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
	$modID = "";
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
    my $dataline = join "\t", $geneid, $genename, $longname, $des, $speMod{$s}, $ncbiID, $ensemblID, $uniproID, $pantherID, $synonym, $diseaseDes, $exprLocDes, $exprStaDes, $varDes, $gIntDes, $mIntDes, $humanOrtho, $mouseOrtho, $ratOrtho, $fishOrtho, $flyOrtho, $wormOrtho, $yeastOrtho;

#my $nameHeader = "AGR Gene ID\tAGR Gene Name\tMOD ID\tNCBI ID\tENSEMBL ID\tUniProtKB ID\tPANTHER ID";
    my $nameline = join "\t", $geneid, $genename, $longname, $modID, $ncbiID, $ensemblID, $uniproID, $pantherID, $synonym;
    
    if ($s eq "HGNC") {
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
close (MIXOUT);
close (MIXNAM);

print "Total Human Genes: $totalHumanGene\n";
print "Total Mouse Genes: $totalMGIGene\n";
print "Total Rat Genes: $totalRGDGene\n";
print "Total Fish Genes: $totalZFINGene\n";
print "Total Fly Genes: $totalFBGene\n";
print "Total Worm Genes: $totalWBGene\n";
print "Total Yeast Genes: $totalSGDGene\n";
