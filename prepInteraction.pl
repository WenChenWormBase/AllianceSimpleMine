#!/usr/bin/perl -w
use strict;

print "Create AGR SimpleMine source data table ...\n";

my @spe = ("HUMAN", "MGI", "RGD", "ZFIN", "FB", "WB", "SGD");
my %taxidSpe = ("taxid:9606" => "HUMAN", 
	       "taxid:10090" => "MGI", 
	       "taxid:10116" => "RGD", 
	       "taxid:7955" => "ZFIN", 
	       "taxid:7227" => "FB", 
	       "taxid:6239" => "WB", 
	       "taxid:559292" => "SGD");
my $s; #species name

my ($line, $columns);
my @tmp;
my $i; #this index geneList
my $totalGenes;
my @geneList;
my ($geneid, $genename, $des, $taxid, $entrez);
my %duplicateInfo; #check if the information was aready taken;
my $checkstring;

#---- Screen interactor files, get entrez ID - public name table ---
my %entrezName; #given NCBI entrez ID, find the gene name
my %nameEntrez; #given the gene name, find the NCBI entrez ID
my @stuff;
my $totalEntrezGenes = 0;
my $speGeneName; #a unique name combining species taxid and gene name
    
open (MINT, "molIntFile") || die "can't open molIntFile!";
while ($line = <MINT>) {
    chomp ($line);
    next unless ($line =~ /^entrez/);    
    @tmp = ();
    @tmp = split /\t/, $line;

    @stuff= ();
    ($stuff[0], $entrez) = split ":", $tmp[0];
    ($stuff[0], $stuff[1], $genename) = split ":", $tmp[2];
    $taxid = $tmp[9];
    if ($taxidSpe{$taxid}) {
	next unless !($entrezName{$entrez});
	$entrezName{$entrez} = $genename;
	$speGeneName = join ":", $taxidSpe{$taxid}, $genename;
	$nameEntrez{$speGeneName} = $entrez;
	$totalEntrezGenes++;
    } else {
	print "ERROR! Cannot find species info for $taxid!\n"
    }

    @stuff= ();
    ($stuff[0], $entrez) = split ":", $tmp[1];
    ($stuff[0], $stuff[1], $genename) = split ":", $tmp[3];
    $taxid = $tmp[9];
    if ($taxidSpe{$taxid}) {
	next unless !($entrezName{$entrez});
	$entrezName{$entrez} = $genename;
	$speGeneName = join ":", $taxidSpe{$taxid}, $genename;
	$nameEntrez{$speGeneName} = $entrez;
	$totalEntrezGenes++;
    } else {
	print "ERROR! Cannot find species info for $taxid!\n"
    }
       
}
close (MINT);
print "Found a total of $totalEntrezGenes Entrez gene names from molIntFile.\n";

     
open (GINT, "geneticIntFile") || die "can't open geneticIntFile!";
while ($line = <GINT>) {
    chomp ($line);
    next unless ($line =~ /^entrez/);    
    @tmp = ();
    @tmp = split /\t/, $line;

    @stuff= ();
    ($stuff[0], $entrez) = split ":", $tmp[0];
    ($stuff[0], $stuff[1], $genename) = split ":", $tmp[2];
    $taxid = $tmp[9];
    if ($taxidSpe{$taxid}) {
	next unless !($entrezName{$entrez});
	$entrezName{$entrez} = $genename;
	$speGeneName = join ":", $taxidSpe{$taxid}, $genename;
	$nameEntrez{$speGeneName} = $entrez;
	$totalEntrezGenes++;
    } else {
	print "ERROR! Cannot find species info for $taxid!\n"
    }

    @stuff= ();
    ($stuff[0], $entrez) = split ":", $tmp[1];
    ($stuff[0], $stuff[1], $genename) = split ":", $tmp[3];
    $taxid = $tmp[9];
    if ($taxidSpe{$taxid}) {
	next unless !($entrezName{$entrez});
	$entrezName{$entrez} = $genename;
	$speGeneName = join ":", $taxidSpe{$taxid}, $genename;
	$nameEntrez{$speGeneName} = $entrez;
	$totalEntrezGenes++;
    } else {
	print "ERROR! Cannot find species info for $taxid!\n"
    }
       
}
close (GINT);
print "Found a total of $totalEntrezGenes Entrez gene names from molIntFile and geneticIntFile.\n";
   
#---------------- Get Gene Description -----------------------------------
open (DES, "desFile") || die "can't open desFile!";
$i = 0;

my %gName; #public name for each gene id
my %gDes; #Description for each gene id
my %gSpe; #species for each gene id
my %speGeneNameID; #given the species and gene name $speGeneName, find the AGR gene id
my $totalHumanGene = 0;
my $totalMGIGene = 0;
my $totalRGDGene = 0;
my $totalZFINGene = 0;
my $totalFBGene = 0;
my $totalWBGene = 0;
my $totalSGDGene = 0;

while ($line = <DES>) {
	chomp ($line);
	if (($line =~ /^HGNC:/)||($line =~ /^MGI:/)||($line =~ /^RGD:/)||($line =~ /^ZFIN:/) ||($line =~ /^FB:/)||($line =~ /^WB:/)||($line =~ /^SGD:/)) {
	    ($geneid, $genename) = split /\s+/, $line;
	    $gName{$geneid} = $genename;
	    $geneList[$i] = $geneid;
	    $i++;
	    if ($line =~ /^HGNC:/) { 
		$totalHumanGene++;
                $gSpe{$geneid} = "HUMAN";
	    } elsif ($line =~ /^MGI:/) { 
		$totalMGIGene++; 
                $gSpe{$geneid} = "MGI";
		
	    } elsif ($line =~ /^RGD:/) { 
		$totalRGDGene++;
                $gSpe{$geneid} = "RGD";		
	    } elsif ($line =~ /^ZFIN:/) { 
		$totalZFINGene++; 
                $gSpe{$geneid} = "ZFIN";
	    } elsif ($line =~ /^FB:/) { 
		$totalFBGene++; 
                $gSpe{$geneid} = "FB";
	    } elsif ($line =~ /^WB:/) { 
		$totalWBGene++; 
                $gSpe{$geneid} = "WB";
	    } elsif ($line =~ /^SGD:/) { 
		$totalSGDGene++; 
                $gSpe{$geneid} = "SGD";
	    }		
	    $speGeneName = join ":", $gSpe{$geneid}, $genename;
	    if ($speGeneNameID{$speGeneName}) {
		if ($speGeneNameID{$speGeneName} ne $geneid) {
		    print "Warning! $geneid $speGeneName already has a different geneid $speGeneNameID{$speGeneName}!\n";
		}
	    } else {
		$speGeneNameID{$speGeneName} = $geneid;
	    }
	} elsif ($line ne "") {
	    $gDes{$geneid} = $line;
	} 
}    
$totalGenes = $i;
close (DES);
print "$totalGenes genes found with description\n";
#---------------------- Done with Gene Description --------------


#------------ Get Genetic Interaction -----------------
my %gGeneticInt; #genetic interaction for each gene
open (GINT, "geneticIntFile") || die "can't open geneticIntFile!";
#open (GINT, "humanGeneticIntFile") || die "can't open humanGeneticIntFile!";
#open (MINT, "molIntFile") || die "can't open molIntFile!";
my ($gInteractor, $intA, $intB, $tempstring);
my @splitName;

my $totalGeneGInt = 0; 
while ($line = <GINT>) {
    chomp ($line);
    next unless (($line =~ /^wormbase/) || ($line =~ /^flybase/)||($line =~ /^entrez/));

    @tmp = ();
    @tmp = split /\t/, $line;

    #convert the names into AGR ID
    $intA = "";
    $intB = "";
    if (($line =~ /^wormbase/) || ($line =~ /^flybase/)) {
       @splitName = ();
       @splitName = split ":", $tmp[0];
       if (($splitName[0]) =~ /flybase/) {
		$splitName[0] = "FB";
       } elsif  (($splitName[0]) =~ /wormbase/) {
		$splitName[0] = "WB";
       }	
       $intA = join ":", @splitName;

       @splitName = ();
       @splitName = split ":", $tmp[1];
       if (($splitName[0]) =~ /flybase/) {
		$splitName[0] = "FB";
       } elsif  (($splitName[0]) =~ /wormbase/) {
		$splitName[0] = "WB";
       }	
       $intB = join ":", @splitName;
    } elsif ($line =~ /^entrez/) {
	#biogrid data
	#next unless ($line =~ //);
	# find $intA
	@stuff= ();
	@stuff = split /\|/, $tmp[2];
	$tempstring = $stuff[1];
	@stuff= ();
	@stuff = split ":", $tempstring;	
	$genename = $stuff[1];
	$taxid = $tmp[9];
	if ($taxidSpe{$taxid}) {
	    #next unless !($entrezName{$entrez});
	    #$entrezName{$entrez} = $genename;
	    $speGeneName = join ":", $taxidSpe{$taxid}, $genename;
	    if ($speGeneNameID{$speGeneName}) {
		$geneid = $speGeneNameID{$speGeneName};
		$intA = $geneid;
	    } else {
		#print "ERROR! Cannot find gene id for $speGeneName!\n";
	    }
	} else {
	    print "ERROR! Cannot find species info for $taxid!\n"
	}

	# find $intB
	@stuff= ();
	@stuff = split /\|/, $tmp[3];
	$tempstring = $stuff[1];
	@stuff= ();
	@stuff = split ":", $tempstring;	
	$genename = $stuff[1];
	$taxid = $tmp[10];
	if ($taxidSpe{$taxid}) {
	    #next unless !($entrezName{$entrez});
	    #$entrezName{$entrez} = $genename;
	    $speGeneName = join ":", $taxidSpe{$taxid}, $genename;
	    if ($speGeneNameID{$speGeneName}) {
		$geneid = $speGeneNameID{$speGeneName};
		$intB = $geneid;
	    } else {
		#print "ERROR! Cannot find gene id for $speGeneName!\n";
	    }
	} else {
	    print "ERROR! Cannot find species info for $taxid!\n"
	}	
    }
    #----- done converting -------

    #------ start to process interaction data -----
    
    next unless (($intA ne "")&&($intB ne ""));
    $geneid = $intA;
    $gInteractor = $intB;
    $checkstring = join "gInt", $geneid, $gInteractor;
    next unless !($duplicateInfo{$checkstring}); 
    #if  ($geneid =~ /^FB/) {
    #	print "$checkstring\n";
    #}
    next unless $gName{$gInteractor}; #must be able to find interactor name
    $tempstring = join ":", $gInteractor, $gName{$gInteractor};
    if ($gGeneticInt{$geneid}) {
       $gGeneticInt{$geneid} =  join " \| ",  $gGeneticInt{$geneid}, $tempstring; 
    } else {
	$gGeneticInt{$geneid} = $tempstring;
	$totalGeneGInt++;
    }    
    $duplicateInfo{$checkstring} = 1; #record this pair of genetic interaction


    $geneid = $intB;
    $gInteractor = $intA;  
    $checkstring = join "gInt", $geneid, $gInteractor;
    next unless !($duplicateInfo{$checkstring});
    next unless $gName{$gInteractor}; #must be able to find interactor name
    $tempstring = join ":", $gInteractor, $gName{$gInteractor};
    if ($gGeneticInt{$geneid}) {
       $gGeneticInt{$geneid} =  join " \| ",  $gGeneticInt{$geneid}, $tempstring; 
    } else {
	$gGeneticInt{$geneid} = $tempstring;
	$totalGeneGInt++;
    }    
    $duplicateInfo{$checkstring} = 1; #record this pair of genetic interaction
    
}    
close (GINT);
print "Found Genetic interaction info for $totalGeneGInt genes.\n";
#close (MINT);
#-------- Done getting genetic interaction ----


#------------ Get Molecular Interaction -----------------
my %gMolInt; #molecular interaction for each gene
open (MINT, "molIntFile") || die "can't open molIntFile!";
my $mInteractor;

my $totalGeneMInt = 0;

while ($line = <MINT>) {
    chomp ($line);
    next unless (($line =~ /^wormbase/) || ($line =~ /^flybase/)||($line =~ /^entrez/));

    @tmp = ();
    @tmp = split /\t/, $line;

    #convert the names into AGR ID
    $intA = "";
    $intB = "";
    if (($line =~ /^wormbase/) || ($line =~ /^flybase/)) {
       @splitName = ();
       @splitName = split ":", $tmp[0];
       if (($splitName[0]) =~ /flybase/) {
		$splitName[0] = "FB";
       } elsif  (($splitName[0]) =~ /wormbase/) {
		$splitName[0] = "WB";
       }	
       $intA = join ":", @splitName;

       @splitName = ();
       @splitName = split ":", $tmp[1];
       if (($splitName[0]) =~ /flybase/) {
		$splitName[0] = "FB";
       } elsif  (($splitName[0]) =~ /wormbase/) {
		$splitName[0] = "WB";
       }	
       $intB = join ":", @splitName;
    } elsif ($line =~ /^entrez/) {
	#biogrid data
	#next unless ($line =~ //);
	# find $intA
	@stuff= ();
	@stuff = split /\|/, $tmp[2];
	$tempstring = $stuff[1];
	@stuff= ();
	@stuff = split ":", $tempstring;	
	$genename = $stuff[1];
	$taxid = $tmp[9];
	if ($taxidSpe{$taxid}) {
	    #next unless !($entrezName{$entrez});
	    #$entrezName{$entrez} = $genename;
	    $speGeneName = join ":", $taxidSpe{$taxid}, $genename;
	    if ($speGeneNameID{$speGeneName}) {
		$geneid = $speGeneNameID{$speGeneName};
		$intA = $geneid;
	    } else {
		#print "ERROR! Cannot find gene id for $speGeneName!\n";
	    }
	} else {
	    print "ERROR! Cannot find species info for $taxid!\n"
	}

	# find $intB
	@stuff= ();
	@stuff = split /\|/, $tmp[3];
	$tempstring = $stuff[1];
	@stuff= ();
	@stuff = split ":", $tempstring;	
	$genename = $stuff[1];
	$taxid = $tmp[10];
	if ($taxidSpe{$taxid}) {
	    #next unless !($entrezName{$entrez});
	    #$entrezName{$entrez} = $genename;
	    $speGeneName = join ":", $taxidSpe{$taxid}, $genename;
	    if ($speGeneNameID{$speGeneName}) {
		$geneid = $speGeneNameID{$speGeneName};
		$intB = $geneid;
	    } else {
		#print "ERROR! Cannot find gene id for $speGeneName!\n";
	    }
	} else {
	    print "ERROR! Cannot find species info for $taxid!\n"
	}	
    }
    #----- done converting -------

    #------ start to process interaction data -----
    
    next unless (($intA ne "")&&($intB ne ""));
    $geneid = $intA;
    $mInteractor = $intB;
    $checkstring = join "mInt", $geneid, $mInteractor;
    next unless !($duplicateInfo{$checkstring}); 
    #if  ($geneid =~ /^FB/) {
    #	print "$checkstring\n";
    #}
    next unless $gName{$mInteractor}; #must be able to find interactor name
    $tempstring = join ":", $mInteractor, $gName{$mInteractor};
    if ($gMolInt{$geneid}) {
       $gMolInt{$geneid} =  join " \| ",  $gMolInt{$geneid}, $tempstring; 
    } else {
	$gMolInt{$geneid} = $tempstring;
	$totalGeneMInt++;
    }    
    $duplicateInfo{$checkstring} = 1; #record this pair of genetic interaction


    $geneid = $intB;
    $mInteractor = $intA;  
    $checkstring = join "mInt", $geneid, $mInteractor;
    next unless !($duplicateInfo{$checkstring});
    next unless $gName{$mInteractor}; #must be able to find interactor name
    $tempstring = join ":", $mInteractor, $gName{$mInteractor};
    if ($gMolInt{$geneid}) {
       $gMolInt{$geneid} =  join " \| ",  $gMolInt{$geneid}, $tempstring; 
    } else {
	$gMolInt{$geneid} = $tempstring;
	$totalGeneMInt++;
    }    
    $duplicateInfo{$checkstring} = 1; #record this pair of genetic interaction
    
}    
close (MINT);
print "Found Molecular Interaction info for $totalGeneMInt genes.\n";
#-------- Done getting molecular interaction ----




#----- print the source data table -------------
open (NAM, ">GeneNameDictionary.csv") || die "cannot open GeneName.csv!\n";
print NAM "AGR Gene ID\tAGR Gene Name\tSpecies\tEntrez ID\n";

open (OUT, ">Interaction.csv") || die "cannot open SimpleMineSourceData.csv!\n";
print OUT "AGR Gene ID\tGenetic Interaction\tMolecular\/Physical Interaction\n";

my ($gIntDes, $mIntDes);

foreach $geneid (@geneList) {
    next unless ($gName{$geneid});

    #--- fill in the blank fields ----------
    if ($gName{$geneid}) {
	$genename = $gName{$geneid};
    } else {
	$genename = "N.A.";
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

    
    $entrez = "N.A.";
    print OUT "$geneid\t$gIntDes\t$mIntDes\n";
    print NAM "$geneid\t$genename\t$s\t$entrez\n";
}
close (OUT);
close (NAM);

print "Total Human Genes: $totalHumanGene\n";
print "Total Mouse Genes: $totalMGIGene\n";
print "Total Rat Genes: $totalRGDGene\n";
print "Total Fish Genes: $totalZFINGene\n";
print "Total Fly Genes: $totalFBGene\n";
print "Total Worm Genes: $totalWBGene\n";
print "Total Yeast Genes: $totalSGDGene\n";
