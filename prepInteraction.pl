#!/usr/bin/perl -w
use strict;

print "Create AGR SimpleMine source data table ...\n";

my @spe = ("HUMAN", "MGI", "RGD", "ZFIN", "FB", "WB", "SGD", "XBXL", "XBXT");
my %taxidSpe = ("NCBITaxon:9606" => "HUMAN", 
	        "NCBITaxon:10090" => "MGI", 
	        "NCBITaxon:10116" => "RGD", 
	        "NCBITaxon:7955" => "ZFIN", 
	        "NCBITaxon:7227" => "FB", 
		"NCBITaxon:6239" => "WB",
		"NCBITaxon:8355" => "XBXL",
#		"NCBITaxon:2697049" => "XBXT",		
		"NCBITaxon:8364" => "XBXT",		
	        "NCBITaxon:559292" => "SGD");
my $s; #species name
my ($line, $columns);
my @tmp;

#-----------Get Gene Symbols (long names like MAP-Homologous Protein)-------
open (LNAM, "/home/wen/agrSimpleMine/src/nameFile") || die "can't open nameFile!";
my ($geneid, $mod, $modID, $taxid, $genename, $longname, $synonym);
my %gModID; #Gene ID in its own Model Organism Database
my %gLName; #Long names for each gene id
my %gName; #public name for each gene id
my %gSpe; #species for each gene id
my %gSYNO; #synonym ids
my @geneList;
my $i = 0; #this index geneList
my $totalSynGene = 0;
my @findSyn;

while ($line = <LNAM>) {
    next unless ($line =~ /NCBITaxon/);
    next unless !($line =~ /NCBITaxon:2697049/); 
    chomp ($line);
    @findSyn = ();
    @findSyn = split /\t/, $line;
    $columns = @findSyn;
    if ($columns > 2) {	#there is a gene name
	$geneid = $findSyn[0];
	$geneList[$i] = $geneid; #add this gene to gene list.
	$i++;

	($mod, $modID) = split ':', $geneid;
	$gModID{$geneid} = $modID;
	
	$taxid = $findSyn[1];
	if ($taxidSpe{$taxid}) {
	    $s = $taxidSpe{$taxid};
	    $gSpe{$geneid} = $s;
	} else {
	    print "ERROR! Cannot find species for $taxid!\n";
	}
    	
	$genename = $findSyn[2];
	$gName{$geneid} = $genename; 

	if ($columns > 3) { #there is a long name
	    $longname = $findSyn[3];    
	    $gLName{$geneid} = $longname;
	}
    }
}
close (LNAM);
print "$i genes identified from the Gene Name and Symol file.\n";
#----------Done getting Gene Symbol --------------


#--------------Get Synonyms for Gene Names --------------
open (SYNO, "/home/wen/agrSimpleMine/src/synonymFile") || die "can't open synonymFile!";
$line = <SYNO>;

while ($line = <SYNO>) {
    next unless ($line =~ /NCBITaxon/);
    next unless !($line =~ /NCBITaxon:2697049/); 
    chomp ($line);
    @findSyn = ();
    @findSyn = split /\t/, $line;
    $columns = @findSyn;
    next unless ($columns == 3); 
    $geneid = $findSyn[0];
    $synonym = $findSyn[2]; 
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


#---------------- Get Genomic Resource IDs -----------------------------------
open (XREF, "/home/wen/agrSimpleMine/src/xrefFile") || die "can't open xrefFile!";

my %gUniPro; #UniProtKB ID for gene ID
my %gNCBI; #NCBI ID for gene ID
my %gPANTHER; #PANTHER ID for gene ID
my %gENSEM; #ENSEMBL ID for gene ID
my %gRefSeq; #RefSeq ID for gene ID
my %ncbiGene; #Gene ID for NCBI ID
my %uniproGene; #Gene ID for UniProtKB ID
my %refseqGene; #Gene ID for RefSeq ID
my ($resource, $rID);

while ($line = <XREF>) {
    chomp ($line);
    @tmp = ();
    @tmp = split /\t/, $line;
    $columns = @tmp;
    next unless ($columns == 3);
    $geneid = $tmp[0];

	if ($tmp[2] =~ /^UniProtKB:/) {
	    ($resource, $rID) = split ":", $tmp[2];
	    if ($gUniPro{$geneid}) {
		$gUniPro{$geneid} = join " | ",  $gUniPro{$geneid}, $rID;
	    } else {
		$gUniPro{$geneid} = $rID;
		$uniproGene{$rID} = $geneid;
	    }
	} elsif ($tmp[2] =~ /^NCBI_Gene:/) {
	    ($resource, $rID) = split ":", $tmp[2];
	    $gNCBI{$geneid} = $rID;
	    $ncbiGene{$rID} = $geneid;
	} elsif ($tmp[2] =~ /ENSEMBL:/) {
	    ($resource, $rID) = split ":", $tmp[2];
	    $gENSEM{$geneid} = $rID;
	} elsif ($tmp[2] =~ /PANTHER:/) {
	    ($resource, $rID) = split ":", $tmp[2];
	    $gPANTHER{$geneid} = $rID;
	} elsif ($tmp[2] =~ /RefSeq:/) {
	    ($resource, $rID) = split ":", $tmp[2];
	    $gRefSeq{$geneid} = $rID;
	    $refseqGene{$rID} = $geneid;
	}    
}    
close (XREF);
#---------------------- Done with Cross Reference  --------------

my %duplicateInfo; #check if the information was aready taken;
my $checkstring;

#------------ Get Genetic Interaction -----------------
my %gGeneticInt; #genetic interaction for each gene
open (GINT, "/home/wen/agrSimpleMine/src/geneticIntFile") || die "can't open geneticIntFile!";
my ($gInteractor, $intA, $intB, $tempstring);
my @splitName;

my $totalGeneGInt = 0; 
while ($line = <GINT>) {
    chomp ($line);
    next unless (($line =~ /^wormbase/) || ($line =~ /^flybase/)||($line =~ /^entrez/)||($line =~ /^uniprot/));

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
    } elsif (($line =~ /^entrez/)||($line =~ /^uniprot/)) {
	# find $intA
	($resource, $rID) = split ':', $tmp[0];
	if ($resource =~ /^entrez/) {
	    next unless ($ncbiGene{$rID});	
	    $intA = $ncbiGene{$rID};
	    #print "NCBI:$rID mapped to $intA\n";
	} else {
	    print "$rID\n";
	    next unless ($uniproGene{$rID});	
	    $intA = $uniproGene{$rID};	    
	}
	    
	# find $intB
	($resource, $rID) = split ':', $tmp[1];
	if ($line =~ /^entrez/) {
	    next unless ($ncbiGene{$rID});	
	    $intB = $ncbiGene{$rID};
	} else {
	    next unless ($uniproGene{$rID});	
	    $intB = $uniproGene{$rID};	    
	}
    }
    #----- done converting -------

    #------ start to process interaction data -----    
    next unless (($intA ne "")&&($intB ne ""));
    $geneid = $intA;
    $gInteractor = $intB;
    $checkstring = join "gInt", $geneid, $gInteractor;
    next unless !($duplicateInfo{$checkstring}); 
    next unless $gName{$gInteractor}; #must be able to find interactor name
    $tempstring = join ":", $gInteractor, $gName{$gInteractor};
    #$tempstring = $gInteractor;
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
    #$tempstring = $gInteractor;
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
my %gMolInt; #molecular interaction for each gene, excluding protein-protein interaction
open (MINT, "/home/wen/agrSimpleMine/src/molIntFile") || die "can't open molIntFile!";
open (PINT, ">ppInteraction.csv") || die "can't open ppInteraction.csv!";
my $mInteractor;

my $totalGeneMInt = 0;

while ($line = <MINT>) {
    chomp ($line);
    
    @tmp = ();
    @tmp = split /\t/, $line;

    my $col = 0;
    my $colName = "";
    if ($tmp[0] =~ /interactor/) {
	print PINT "$line\n";
        foreach $colName (@tmp) {
	    #print "$col $colName\n";
	    $col++;
	}     
    }
    
    next unless (($line =~ /^wormbase/)||($line =~ /^flybase/)||($line =~ /^entrez/)||($line =~ /^uniprot/));
    
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

    } elsif (($line =~ /^entrez/) || ($line =~ /^uniprot/)) {
	# find $intA
	($resource, $rID) = split ':', $tmp[0];
	if ($line =~ /^entrez/) {
	    next unless ($ncbiGene{$rID});	
	    $intA = $ncbiGene{$rID};
	} else {
	   # if  ($uniproGene{$rID}) {
		#do nothing
	   # } else {
	   #	print "Cannot find gene id for Uniprot $rID!\n";
	   # }
	    next unless ($uniproGene{$rID});	
	    $intA = $uniproGene{$rID};	    
	    #print "UniProtKB:$rID mapped to $intA\n";
	}
	    
	# find $intB
	($resource, $rID) = split ':', $tmp[1];
	if ($line =~ /^entrez/) {
	    next unless ($ncbiGene{$rID});	
	    $intB = $ncbiGene{$rID};
	} else {
	    next unless ($uniproGene{$rID});	
	    $intB = $uniproGene{$rID};	    
	}
    } 
    #----- done converting -------

    my $type = "mol";
    if (($tmp[20] =~ /protein/)&&($tmp[21] =~ /protein/)) {
	$type = "PP";
	print PINT "$line\n";
    }

    next unless ($type eq "mol");
    
    #------ start to process interaction data -----
    
    next unless (($intA ne "")&&($intB ne ""));
    $geneid = $intA;
    $mInteractor = $intB;
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
close (PINT);
print "Found Molecular Interaction info for $totalGeneMInt genes.\n";
#-------- Done getting molecular interaction ----


#--------- Get Protein-Protein interaction ---
open (MINT, "ppInteraction.csv") || die "can't open ppInteraction.csv!";
$totalGeneMInt = 0;
my %gPPInt; #protein-protein interaction for each gene

while ($line = <MINT>) {
    chomp ($line);
#  next unless (($line =~ /^wormbase/)||($line =~ /^flybase/)||($line =~ /^entrez/)||($line =~ /^uniprot/));        
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

    } elsif (($line =~ /^entrez/) || ($line =~ /^uniprot/)) {
	# find $intA
	($resource, $rID) = split ':', $tmp[0];
	if ($line =~ /^entrez/) {
	    next unless ($ncbiGene{$rID});	
	    $intA = $ncbiGene{$rID};
	} else {
	   # if  ($uniproGene{$rID}) {
		#do nothing
	   # } else {
	   #	print "Cannot find gene id for Uniprot $rID!\n";
	   # }
	    next unless ($uniproGene{$rID});	
	    $intA = $uniproGene{$rID};	    
	    #print "UniProtKB:$rID mapped to $intA\n";
	}
	    
	# find $intB
	($resource, $rID) = split ':', $tmp[1];
	if ($line =~ /^entrez/) {
	    next unless ($ncbiGene{$rID});	
	    $intB = $ncbiGene{$rID};
	} else {
	    next unless ($uniproGene{$rID});	
	    $intB = $uniproGene{$rID};	    
	}
    } 
    #----- done converting -------
    
    #------ start to process interaction data -----
    
    next unless (($intA ne "")&&($intB ne ""));
    $geneid = $intA;
    $mInteractor = $intB;
    $checkstring = join "pInt", $geneid, $mInteractor; #protein-protein interaction
    next unless !($duplicateInfo{$checkstring}); 
    next unless $gName{$mInteractor}; #must be able to find interactor name
    $tempstring = join ":", $mInteractor, $gName{$mInteractor};
    
    if ($gPPInt{$geneid}) {
       $gPPInt{$geneid} =  join " \| ",  $gPPInt{$geneid}, $tempstring; 
    } else {
	$gPPInt{$geneid} = $tempstring;
	$totalGeneMInt++;
    }    
    $duplicateInfo{$checkstring} = 1; #record this pair of genetic interaction

    $geneid = $intB;
    $mInteractor = $intA;  
    $checkstring = join "pInt", $geneid, $mInteractor;
    next unless !($duplicateInfo{$checkstring});
    next unless $gName{$mInteractor}; #must be able to find interactor name
    $tempstring = join ":", $mInteractor, $gName{$mInteractor};

    if ($gPPInt{$geneid}) {
       $gPPInt{$geneid} =  join " \| ",  $gPPInt{$geneid}, $tempstring; 
    } else {
	$gPPInt{$geneid} = $tempstring;
	$totalGeneMInt++;
    }    
    $duplicateInfo{$checkstring} = 1; #record this pair of genetic interaction
    
}    
close (MINT);
print "Found Protein-Protein Interaction info for $totalGeneMInt genes.\n";
#------ Done getting Protein-Protein interaction ----



#----- print the source data table -------------
#open (NAM, ">GeneNameDictionary.csv") || die "cannot open GeneNameDictionary.csv!\n";
#print NAM "AGR Gene ID\tAGR Gene Name\tSpecies\tEntrez ID\n";
#print NAM "AGR Gene ID\tSpecies\tGene Name\tGene Symbol\tNCBI\tUniProt\tENSEMBL\tPANTHER\n";

open (OUT, ">Interaction.csv") || die "cannot open Interaction.csv!\n";
print OUT "AGR Gene ID\tMOD Gene ID\tGene Name\tGene Symbol\tSpecies\tNCBI ID\tEnsembl ID\tUniprot ID\tPanther ID\tRefSeq ID\tSynonym\tGenetic Interaction\tOther Molecular Interaction\tProtein-Protein Interaction\n";

my ($gIntDes, $mIntDes, $pIntDes);
my ($uniproID, $ncbiID, $pantherID, $ensemblID, $refseqID);

foreach $geneid (@geneList) {

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
	$modID = "N.A.";
    }
    if  ($gUniPro{$geneid}) {
	$uniproID = $gUniPro{$geneid};
    } else {
	$uniproID = "N.A.";
    }
    if  ($gNCBI{$geneid}) {
	$ncbiID = $gNCBI{$geneid};
    } else {
	$ncbiID = "N.A.";
    }
    if  ($gPANTHER{$geneid}) {
	$pantherID = $gPANTHER{$geneid};
    } else {
	$pantherID = "N.A.";
    }
    if  ($gENSEM{$geneid}) {
	$ensemblID = $gENSEM{$geneid};
    } else {
	$ensemblID = "N.A.";
    }
    if  ($gRefSeq{$geneid}) {
	$refseqID = $gRefSeq{$geneid};
	#print "$geneid\t$refseqID\n";
    } else {
	$refseqID = "N.A.";
    }    
    if ($gSYNO{$geneid}) {
       $synonym = $gSYNO{$geneid};
    } else {
       $synonym = "N.A.";
    }  
    if ($gSpe{$geneid}) {
	$s = $gSpe{$geneid};
    } else {
	$s = "N.A.";
	print "ERROR! No species information for $geneid!\n";
    }
   
    #--- fill in the blank fields ----------
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

    if ($gPPInt{$geneid}) {
	$pIntDes = $gPPInt{$geneid};
	if ($geneid =~ /WBGene/) {
	    print "$geneid has ppInt $pIntDes\n";
	}
    } else {
	$pIntDes = "N.A.";
    }
    
    print OUT "$geneid\t$modID\t$genename\t$longname\t$s\t$ncbiID\t$ensemblID\t$uniproID\t$pantherID\t$refseqID\t$synonym\t$gIntDes\t$mIntDes\t$pIntDes\n";
    #print NAM "$geneid\t$s\t$genename\t$longname\t$synonym\t$ncbiID\t$uniproID\t$ensemblID\t$pantherID\n";
}
close (OUT);
#close (NAM);
