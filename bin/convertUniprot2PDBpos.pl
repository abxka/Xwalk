#!/usr/bin/perl -w
# Author: Abdullah Kahraman
# Date: 03.03.2011

################################################################################
################################################################################
### Converts UniProt residue numbers into PDB residue numbers.               ###
################################################################################
################################################################################

use strict;
use warnings;

use Getopt::Long;

use LWP::UserAgent;

my (
    # variable for parameters which are read in from commandline
    $help,
    $xlFile,
    $pdbFile,
    $fastaFile,
    $needle,
    $pdbFile2,
    $force,
   );
##############################################################################
### read all needed parameters from commandline ##############################

&GetOptions(
    "xl=s"     => \$xlFile,    # cross-link distance file, with residue numbers from -fasta/-pdb2 
    "pdb=s"    => \$pdbFile,   # PDB file
    "fasta=s"  => \$fastaFile, # FASTA file
    "needle=s" => \$needle,    # path to the needle executable
    "pdb2:s"   => \$pdbFile2,  # 2nd PDB file instead of FASTA file
    "force:s"  => \$force,     # force the aligned residue types to be identical 
    "help!"    => \$help,      # print this help
) or die "\nTry \"$0 -h\" for a complete list of options\n\n";

##############################################################################

# help
if ($help) {printHelp(); exit}

################################################################################
# SETTINGS
###############################################################################

my $warningMsg = "";
my %aa3 = ("GLY" => "G",
	   "ALA" => "A", 
	   "VAL" => "V",
	   "LEU" => "L",
	   "ILE" => "I",
	   "PHE" => "F",
	   "PRO" => "P",
	   "TRP" => "W",
	   "ASN" => "N",
	   "GLN" => "Q",
	   "MET" => "M",
	   "CYS" => "C",
	   "THR" => "T",
	   "TYR" => "Y",
	   "SER" => "S",
	   "ARG" => "R",
	   "LYS" => "K",
	   "HIS" => "H",
	   "ASP" => "D",
	   "GLU" => "E",
	   "ASX" => "B",
	   "GLX" => "Z");

################################################################################
# SUBROUTINES
###############################################################################

###############################################################################
sub printHelp {
    # prints a help about the using and parameters of this scripts 
    # (execute if user types commandline parameter -h)
    # param:  no paramaters
    # return: no return value

    my (
	$usage,
	$sourceCode,
	@rows,
	$row,
	$option,
	$scriptInfo,
	$example,
       );

    $usage = "$0 -xl 6PGD_ECOLI_xl.txt -pdb 2zyaA.pdb -fasta 6PGD_ECOLI.fasta\n";

    print "\nUsage: " .  $usage . "\n";

    print "Valid options are:\n\n";
    open(MYSELF, "$0") or
      die "Cannot read source code file $0: $!\n";
    $sourceCode .= join "", <MYSELF>;
    close MYSELF;
    $sourceCode =~ s/^.+?\&GetOptions\(\n//s;
    $sourceCode =~ s/\n\).+$//s;
    @rows = split /\n/, $sourceCode;
    foreach $row (@rows){
        $option = $row;
	$option =~ s/\s+\"//g;
	$option =~ s/\"\s.+\#/\t\#/g;
	$option =~ s/=./\t<value> [required]/;
	$option =~ s/:./\t<value> [optional]/;
	$option =~ s/!/\t<non value> [optional]/;

	$row =~ s/^.*//;
	print "\t";
	printf("%-1s%-30s%-30s\n", "-",$option,$row);

    } # end of foreach $row (@rows)
    print "\n";
    print "Options may be abreviated, e.g. -h for --help\n\n";

    $example  = "$0";
}
################################################################################
sub mapUniprot2PDBseqNumber{
    my $newXlTable = "";
    (my $m1_1, my $m1_2, my $m2_1, my $m2_2) = &makeAlignment();
    my %aln2PDB = %$m1_1;
    my %aln2uniprot = %$m2_1;
    my %pdb2aln = %$m1_2;
    my %uniprot2aln = %$m2_2;
    (my $mapPDBnum, my $mapPDBname) = &getPDBresNumberName($pdbFile);
    my %mapPDBnum = %$mapPDBnum;
    my %mapPDBname = %$mapPDBname;

    my %mapPDBnum2;
    my %mapPDBname2;
    if(defined $pdbFile2){
	(my $mapPDBnum2, my $mapPDBname2) = &getPDBresNumberName($pdbFile2);
	%mapPDBnum2 = %$mapPDBnum2;
	%mapPDBname2 = %$mapPDBname2;
    }


    my $newDistFile = "";
    open(X, $xlFile) or die "Failed to open $xlFile";
    while(<X>){
	chomp($_);
	next if(/^#/);
	(my $index, my $fileName, my $atom1, my @rest) = split(/\t/);
	(my $resName1, my $resNo1, my $chainId1, my $atomName1) = split(/-/, $atom1);
	if(defined $pdbFile2){
	    foreach my $alnPos (sort keys %mapPDBnum2){
		if($mapPDBnum2{$alnPos} == $resNo1){
		    $resNo1 = $alnPos;
		    last;
		}
	    }
	}
	my $alnPos1 = $uniprot2aln{$resNo1};
	if(!defined $alnPos1){
	    $warningMsg .= "0. First atom in $atom1 does not exist in the PDB structure $pdbFile. ".
		           "Please check your cross-link file $xlFile\n";
	    next;
	}
	if($aln2PDB{$alnPos1} != -1){
	    $resNo1 = $mapPDBnum{$aln2PDB{$alnPos1}};
	    if(defined $force){
		if($mapPDBname{$aln2PDB{$alnPos1}} ne $resName1){
		    print STDERR "1st residue has non-identical mapping: ".
			         "$resName1 vs $mapPDBname{$aln2PDB{$alnPos1}}".
				 "\n";
		    next;
		}
	    }
	    $resName1 = $mapPDBname{$aln2PDB{$alnPos1}};
	}
	my $alnPos2;
	# XL info
	if($#rest >= 0){
	    (my $resName2, my $resNo2, my $chainId2, my $atomName2) = split(/-/, $rest[0]);
	    if(defined $pdbFile2){
		foreach my $alnPos (sort keys %mapPDBnum2){
		    if($mapPDBnum2{$alnPos} == $resNo2){
			$resNo2 = $alnPos;
			last;
		    }
		}
	    }

	    $alnPos2 = $uniprot2aln{$resNo2};
	    if(!defined $alnPos2){
		$warningMsg .= "0. Second atom in $rest[0] does not exist in the PDB structure $pdbFile. ".
		               "Please check your cross-link file $xlFile\n";
		next;
	    }
	    if($aln2PDB{$alnPos2} != -1){
		$resNo2 = $mapPDBnum{$aln2PDB{$alnPos2}};

		if(defined $force){
		    if($mapPDBname{$aln2PDB{$alnPos2}} ne $resName2){
			print STDERR "2nd residue has non-identical mapping: ".
			    "$resName2 vs $mapPDBname{$aln2PDB{$alnPos2}}".
			    "\n";
			next;
		    }
		}

		$resName2 = $mapPDBname{$aln2PDB{$alnPos2}};
	    }
	    if($aln2PDB{$alnPos1} != -1 && $aln2PDB{$alnPos2} != -1){
		$newDistFile .= "$index\t$fileName\t$resName1-$resNo1-$chainId1-$atomName1\t".
		                "$resName2-$resNo2-$chainId2-$atomName2";
		for(my $t=1; $t<=$#rest; $t++){
		    $newDistFile .= "\t$rest[$t]";
		}
		$newDistFile .= "\n";
	    }
	}
	# monolink info
	elsif($aln2PDB{$alnPos1} != -1 && $#rest == -1){
	    $newDistFile .= "$index\t$fileName\t$resName1-$resNo1-$chainId1-$atomName1\n";
	}

	# Error/Warning msg
	if($aln2PDB{$alnPos1} == -1){
	    if($#rest >= 0 && $aln2PDB{$alnPos2} == -1){
		$warningMsg .= "1. Both UniProt residue numbers $alnPos1 and $alnPos2 could ".
	                       "not be found in the PDB structure.\n";
	    } else{
		# monoLink not found msg
		$warningMsg .= "2. The UniProt residue number $alnPos1 could ".
		               "not be found in the PDB structure.\n";
	    }
	}
	# second atom could not be found
	elsif($#rest >= 0 && $aln2PDB{$alnPos2} == -1){
	    $warningMsg .= "3. The UniProt residue number $alnPos2 could ".
		           "not be found in the PDB structure.\n";
	}
    }
    close(X);
    return $newDistFile;
}
################################################################################
sub makeAlignment(){
    my $fasta1 = &createFastaFile4PDBinput($pdbFile);
    my $fasta2;
    if(defined $fastaFile){
	$fasta2 = $fastaFile;
    } else {
	$fasta2 = &createFastaFile4PDBinput($pdbFile2);
    }
    my $aln = $pdbFile;
    $aln =~ s/.*\///;
    $aln =~ s/(.*)\..*/$1.aln/;
    my $command = "$needle $fasta1 $fasta2 -gapopen 10 -gapextend 0.5 $aln";
    print STDERR "$command\n";
    system(`$command >& /dev/null`) or 
	die("Error while performing Needleman-Wunsch alignment: $!");

    return &readAlignment($aln);
}
################################################################################
sub createFastaFile4PDBinput(){
    (my $pdbFile) = @_;

    my $id = $pdbFile;
    $id =~ s/.*\///;
    $id =~ s/(.*)\..*/$1/;
    my $fasta = "$id.fasta";
    my @a = split(/\n/,`grep ^ATOM $pdbFile`);
    open(F,">$fasta") or 
        die ("Couldn't create sequence file \"$fasta\" from PDB file: $!");
    print F ">$id.pdb\n";
    my %uniqAA;
    foreach my $l (@a){
	if($l=~/^ATOM/){
	    my $resName = substr($l, 17, 3);
	    my $resNo   = substr($l, 22, 4);
	    my $chainId = substr($l, 21, 1);
	    my $iCode = substr($l, 26, 1);
	    next if($iCode ne " ");
	    my $aa = "$resName-$resNo-$chainId";
	    if(!exists $uniqAA{$aa}){
		$uniqAA{$aa} = $aa;
		if(exists $aa3{$resName}){
		    print F $aa3{$resName};
		}
		else{
		    print F "X";
		}
	    }
	}
    }
    print F "\n";
    close(F);

    return $fasta;
}
################################################################################
sub readAlignment(){
    (my $aln) = @_;

    my $id = $pdbFile;
    $id =~ s/.*\///;
    $id =~ s/(.*)\..*/$1/;

    open(F,$aln) or die("Failed to open alignment file \"$aln\": $!");

    my $h1;
    my $h2;
    while(<F>){
	   chomp($_);
	   if(/^[A-Za-z0-9]/){
	       my $id1=substr("$id.pdb",0,13);
	       if(/^$id1/){
		      $_=~s/.*\d\s([A-Za-z\-])/$1/;
		      $_=~s/\s+\d+//;
		      $h1 .= $_;
	       }
	       else{
		      $_=~s/.*\d\s([A-Za-z\-])/$1/;
		      $_=~s/\s+\d+//;
		      $h2 .= $_;
	       }
	   }
    }
    close(F);
    if(length($h1) != length($h2)){
	   die("Error while reading alignment. Alignment sequence lengths do not match!\n");
    }

    my %map1_1; # alignment position vs sequence position
    my %map1_2; # sequence position vs alignment position
    my %map2_2; # alignment position vs sequence position
    my %map2_1; # sequence position vs alignment position
    my $j1=0; # position in sequence
    my $j2=0; # position in sequence
    for(my $i=0; $i < length($h1); $i++){
	   if(substr($h1, $i, 1) ne "-"){
	       $j1++;
	       $map1_1{$i+1} = $j1;
	       $map1_2{$j1} = $i+1;
	    } else {
	       $map1_1{$i+1} = -1;
    	}
	
    	if(substr($h2, $i, 1) ne "-"){
	        $j2++; 
	        $map2_1{$i+1} = $j2;
	       $map2_2{$j2} = $i+1;
        } else {
	        $map2_1{$i+1} = -1;
	    }
    }
    return(\%map1_1, \%map1_2, \%map2_1, \%map2_2);
}
################################################################################
sub getPDBresNumberName(){
    (my $pdbFile) = @_;
    my @a = split(/\n/,`grep ^ATOM $pdbFile`);
    my %uniqAA;
    my %mapNo;
    my %mapName;
    foreach my $l (@a){
	   if($l=~/^ATOM/){
	       my $resName = substr($l, 17, 3);
	       $resName =~ s/\s//g;
	       my $resNo   = substr($l, 22, 4);
	       $resNo =~ s/\s//g;
	       my $chainId = substr($l, 21, 1);
	       $chainId =~ s/\s//g;
	       my $iCode = substr($l, 26, 1);
	       next if($iCode ne " ");
	       my $aa = "$resName-$resNo-$chainId-$iCode";
	       if(!exists $uniqAA{$aa}){
		   $uniqAA{$aa} = $aa;
		   $mapNo{keys (%uniqAA)} = $resNo;
		   $mapName{keys (%uniqAA)} = $resName;
	       }
	   }
    }
    return(\%mapNo, \%mapName);
}
################################################################################
# MAIN #########################################################################
################################################################################

# check whether all parameter are real files
if(defined $xlFile and defined $pdbFile and (defined $fastaFile or defined $pdbFile2) and
   defined $needle){
    die "$xlFile does not exist\n" unless(-e $xlFile);
    die "$pdbFile does not exist\n" unless(-e $pdbFile);

    if(defined $fastaFile) {
	die "$fastaFile does not exist\n" unless(-e $fastaFile);
    }
    else{
	unless(defined $pdbFile2){
	    die "Please specify -fasta or -pdb2\n";
	}
	else {
	    die "$pdbFile2 does not exist\n" unless(-e $pdbFile2);
	}
    }
    die "Please specify the location of the needle executable\n" unless(defined $needle);
    

# create fasta file from pdb file
    my $newDistFile = &mapUniprot2PDBseqNumber();
    print STDERR $warningMsg;
    print $newDistFile;
    
}
else {
    die "\nTry $0 -help to get a full list of parameters.\n\n";
}

