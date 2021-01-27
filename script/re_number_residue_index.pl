#! /usr/bin/perl -w

 require 5.003; # need this version of Perl or newer
 use English; # use English names, not cryptic ones
 use FileHandle; # use FileHandles instead of open(),close()
 use Carp; # get standard error / warning messages
 use strict; # force disciplined use of variables
 use Cwd;
 use Cwd 'abs_path';
 use Scalar::Util qw(looks_like_number);
 sub get_len($);
 if(@ARGV<2)
 {
     print "This script will just re-number the residue and atom index!\n";
     print "For example:\n";
     print "perl $0 ../T0787_refine/T0787_refined3.pdb ../T0787_refine/T0787_refined3.pdb.new\n";
     exit(0);
 }
##############standard Amino Acids (3 letter <-> 1 letter)#######
my(%amino)=();
$amino{"ALA"} = 'A';
$amino{"CYS"} = 'C';
$amino{"ASP"} = 'D';
$amino{"GLU"} = 'E';
$amino{"PHE"} = 'F';
$amino{"GLY"} = 'G';
$amino{"HIS"} = 'H';
$amino{"ILE"} = 'I';
$amino{"LYS"} = 'K';
$amino{"LEU"} = 'L';
$amino{"MET"} = 'M';
$amino{"ASN"} = 'N';
$amino{"PRO"} = 'P';
$amino{"GLN"} = 'Q';
$amino{"ARG"} = 'R';
$amino{"SER"} = 'S';
$amino{"THR"} = 'T';
$amino{"VAL"} = 'V';
$amino{"TRP"} = 'W';
$amino{"TYR"} = 'Y';

### here is some ambiguous amino acid
#$amino{"ASX"} = 'B';
#$amino{"GLX"} = 'Z';
#$amino{"XLE"} = 'J';
#$amino{"XAA"} = 'X';
###################################################################
############## 1 letter to 3 letter ################
my(%amino_1_to_3) = ();
my($new_key);

foreach $new_key (keys %amino)  
{
	$amino_1_to_3{$amino{$new_key}} = $new_key;
}
####################################################
 my($addr_input)=$ARGV[0];
 my($addr_output)=$ARGV[1];

 my($file,$path_one_target,$path_out);
 my(@files,@tem);
 my($IN,$OUT,$line,$i);
 my($tm_chain);
 

 ######## first get the TMalign output #############
 ############# test #############
 #print $alignments[0]."".$alignments[1]."".$alignments[2]."";
 ################################


 ######## now, make a new template filling the aligned residues ############
 my($current_index) = "NULL";
 $i=-1;
 my($residue_index)= 1;
 my($atom_index)=0;
 my($tmp_index);
 
 $OUT = new FileHandle ">$addr_output";
 $IN = new FileHandle "$addr_input";
 while(defined($line=<$IN>))
 {
	 @tem = split(/\s+/,$line);
	 if($tem[0] ne "ATOM")
	 {
		 next;
	 }
	 $tmp_index = substr($line,22,4);
	 $tmp_index=int($tmp_index);


substr($line,21,1) = " ";                 # set the chain id empty

	 if($current_index ne $tmp_index)
	 {# this is a new residue
		 $i++; 
                 
                 $current_index = $tmp_index;
                 $atom_index++;
		 substr($line,6,5) = sprintf("%5.0f",$residue_index);
                 $residue_index++;
                 substr($line,22,4) = sprintf("%4.0f",$atom_index);
                 
		 print $OUT $line;
 #                print $OUT_align $hash_res{$i};         
	 }
	 else
	 {#
                 substr($line,6,5) = sprintf("%5.0f",$residue_index);
                 $residue_index++;
                 substr($line,22,4) = sprintf("%4.0f",$atom_index);
                 print $OUT $line;

	 }
        
 }
 
 $IN->close();
 $OUT->close();


