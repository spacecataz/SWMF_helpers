#!/usr/bin/perl -w
#cleanup.pl

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#This script "cleans up" an swmf run directory by removing output files,
#log files and other unnecessary items.  When called, the script asks 
#the user what files should be removed.
#
#Options:
#  -all  Overrides the prompting and performs full cleanup.
#
#DTW 2007
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

use strict;

# Declare subroutines
sub print_help;

my $version = 1.1;

#------------------------------------------------------------------------
#Handle all arguments.
#------------------------------------------------------------------------
my $allflag = 0;

foreach(@ARGV){
    #Print help if -h argument is present.
    print_help if (($_ eq '-h') or ($_ eq '-H'));

    #If the all argument is present, set the allflag.
    if($_ eq '-all'){
	$allflag = 1;
	print "Performing full clean.\n";
    }

}


#------------------------------------------------------------------------
#Remove all emacs cache files (tilde files).
#------------------------------------------------------------------------
while(<*~>){unlink $_ or warn "Can't remove file $_ : $!\n"}


#------------------------------------------------------------------------
#Examine run directory and job files to get list of deletable crap.
#------------------------------------------------------------------------

#Create lists to hold names of files to be deleted.
my @loglist;
my @pbslist;
my @plotlist;
my @restlist;

#Loop through all job.* files and examine them.
foreach(<job.*>){
    open(FH, '<', $_) or die "Couldn't read job file $_: $!\n";
    
    #Logfiles are named in the job script; pbs logfiles are named after
    #the PBS job name.  Let's extract both from every file.
    while(<FH>){
	if(/^\#PBS -N\s*(\S+)/){
	    next if($1 eq 'SWMF');
	    push @pbslist, $1.'.e*';
	    push @pbslist, $1.'.o*';
	}
	if(/mpirun.+(\>|\>\&)\s*(\S+)/){
	    #Check for syntax that appends the date to the log file.
	    #If it's there, we have to rely on wild cards to delete 
	    #the correct files.
	    if($2 =~ /\`/){
		push @loglist, 'log.*';
	    }else{
		push @loglist, $2;
	    }
	}
    }
    close(FH);
}

#Examine pwd, figure out what components are installed.
#Add any restart directories to the list.
while(<*>){
    if($_ eq 'GM'){push @plotlist, 'GM/IO2/*'}
    if($_ eq 'IM'){push @plotlist, 'IM/plots/*'}
    if($_ eq 'IE'){push @plotlist, 'IE/ionosphere/*'}
    if($_ eq 'PS'){push @plotlist, 'PS/Output/*'}
    if($_ eq 'UA'){push @plotlist, 'UA/data/*'}
    if(/^RESTART_/){push @restlist, $_}
}


#------------------------------------------------------------------------
#Review clean-up list.
#------------------------------------------------------------------------

#Check the results of the analysis:
print "\nOutput directories to clean:\n";
print "$_; " foreach(@plotlist);
print "\nPBS Logfiles to remove:\n";
print "$_; " foreach(@pbslist);
print "\nLogfiles to be removed:\n";
print "$_; " foreach(@loglist);
print "\nRestart directories to be removed:\n";
print "$_; " foreach(@restlist);


#------------------------------------------------------------------------
#Delete files; ask first if -all was missing from argument list.
#------------------------------------------------------------------------
my @dellist;

#Add to list names of files that should -always- be deleted.
push @dellist, 'SWMF.SUCCESS';

if($allflag){
    push @dellist, @plotlist, @pbslist, @loglist;
    system("rm -rf @restlist");
}else{
    print "\nDelete all plot and output files? (y or n);";
    chomp(my $ans = (<STDIN>));
    push @dellist, @plotlist unless($ans ne 'y');

    print "Delete all PBS log files? (y or n);";
    chomp($ans = (<STDIN>));
    push @dellist, @pbslist unless($ans ne 'y');

    print "Delete all SWMF log files? (y or n);";
    chomp($ans = (<STDIN>));
    push @dellist, @loglist unless($ans ne 'y');

    print "Delete all Restart directories? (y or n)";
    chomp($ans = (<STDIN>));
    system("rm -rf @restlist") unless ($ans ne 'y');
}

print "CLEANING UP!\n";
unlink (<@dellist>) or warn "Cannot remove files: $!\n";

#===============================================================

sub print_help
{
    print "\n    cleanup.pl version $version\n";
	print "    DTW, 2007\n\n";

    print<<EOF;
    This script quickly cleans up the current run directory
    so that another simulation can be run with out output and 
    restart files from the previous simulation getting in the 
    way.  It removes restarts, system logs, SWMF logs, and 
    output files from the folders in the PWD.  It will ask for
    permission to delete each group so that you can remove a 
    subset of the list to be deleted rather than blindly wiping
    the files.

    USAGE--
         Run this script only in the run directory of the SWMF.
	 If the -all flag is not present, the script will ask
	 for permission to remove each group of files.

	 -H or -h -- Help.  Print this help screen.

	 -all -- Do not ask before removing all files.

EOF

exit
}
