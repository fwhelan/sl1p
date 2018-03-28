# /usr/bin/perl -w
#Author: Fiona J Whelan <whelanfj@mcmaster.ca>
#Last Modified Date: March 21 2017
use strict; use warnings;

if(!(exists $ENV{'SL1P_BIN'} && defined $ENV{'SL1P_BIN'})) {
	print "[sl1p]SL1P_BIN undefined; please set environment variable and try again.\n";
	exit;
}

my $bin = $ENV{'SL1P_BIN'};
chdir($bin);
#install local cutadapt
print "[sl1p]###################INSTALLING CUTADAPT v1.8.1###################\n";
if (! -e 'cutadapt-1.8.1.tar.gz') {
	`wget https://pypi.python.org/packages/source/c/cutadapt/cutadapt-1.8.1.tar.gz`;
}
`tar -zxvf cutadapt-1.8.1.tar.gz`;
#install local AbundantOTU+
print "[sl1p]###################INSTALLING ABUNDANTOTU+ v0.93b###################\n";
if (! -e 'get.php?justdoit=yes&software=AbundantOTUplus0.93b.tar.gz') {
	print "[sl1p]wget http://omics.informatics.indiana.edu/mg/get.php\?justdoit\=yes\&software\=AbundantOTUplus0.93b.tar.gz\n";
	`wget http://omics.informatics.indiana.edu/mg/get.php?justdoit=yes\\&software=AbundantOTUplus0.93b.tar.gz`;
}
print "[sl1p]tar -zxvf get.php?justdoit=yes&software=AbundantOTUplus0.93b.tar.gz\n";
`tar -zxvf get.php?justdoit=yes\\&software=AbundantOTUplus0.93b.tar.gz`;
chdir('AbundantOTU+0.93b');
`./install`;
chdir('..');
#install local sickle
print "[sl1p]###################INSTALLING SICKLE###################\n";
if (! -e 'master.zip') {
	print "wget https://github.com/najoshi/sickle/archive/master.zip\n";
	`wget https://github.com/najoshi/sickle/archive/master.zip`;
}
print "[sl1p]unzip master.zip\n";
`unzip master.zip`;
print "make sickle\n";
chdir('sickle-master');
`make`;
chdir('..');
`mv sickle-master/ sickle`;
chdir('..');
