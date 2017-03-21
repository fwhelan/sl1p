# /usr/bin/perl -w
#Author: Fiona J Whelan <whelanfj@mcmaster.ca>
#Last Modified Date: March 21 2017
use strict; use warnings;

if(!(exists $ENV{'SL1P_BIN'} && defined $ENV{'SL1P_BIN'})) {
	print "SL1P_BIN undefined; please set environment variable and try again.\n";
	exit;
}

$bin = $ENV{'SL1P_BIN'};
chdir($bin);
#install local cutadapt
print "###################INSTALLING CUTADAPT v1.8.1###################\n";
if (! -e 'cutadapt-1.8.1.tar.gz') {
	`wget https://pypi.python.org/packages/source/c/cutadapt/cutadapt-1.8.1.tar.gz`;
}
`tar -zxvf cutadapt-1.8.1.tar.gz`;
#install local AbundantOTU+
print "###################INSTALLING ABUNDANTOTU+ v0.93b###################\n";
if (! -e 'get.php?justdoit=yes&software=AbundantOTUplus0.93b.tar.gz') {
	print "wget http://omics.informatics.indiana.edu/mg/get.php\?justdoit\=yes\&software\=AbundantOTUplus0.93b.tar.gz\n";
	`wget http://omics.informatics.indiana.edu/mg/get.php?justdoit=yes\\&software=AbundantOTUplus0.93b.tar.gz`;
}
print "tar -zxvf get.php?justdoit=yes&software=AbundantOTUplus0.93b.tar.gz\n";
`tar -zxvf get.php?justdoit=yes\\&software=AbundantOTUplus0.93b.tar.gz`;
chdir('AbundantOTU+0.93b');
`./install`;
chdir('..');
#install local sickle
print "###################INSTALLING SICKLE###################\n";
if (! -e 'master.zip') {
	print "wget https://github.com/najoshi/sickle/archive/master.zip\n";
	`wget https://github.com/najoshi/sickle/archive/master.zip`;
}
print "tar -zxvf master.zip\n";
`tar -zxvf master.zip`;
chdir('..');
