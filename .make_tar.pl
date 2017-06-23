#!/usr/bin/perl
#
# SYNOPSIS:
#   make_html.pl [OPTS] 
#
# DESCRIPTION:
#   Create an html page that provides a link to the latest generated
#   BioNetGen distribution.
#
# OPTIONS:
#   --platform   PLAT  : Choices are Linux, MacOSX, Windows


use strict;
use warnings;

# Perl Core Modules
use FindBin;
use File::Spec;
use Getopt::Long;
use Cwd ("getcwd");
use Config;
use File::Path qw(remove_tree);


# distribution version (default undefined)
my $version = '';
# platform choices are: (MacOSX or Linux or Windows)
my $platform = '';
# path to version file
my $path_to_version_file = '.';


GetOptions( 'help|h'        => sub { display_help(); exit(0); },
            'platform=s'    => \$platform);


    print " In perl  platform = ".$platform."\n";

my $zip_type  = '';
my $travis_os = '';
if ($platform eq "linux") {
  $zip_type = ".tar.gz";  $travis_os = "Linux";
    
  my $archive_file = "./dist/Atomizer-source-".$platform.$zip_type;
  print "\nCreating Atomizer source archive:\n";
  system("tar cvzf ${archive_file}  SBMLparser XMLExamples config gml2sbgn reactionDefinitions stats test utils Makefile requirements.txt twistedServer.py ");
  system("ls -lt dist");

} else {
  if ($platform eq "osx") {
  $zip_type = ".tar.gz";  $travis_os = "MacOSX";
    
  my $archive_file = "./dist/Atomizer-source-".$platform.$zip_type;
  print "\nCreating Atomizer source archive:\n";
  system("tar cvzf ${archive_file}  SBMLparser XMLExamples config gml2sbgn reactionDefinitions stats test utils Makefile requirements.txt twistedServer.py ");
  system("ls -lt dist");

  } else {
  
    if ($platform eq "Win32") {
      $zip_type = ".zip";  $travis_os = "Win32";
      system('copy  .\Win32_archive\sbmlTranslator.exe  .\dist\sbmlTranslator-Win32.exe');      
      my $archive_file = "./dist/Atomizer-source-Win32".$zip_type;
      
      print "\nCreating Atomizer-".$platform.".exe source archive:\n";
      system("7z a  ${archive_file} SBMLparser XMLExamples config gml2sbgn reactionDefinitions stats test utils Makefile requirements.txt twistedServer.py ");
      system("dir ${archive_file}");
      
    } else {
      if ($platform eq "Win64") {
        $zip_type = ".zip";  $travis_os = "Win64";
        system('copy  .\dist\sbmlTranslator.exe  .\dist\sbmlTranslator-Win64.exe');
        my $archive_file = "./dist/Atomizer-source-Win64".$zip_type;
        
        print "\nCreating Atomizer-".$platform.".exe source archive:\n";
        system("7z a  ${archive_file} SBMLparser XMLExamples config gml2sbgn reactionDefinitions stats test utils Makefile requirements.txt twistedServer.py ");
        system("dir ${archive_file}");
        
      } else {
        print "Invalid platform: ".$platform."\n";
        exit;
      }
    }
  }
}


exit;

# ########################################################################
#   HELP 
# ########################################################################
# display help menu
sub display_help
{
print <<END_HELP
make_html.pl
SYNOPSIS:
   make_html.pl [OPTS] 
DESCRIPTION:
   Create an html page that provides a link to the latest generated
   BioNetGen distribution.
OPTIONS:
   --platform   PLAT  : Choices are Linux, MacOSX, Windows
END_HELP

}
