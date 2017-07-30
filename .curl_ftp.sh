
if [ "${TRAVIS_OS_NAME}" = "linux" ]; then
  curl -v -T ./dist/sbmlTranslator-$TRAVIS_OS_NAME               -u $FTP_USER_4:$FTP_PASSWORD ftp://ftp.drivehq.com/d_data/d_travis/  
else
  curl -v -T ./dist/sbmlTranslator-$TRAVIS_OS_NAME               -u $FTP_USER_5:$FTP_PASSWORD ftp://ftp.drivehq.com/d_data/d_travis/  
fi


#    The FTP_USER and FTP_PASSWORD strings are actually defined in the env: global: -secure: strings in travis.yml.  Those strings
#  Were generated using a program that may be installed on Linux.  Documentation can be found at: https://docs.travis-ci.com/user/encryption-keys/
#  To make a long story short, the program may be installed on Linux with:      gem install travis
#  the command to generate the FTP_PASSWORD string is
#
#              travis encrypt  FTP_USER="password" --skip-version-check -r RuleWorld/nfsim
#
#  The password should be surroundedd by double quotes, and the -r parameter indicates the repository for which the password
#  will be used.
#

#  The  -r  parameter is important.  The parameter should be set as follows for the various build processes:
#     Atomizer:    -r  RuleWorld/atomizer
#     NFsim:       -r  RuleWorld/nfsim
#     BioNetGen:   -r  RuleWorld/bionetgen
#  There is no need to specify a subdirectory when setting -r  For example, much of the build process for BioNetGen
#  occurs in the  bng2  subdirectory of BioNetGen.  But you will get an error if you use   -r RuleWorld/bionetgen/bng2  
#  to generate a password that will be used in that subdirectory.
