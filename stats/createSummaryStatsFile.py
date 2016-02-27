import readBNGXML as bxml
import fnmatch
import os
import cPickle as pickle

def getFiles(directory,extension):
    """
    Gets a list of bngl files that could be correctly translated in a given 'directory'

    Keyword arguments:
    directory -- The directory we will recurseviley get files from
    extension -- A file extension filter
    """
    matches = []
    for root, dirnames, filenames in os.walk(directory):
        for filename in fnmatch.filter(filenames, '*.{0}'.format(extension)):
            filepath = os.path.abspath(os.path.join(root, filename))
            matches.append([filepath, os.path.getsize(os.path.join(root, filename))])

    #sort by size
    #matches.sort(key=lambda filename: filename[1], reverse=False)
    
    matches = [x[0] for x in matches]

    return matches

def extractStatsFromFile(fileName):
    fileStats = {}
    structures = bxml.parseFullXML(fileName)
    fileStats['compression'] = 1 - len(structures['molecules'])*1.0/len(structures['observables']) \
                               if len(structures['observables']) > 0 else 0
    fileStats['index'] = fileName.split('/')[-1].split('.')[0]
    fileStats['nreactions'] = len(structures['rules'])
    fileStats['nspecies'] = len(structures['observables'])
    fileStats['atomization'] = len([x for x in structures['rules'] if x[0].actionslen(structures[''])
    print structures['rules']
    return fileStats

if __name__ == "__main__":
    directory = 'non_curated'
    files = getFiles(directory, 'xml.xml')
    fileStats = []
    for bfile in files:
        fileStats.append(extractStatsFromFile(bfile))
    with open(os.path.join(directory, 'sortedD.dump'), 'wb') as f:
        pickle.dump(fileStats, f)

