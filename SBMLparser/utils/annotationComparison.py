import annotationExtractor as annEx
import argparse
import fnmatch
import argparse
import os
import progressbar


def defineConsole():
    parser = argparse.ArgumentParser(description='SBML to BNGL translator')
    parser.add_argument('-m1','--model1',type=str,help='input SBML model1',required=True)
    parser.add_argument('-m2','--model2',type=str,help='input SBML model2',required=True)
    #parser.add_argument('-o','--output-file',type=str,help='output SBML file',required=True)
    return parser    


def getFiles(directory, extension):
    """
    Gets a list of <*.extension> files. include subdirectories and return the absolute 
    path. also sorts by size.
    """
    matches = []
    for root, dirnames, filenames in os.walk(directory):
        for filename in fnmatch.filter(filenames, '*.{0}'.format(extension)):
            matches.append([os.path.join(os.path.abspath(root), filename),os.path.getsize(os.path.join(root, filename))])

    #sort by size
    matches.sort(key=lambda filename: filename[1], reverse=False)
    
    matches = [x[0] for x in matches]

    return matches

def annotationComparison(model1 ,model2):
    annotationExtractor = annEx.AnnotationExtractor(model1)
    modelAnnotations1 = annotationExtractor.getModelAnnotations()
    #elementalMolecules = [x for x in annotationExtractor.sct if annotationExtractor.sct[x] == []]
    complexMolecules1 = [x for x in annotationExtractor.sct if annotationExtractor.sct[x] != []]

    annotationDict1 = annotationExtractor.getAnnotationSystem()

    annotationExtractor = annEx.AnnotationExtractor(model2)
    modelAnnotations2 = annotationExtractor.getModelAnnotations()
    #elementalMolecules = [x for x in annotationExtractor.sct if annotationExtractor.sct[x] == []]
    complexMolecules2 = [x for x in annotationExtractor.sct if annotationExtractor.sct[x] != [] and x in complexMolecules1]
    annotationDict2 = annotationExtractor.getAnnotationSystem()

    error = 0
    for entry in  complexMolecules2:
        if not set([x for x in annotationDict2[entry]['BQB_HAS_PART'] if 'uniprot' in x]).issubset(set([x for x in annotationDict1[entry]['BQB_HAS_PART'] if 'uniprot' in x])):
            error+=1

        if not set([x for x in annotationDict2[entry]['BQB_IS_VERSION_OF'] if 'uniprot' in x]).issubset(set([x for x in annotationDict1[entry]['BQB_IS_VERSION_OF'] if 'uniprot' in x])):
            error += 1
    if error > 0:
        print model1.split('/')[1], error*1.0/len(complexMolecules2)


def batchAnnotationComparison(removedAnnotationsDir, referenceDir):
    referenceFiles = getFiles(referenceDir,'xml')
    progress = progressbar.ProgressBar()

    for fileIdx in progress(range(len(referenceFiles))):
        file = referenceFiles[fileIdx]
        annotationComparison(os.path.join(removedAnnotationsDir, file.split('/')[-1]), file)

if __name__ == "__main__":
    batchAnnotationComparison('annotationsRemoved', '../XMLExamples/curated')
    '''
    parser = defineConsole()
    namespace = parser.parse_args()

    annotationExtractor = annEx.AnnotationExtractor(namespace.model1)
    modelAnnotations1 = annotationExtractor.getModelAnnotations()
    #elementalMolecules = [x for x in annotationExtractor.sct if annotationExtractor.sct[x] == []]
    complexMolecules1 = [x for x in annotationExtractor.sct if annotationExtractor.sct[x] != []]

    annotationDict1 = annotationExtractor.getAnnotationSystem()

    annotationExtractor = annEx.AnnotationExtractor(namespace.model2)
    modelAnnotations2 = annotationExtractor.getModelAnnotations()
    #elementalMolecules = [x for x in annotationExtractor.sct if annotationExtractor.sct[x] == []]
    complexMolecules2 = [x for x in annotationExtractor.sct if annotationExtractor.sct[x] != [] and x in complexMolecules1]
    annotationDict2 = annotationExtractor.getAnnotationSystem()

    for entry in  complexMolecules2:
        if not set([x for x in annotationDict2[entry]['BQB_HAS_PART'] if 'uniprot' in x]).issubset(set([x for x in annotationDict1[entry]['BQB_HAS_PART'] if 'uniprot' in x])):
            print '----- HAS PART'
            print entry
            print annotationDict1[entry]
            print annotationDict2[entry]

        if not set([x for x in annotationDict2[entry]['BQB_IS_VERSION_OF'] if 'uniprot' in x]).issubset(set([x for x in annotationDict1[entry]['BQB_IS_VERSION_OF'] if 'uniprot' in x])):
            print '----- IS VERSION OF'
            print entry
            print annotationDict1[entry]
            print annotationDict2[entry]

    '''