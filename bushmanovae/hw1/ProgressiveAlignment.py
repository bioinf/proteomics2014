__author__ = 'letovesnoi'

import PipelineUtils
import FastaUtils
import ParsingMatricesUtils
import Clustering
import AlignmentsUtils
import ProfilesUtils

acids = ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K',
         'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V', 'B', 'Z', 'X', '*']

def main():
    global acids

    args = PipelineUtils.getArgs()

    matrix = ParsingMatricesUtils.parseMatrix(args.matrix)

    sequences = FastaUtils.readFASTA(args.fasta)

    EDM = Clustering.constructEDM(sequences, matrix)
    #print EDM

    tree = []
    alignments = []
    profiles = []
    numbers = []

    for seq in sequences:
        alignments.append([seq])
        profiles.append(ProfilesUtils.getProfileSequence(seq, acids))
        numbers.append(1)

    while len(EDM) > 0:
        EDM, sequences, numbers, tree, iVertex1, iVertex2 = Clustering.updateEDMSTI(EDM, sequences, tree, numbers, args.type)

        #print tree

        profile1 = profiles[iVertex1]
        profile2 = profiles[iVertex2]
        number1 = numbers[iVertex1]
        number2 = numbers[iVertex2]
        alignment1 = alignments[iVertex1]
        alignment2 = alignments[iVertex2]

        tracebackProfiles = AlignmentsUtils.getTracebackProfiles(matrix, profile1, profile2)


        profiles.pop(min(iVertex1, iVertex2))
        profiles.pop(max(iVertex1, iVertex2) - 1)
        numbers.pop(min(iVertex1, iVertex2))
        numbers.pop(max(iVertex1, iVertex2) - 1)
        alignments.pop(min(iVertex1, iVertex2))
        alignments.pop(max(iVertex1, iVertex2) - 1)

        tmpProfile = {}
        for acid in acids:
            tmpProfile[acid] = []
        tmpAlignment = []
        for i in range(len(alignment1) + len(alignment2)):
            tmpAlignment.append('')

        profiles.append(ProfilesUtils.getGlueProfile(acids, tmpProfile, profile1, profile2, number1, number2, tracebackProfiles,
                                                     len(profile1['*']) - 1, len(profile2['*']) - 1))

        alignments.append(AlignmentsUtils.getGlueAlignment(tracebackProfiles, tmpAlignment, alignment1, alignment2,
                                                           len(alignment1[0]) - 1, len(alignment2[0]) - 1))
        numbers.append(number1 + number2)

    likeNewick = sequences[0]
    alignment = alignments[0]

    with open(args.out, 'w') as fout:
        fout.write(likeNewick + '\n\n')
        #print(likeNewick)
        for seq in alignment:
            fout.write(seq + '\n')
            #print seq

main()