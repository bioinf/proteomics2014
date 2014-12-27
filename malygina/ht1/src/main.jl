require("ArgParse")

using ArgParse

push!(LOAD_PATH, dirname(@__FILE__()))

importall Clustering
importall ProfileAligner
importall DataReader
importall DataWriter

function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table s begin
        "--score-matrix", "-s"
            help = "an option with an argument"
            arg_type = String
            #required = true
        "--clustering", "-c"
            help = "clustering type: N for neighbour joining, U for UPGMA, W for WPGMA"
            arg_type = String
            default = "N"
            #required = false
        "--verbose"
          help = "show debug messages while processing data"
          action = :store_true
        "fasta-file"
            help = "file in .fasta format with a set of protein strings"
            required = true
        "output-file"
            help = "file in .fasta format for saving results"
            required = true
    end

    return parse_args(s)
end

function scoreFunc(p1 :: Profile{Float64}, p2 :: Profile{Float64})
  scoreprofiles(p1, p2)
end

function mergeFunc(p1 :: Profile{Float64}, p2 :: Profile{Float64})
  align(p1, p2)
end

strToProfiles(strings :: Vector{FastaRecord}) = [Profile{Float64}(record.sequence, record.description) :: Profile{Float64} for record in strings]


function main()
    parsed_args = parse_commandline()

    input_file = parsed_args["fasta-file"]
    output_file = parsed_args["output-file"]
    clustering = parsed_args["clustering"]
    score_matrix_file = parsed_args["score-matrix"]
    score_matrix = readMatrix(score_matrix_file)
    println("Score matrix loaded...")

    ProfileAligner.setScoringMatrix(score_matrix)
    fasta_sequences = readSequences(input_file)
    println("FASTA sequences loaded...")
    seq = strToProfiles( fasta_sequences)
    println("Sequences converted to profiles successfully, starting tree construction...")
    result = (clustering == "N" ? NeighbourJoining(seq, scoreFunc, mergeFunc) :
      clustering == "W" ? WPGMA(seq, scoreFunc, mergeFunc) : UPGMA(seq, scoreFunc, mergeFunc) )
    println("Tree constructed, writing to output_file...")
    writeSequences(output_file, getstrings(result))
    println("All done, go and do smthng else")
end

main()
