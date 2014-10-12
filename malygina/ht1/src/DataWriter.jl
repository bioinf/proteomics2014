module DataWriter
  export writeSequences

  function writeSequences(output_file_name :: String, sequences :: Vector{ASCIIString})

    output_file = open(output_file_name, "w")
    for s in sequences
      write(output_file, s)
      write(output_file, '\n')
    end
    close(output_file)
  end

end
