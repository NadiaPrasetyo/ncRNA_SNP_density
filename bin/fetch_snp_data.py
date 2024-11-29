import pyBigWig

# Path to your BigBed file
bigbed_file = "path/to/your/file.bb"

# Number of lines to print
lines_to_print = 5

try:
    # Open the BigBed file
    bb = pyBigWig.open(bigbed_file)
    
    if bb.isBigBed():
        print("Opened BigBed file successfully.")
        print(f"Header Info: {bb.header()}")
        
        # Iterate over the entries
        entries_printed = 0
        for chrom in bb.chroms().keys():  # List of chromosomes
            intervals = bb.entries(chrom)
            if intervals:
                for interval in intervals:
                    print(f"Chromosome: {chrom}, Interval: {interval}")
                    entries_printed += 1
                    if entries_printed >= lines_to_print:
                        break
            if entries_printed >= lines_to_print:
                break
    else:
        print("The file is not a BigBed file.")
    
    # Close the file
    bb.close()

except Exception as e:
    print(f"An error occurred: {e}")
