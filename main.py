import sys, rearrange

def main():
    rearrange.readSamFile(sys.argv[1], sys.argv[2], sys.argv[3])
    print('SAM file filtering and rearrangement complete.')

if __name__ == "__main__":
    main()
