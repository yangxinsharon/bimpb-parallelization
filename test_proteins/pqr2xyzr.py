import sys
import re

def main():
    # Run "python pqr2xyzr.py <PDBID>" with replacing the <PDBID> with the PDB ID of your protein.
    #  This is case sensisitve 
    PDBID = sys.argv[1]
    pqr_file = PDBID+'.pqr'
    xyzr_file = PDBID+'.xyzr'
    readpqr = open(pqr_file,'r')
    wrtxyzr = open(xyzr_file,'w')

    pqrlines = readpqr.readlines()
    for line in pqrlines:
        if re.match(r'^(ATOM).+$', line):
            items = line.split()
            x = items[5]
            y = items[6]
            z = items[7]
            r = items[9]
            wrtxyzr.write(' {:>7,.3f} {:>7,.3f} {:>7,.3f} {:>7,.4f} '.format(float(x),float(y),float(z),float(r)))
            wrtxyzr.write('\n')


if __name__=="__main__":
    main()
