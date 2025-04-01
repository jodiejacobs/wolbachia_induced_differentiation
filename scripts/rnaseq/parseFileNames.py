import argparse
import pandas as pd

parser = argparse.ArgumentParser()
parser.add_argument('-i')
parser.add_argument('-e')

args = parser.parse_args()

filenames = args.i
extension = args.e


names = []
prefix = ""

file = open (filenames, 'r')

for line in file:
    if len(line.strip()) != 0:
        l = line.strip('\n')
        
        if line[:2] == "./":
            prefix = l[2:-1]
        
        elif l[-len(extension):] == extension:
            
            names.append(l.strip(extension))

print(names)
# for n in names:
#     print(n)





