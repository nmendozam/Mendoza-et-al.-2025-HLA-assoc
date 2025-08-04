import argparse

parser = argparse.ArgumentParser(description='Modify last column of a file.')
parser.add_argument('input_file', type=str, help='input file name')
parser.add_argument('output_file', type=str, help='output file name')
args = parser.parse_args()

with open(args.input_file, 'r') as f:
    lines = f.readlines()
    for i in range(len(lines)):
        line = lines[i].strip().split()
        lines[i] = ' '.join(line[0:2]) + ' ' + ' '.join(line[-2:]) + '\n'

with open(args.output_file, 'w') as f:
    f.writelines("FID IID SEX LLI\n")
    f.writelines(lines)
