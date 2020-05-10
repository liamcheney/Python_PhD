import pandas as pd

def main():

    df = pd.read_csv('/Users/liamcheney/Desktop/Book2.csv')
    infile = open('/Users/liamcheney/Desktop/Untitled.txt','r').read().splitlines()

    sub = df[df['ST'].isin(infile)]
    sub = sub.sort_values('Contig Count', ascending=True)
    sub.to_csv('/Users/liamcheney/Desktop/Book3.csv')


    # infile = open('/Users/liamcheney/Desktop/assembly_result.xml').read()
    # infiles = infile.split('\n\n')
    # for file in infiles:
    #     lines = file.split('\n')
    #     GC = ''
    #     BI = ''
    #     BP = ''
    #     collcetion_date = ''
    #     for line in lines:
    #         if '<AssemblyAccession>' in line:
    #             GC = line.split('>')[1].split('<')[0]
    #         if '<BioSampleAccn>' in line:
    #             BI = line.split('>')[1].split('<')[0]
    #         if '<BioprojectAccn>' in line:
    #             BP = line.split('>')[1].split('<')[0]
    #         if ('collection' in line) and ('date' in line):
    #             collcetion_date = line.split('>')[1].strip().split('<')[0]
    #
    #
    #     print(GC, BI, BP, collcetion_date)

if __name__ == '__main__':
    main()