#i did not write.
#aims to split a large fasta file with many sequences into smaller files with less sequences.

from Bio import SeqIO

def batch_iterator(iterator, batch_size):
    entry = True  # Make sure we loop once

    while entry:
        batch = []

        while len(batch) < batch_size:
            try:
                entry = next(iterator)
            except StopIteration:
                entry = False

            if not entry:
                # End of file
                break

            batch.append(entry)

        if batch:
            yield batch

record_iter = SeqIO.parse('/Users/liam/Desktop/reference/GCF_000006745.1_ASM674v1_protein.faa', 'fasta')

for i, batch in enumerate(batch_iterator(record_iter, 49), start=1):
    filename = 'group_{}.fasta'.format(i)
    count = SeqIO.write(batch, filename, 'fasta')
    with open('/Users/liam/Desktop/Split_reference_faa/GCA_000006745_' + filename, 'w') as handle:
        SeqIO.write(batch,handle, 'fasta')

    print('Wrote {} records to {}'.format(count, filename))
