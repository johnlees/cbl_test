use cbl::CBL;
use needletail::parse_fastx_file;
use needletail::Sequence;
use cbl::kmer::IntKmer;
use cbl::kmer::Kmer;

// define the parameters K and T
const K: usize = 25;
type T = u64; // T must be large enough to store $2k + \lg(2k)$ bits

fn main() {
    // create a CBL index with parameters K and T
    let mut cbl = CBL::<K, T>::new();

    let mut reader = parse_fastx_file("ERR556974.fastq.gz").unwrap();
    // for each sequence of the FASTA/Q file
    let mut cnt = 0;
    while let Some(record) = reader.next() {
        let seqrec = record.expect("Invalid record");
        let rc = seqrec.reverse_complement();
        for (_, kmer, _) in seqrec.canonical_kmers(K as u8, &rc) {
            let int_kmer = IntKmer::<K, T>::from_nucs(kmer);
            if cbl.contains(int_kmer) {
                cnt += 1;
            } else {
                cbl.insert(int_kmer);
            }
        }
    }
    println!("{} non-singletons", cnt);
}
