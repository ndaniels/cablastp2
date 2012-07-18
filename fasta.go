package cablastp

import (
	"io"

	"code.google.com/p/biogo/io/seqio/fasta"
)

// ReadOriginalSeqs reads a FASTA formatted file and returns a slice of
// original sequences.
func ReadOriginalSeqs(fileName string) ([]*OriginalSeq, error) {
	reader, err := fasta.NewReaderName(fileName)
	if err != nil {
		return nil, err
	}

	sequences := make([]*OriginalSeq, 0)
	for i := 0; true; i++ {
		seq, err := reader.Read()
		if err == io.EOF {
			break
		}
		if err != nil {
			return nil, err
		}
		sequences = append(sequences, NewBiogoOriginalSeq(i, seq))
	}

	return sequences, nil
}

// ReadReferenceSeqs reads a FASTA formatted file and returns a slice of
// reference sequences.
func ReadReferenceSeqs(fileName string) ([]*ReferenceSeq, error) {
	reader, err := fasta.NewReaderName(fileName)
	if err != nil {
		return nil, err
	}

	sequences := make([]*ReferenceSeq, 0)
	for i := 0; true; i++ {
		seq, err := reader.Read()
		if err == io.EOF {
			break
		}
		if err != nil {
			return nil, err
		}
		sequences = append(sequences, NewBiogoReferenceSeq(i, seq))
	}

	return sequences, nil
}
