package cablastp

import (
	"compress/gzip"
	"io"
	"os"
	"strings"
	"bytes"

	"github.com/TuftsBCB/io/fasta"
	"github.com/TuftsBCB/seq"
)

// ReadOriginalSeq is the value sent over `chan ReadOriginalSeq` when a new
// sequence is read from a fasta file.
type ReadOriginalSeq struct {
	Seq *OriginalSeq
	Err error
}

// ReadOriginalSeqs reads a FASTA formatted file and returns a channel that
// each new sequence is sent to.
func ReadOriginalSeqs(
	fileName string,
	ignore []byte,
) (chan ReadOriginalSeq, error) {
	var f io.Reader
	var err error

	f, err = os.Open(fileName)
	if err != nil {
		return nil, err
	}
	if strings.HasSuffix(fileName, ".gz") {
		f, err = gzip.NewReader(f)
		if err != nil {
			return nil, err
		}
	}

	reader := fasta.NewReader(f)
	seqChan := make(chan ReadOriginalSeq, 200)
	go func() {
		for i := 0; true; i++ {
			sequence, err := reader.Read()
			if err == io.EOF {
				close(seqChan)
				break
			}
			if err != nil {
				seqChan <- ReadOriginalSeq{
					Seq: nil,
					Err: err,
				}
				close(seqChan)
				break
			}
			for i, residue := range sequence.Residues {
				for _, toignore := range ignore {
					if toignore == byte(residue) {
						sequence.Residues[i] = 'X'
						break
					}
				}
			}
			seqChan <- ReadOriginalSeq{
				Seq: NewFastaOriginalSeq(i, sequence),
				Err: nil,
			}
		}
	}()
	return seqChan, nil
}

func ReduceQuerySeqs(
	query *bytes.Reader) (*bytes.Reader, error) {
		buf := new(bytes.Buffer)
		f := fasta.NewWriter(buf)
		reader := fasta.NewReader(query)
		for i := 0; true; i++ {
			sequence, err := reader.Read()
			if err == io.EOF {
				break
			}
			if err != nil {
				return nil, err
			}
			rs := Reduce(sequence.Bytes())
			n := sequence.Name

			result := seq.NewSequenceString(n, string(rs))
			f.Write(result)
		}
		
		
		return bytes.NewReader(buf.Bytes()), nil
	}
