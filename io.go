package cablastp

import (
	"bytes"
	// "compress/gzip"
	"encoding/binary"
	"encoding/csv"
	"fmt"
	"io"
	"os"
	"strconv"
	"time"

	"github.com/TuftsBCB/io/fasta"
)

func (coarsedb *CoarseDB) readFasta() error {
	Vprintf("\t\tReading %s...\n", FileCoarseFasta)
	timer := time.Now()

	fastaReader := fasta.NewReader(coarsedb.FileFasta)
	for i := 0; true; i++ {
		seq, err := fastaReader.Read()
		if err == io.EOF {
			break
		}
		if err != nil {
			return err
		}
		coarsedb.Seqs = append(coarsedb.Seqs, NewFastaCoarseSeq(i, seq))
	}
	coarsedb.seqsRead = len(coarsedb.Seqs)

	Vprintf("\t\tDone reading %s (%s).\n", FileCoarseFasta, time.Since(timer))
	return nil
}

func (coarsedb *CoarseDB) saveFasta() (err error) {
	Vprintf("Writing %s...\n", FileCoarseFasta)
	Vprintf("Writing %s...\n", FileCoarseFastaIndex)
	timer := time.Now()

	byteOff := int64(0)
	buf := new(bytes.Buffer)

	if coarsedb.fastaIndexSize > 0 {
		info, err := coarsedb.FileFasta.Stat()
		if err != nil {
			return err
		}
		byteOff = info.Size()
	}

	for i := coarsedb.seqsRead; i < len(coarsedb.Seqs); i++ {
		buf.Reset()
		fmt.Fprintf(buf, "> %d\n%s\n", i, string(coarsedb.Seqs[i].Residues))

		if _, err = coarsedb.FileFasta.Write(buf.Bytes()); err != nil {
			return
		}
		err = binary.Write(coarsedb.FileFastaIndex, binary.BigEndian, byteOff)
		if err != nil {
			return
		}

		byteOff += int64(buf.Len())
	}

	Vprintf("Done writing %s (%s).\n", FileCoarseFasta, time.Since(timer))
	Vprintf("Done writing %s (%s).\n", FileCoarseFastaIndex, time.Since(timer))
	return nil
}

func (coarsedb *CoarseDB) readLinks() error {
	Vprintf("\t\tReading %s...\n", FileCoarseLinks)
	timer := time.Now()

	var cnt int32
	for coarseSeqId := 0; true; coarseSeqId++ {
		if binary.Read(coarsedb.FileLinks, binary.BigEndian, &cnt) != nil {
			break
		}
		for i := int32(0); i < cnt; i++ {
			newLink, err := coarsedb.readLink()
			if err != nil {
				return err
			}
			coarsedb.Seqs[coarseSeqId].addLink(newLink)
		}
	}

	Vprintf("\t\tDone reading %s (%s).\n", FileCoarseLinks, time.Since(timer))
	return nil
}

func (coarsedb *CoarseDB) readLink() (_ *LinkToCompressed, err error) {
	br := func(data interface{}) error {
		return binary.Read(coarsedb.FileLinks, binary.BigEndian, data)
	}

	var orgSeqId uint32
	var coarseStart, coarseEnd uint16

	if err = br(&orgSeqId); err != nil {
		return
	}
	if err = br(&coarseStart); err != nil {
		return
	}
	if err = br(&coarseEnd); err != nil {
		return
	}
	return NewLinkToCompressed(orgSeqId, coarseStart, coarseEnd), nil
}

func (coarsedb *CoarseDB) saveLinks() (err error) {
	Vprintf("Writing %s...\n", FileCoarseLinks)
	Vprintf("Writing %s...\n", FileCoarseLinksIndex)
	timer := time.Now()

	byteOff := int64(0)
	buf := new(bytes.Buffer)

	bw := func(data interface{}) error {
		return binary.Write(buf, binary.BigEndian, data)
	}
	for _, seq := range coarsedb.Seqs {
		// Reset the buffer so it's empty. We want it to only contain
		// the next set of links.
		buf.Reset()

		// Get the number of links.
		// This is necessary for reading so we know how many links to read.
		cnt := int32(0)
		for link := seq.Links; link != nil; link = link.Next {
			cnt++
		}
		if err = bw(cnt); err != nil {
			return
		}

		// For each link, write its original sequence id, and the coarse
		// sequence start/end.
		for link := seq.Links; link != nil; link = link.Next {
			if err = bw(link.OrgSeqId); err != nil {
				return
			}
			if err = bw(link.CoarseStart); err != nil {
				return
			}
			if err = bw(link.CoarseEnd); err != nil {
				return
			}
		}

		// Write the bytes to the links file.
		if _, err = coarsedb.FileLinks.Write(buf.Bytes()); err != nil {
			return
		}

		// Now write the byte offset that points to the start of this
		// set of links.
		err = binary.Write(coarsedb.FileLinksIndex, binary.BigEndian, byteOff)
		if err != nil {
			return
		}

		// Set the byte offset to be at the end of this set of links.
		byteOff += int64(buf.Len())
	}

	Vprintf("Done writing %s (%s).\n", FileCoarseLinks, time.Since(timer))
	Vprintf("Done writing %s (%s).\n", FileCoarseLinksIndex, time.Since(timer))
	return nil
}

func (coarsedb *CoarseDB) saveLinksPlain() error {
	Vprintf("Writing %s...\n", FileCoarsePlainLinks)
	timer := time.Now()

	csvWriter := csv.NewWriter(coarsedb.plainLinks)
	record := make([]string, 0, 10)
	for _, seq := range coarsedb.Seqs {
		record = record[:0]
		for link := seq.Links; link != nil; link = link.Next {
			record = append(record,
				fmt.Sprintf("%d", link.OrgSeqId),
				fmt.Sprintf("%d", link.CoarseStart),
				fmt.Sprintf("%d", link.CoarseEnd))
		}
		if err := csvWriter.Write(record); err != nil {
			return err
		}
	}
	csvWriter.Flush()

	Vprintf("Done writing %s (%s).\n", FileCoarsePlainLinks, time.Since(timer))
	return nil
}

func (comdb *CompressedDB) ReadSeq(
	coarsedb *CoarseDB, orgSeqId int) (OriginalSeq, error) {

	off, err := comdb.orgSeqOffset(orgSeqId)
	if err != nil {
		return OriginalSeq{}, err
	}

	newOff, err := comdb.File.Seek(off, os.SEEK_SET)
	if err != nil {
		return OriginalSeq{}, err
	} else if newOff != off {
		return OriginalSeq{},
			fmt.Errorf("Tried to seek to offset %d in the compressed "+
				"database, but seeked to %d instead.", off, newOff)
	}
	return comdb.ReadNextSeq(coarsedb, orgSeqId)
}

func (comdb *CompressedDB) ReadNextSeq(
	coarsedb *CoarseDB, orgSeqId int) (OriginalSeq, error) {

	csvReader := csv.NewReader(comdb.File)
	csvReader.LazyQuotes = true
	csvReader.Comma = ','
	csvReader.FieldsPerRecord = -1

	record, err := csvReader.Read()
	if err == io.EOF && len(record) == 0 {
		return OriginalSeq{}, fmt.Errorf("[csv reader]: id out of range")
	} else if err != nil && err != io.EOF {
		return OriginalSeq{}, fmt.Errorf("[csv reader]: %s", err)
	}

	cseq, err := readCompressedSeq(orgSeqId, record)
	if err != nil {
		return OriginalSeq{}, err
	}
	return cseq.Decompress(coarsedb)
}

func readCompressedSeq(id int, record []string) (CompressedSeq, error) {
	cseq := CompressedSeq{
		Id:    id,
		Name:  string([]byte(record[0])),
		Links: make([]LinkToCoarse, (len(record)-1)/4),
	}

	for i := 1; i < len(record); i += 4 {
		coarseSeqId64, err := strconv.Atoi(record[i+0])
		if err != nil {
			return CompressedSeq{}, nil
		}
		coarseStart64, err := strconv.Atoi(record[i+1])
		if err != nil {
			return CompressedSeq{}, nil
		}
		coarseEnd64, err := strconv.Atoi(record[i+2])
		if err != nil {
			return CompressedSeq{}, nil
		}
		lk := NewLinkToCoarseNoDiff(
			uint(coarseSeqId64), uint(coarseStart64), uint(coarseEnd64))

		cseq.Add(lk)
	}
	return cseq, nil
}

func (comdb *CompressedDB) orgSeqOffset(id int) (seqOff int64, err error) {
	tryOff := int64(id) * 8
	realOff, err := comdb.Index.Seek(tryOff, os.SEEK_SET)
	if err != nil {
		return 0, err
	} else if tryOff != realOff {
		return 0,
			fmt.Errorf("Tried to seek to offset %d in the compressed index, "+
				"but seeked to %d instead.", tryOff, realOff)
	}
	err = binary.Read(comdb.Index, binary.BigEndian, &seqOff)
	return
}

func nextSeqToWrite(
	nextIndex int, saved []CompressedSeq) (*CompressedSeq, []CompressedSeq) {

	for i, cseq := range saved {
		if cseq.Id < nextIndex {
			panic(fmt.Sprintf("Cannot keep sequences (%d) earlier than the "+
				"next index (%d)\n\n%s\n", cseq.Id, nextIndex, cseq))
		}
		if cseq.Id == nextIndex {
			cseq_ := cseq
			saved = append(saved[:i], saved[i+1:]...)
			return &cseq_, saved
		}
	}
	return nil, saved
}

func (comdb *CompressedDB) writer() {
	var record []string
	var err error
	var cseq *CompressedSeq

	byteOffset := int64(0)
	buf := new(bytes.Buffer)
	csvWriter := csv.NewWriter(buf)
	csvWriter.Comma = ','
	csvWriter.UseCRLF = false

	saved := make([]CompressedSeq, 0, 1000)
	nextIndex := comdb.NumSequences()

	// If we're appending to the index, set the byteOffset to be at the end
	// of the current compressed database.
	if comdb.indexSize > 0 {
		info, err := comdb.File.Stat()
		if err != nil {
			fmt.Fprintf(os.Stderr, "%s\n", err)
			os.Exit(1)
		}
		byteOffset = info.Size()
	}

	for possible := range comdb.writerChan {
		// We have to preserve the order of compressed sequences, so we don't
		// write anything until we have the next sequence that we expect.
		if possible.Id < nextIndex {
			panic(fmt.Sprintf("BUG: Next sequence expected is '%d', but "+
				"we have an earlier sequence: %d", nextIndex, possible.Id))
		}
		saved = append(saved, possible)

		cseq, saved = nextSeqToWrite(nextIndex, saved)
		for cseq != nil {
			// Reset the buffer so it's empty. We want it to only contain
			// the next record we're writing.
			buf.Reset()

			// Allocate memory for creating the next record.
			// A record is a sequence name followed by four-tuples of links:
			// (coarse-seq-id, coarse-start, coarse-end, diff).
			record = make([]string, 0, 1+4*len(cseq.Links))
			record = append(record, cseq.Name)
			for _, link := range cseq.Links {
				record = append(record,
					fmt.Sprintf("%d", link.CoarseSeqId),
					fmt.Sprintf("%d", link.CoarseStart),
					fmt.Sprintf("%d", link.CoarseEnd),
					link.OrigSeq)
			}

			// Write the record to our *buffer* and flush it.
			if err = csvWriter.Write(record); err != nil {
				fmt.Fprintf(os.Stderr, "%s\n", err)
				os.Exit(1)
			}
			csvWriter.Flush()

			// Pass the bytes on to the compressed file.
			if _, err = comdb.File.Write(buf.Bytes()); err != nil {
				fmt.Fprintf(os.Stderr, "%s\n", err)
				os.Exit(1)
			}

			// Now write the byte offset that points to the start of this record
			err = binary.Write(comdb.Index, binary.BigEndian, byteOffset)
			if err != nil {
				fmt.Fprintf(os.Stderr, "%s\n", err)
				os.Exit(1)
			}

			// Increment the byte offset to be at the end of this record.
			byteOffset += int64(buf.Len())

			nextIndex++
			cseq, saved = nextSeqToWrite(nextIndex, saved)
		}
	}
	comdb.Index.Close()
	comdb.File.Close()
	comdb.writerDone <- struct{}{}
}
