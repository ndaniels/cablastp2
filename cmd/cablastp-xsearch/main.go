package main

import (
	"bytes"
	"encoding/xml"
	"flag"
	"fmt"
	"io"
	"io/ioutil"
	"log"
	"os"
	"os/exec"
	"path"
	"runtime"
	"runtime/pprof"
	// "strings"

	"github.com/TuftsBCB/io/fasta"
	"github.com/TuftsBCB/seq"

	"github.com/ndaniels/cablastp2"
)

// A BLAST database is created on the reference sequence after compression.
// The search program will blast the input query sequences against this
// database (with a relaxed e-value), and expand the hits using links into an
// in memory FASTA file. This FASTA file is passed to the stdin of a
// `makeblastdb` command, which outputs the fine BLAST database. Finally, the
// query sequences are blasted against this new database, and the hits are
// returned unmodified.

var (
	// A default configuration.
	dbConf = cablastp.DefaultDBConf

	// Flags that affect the operation of search.
	// Flags that control algorithmic parameters are stored in `dbConf`.
	flagMakeBlastDB    = "makeblastdb"
	flagBlastp         = "blastp"
	flagBlastn         = "blastn"
	flagGoMaxProcs     = runtime.NumCPU()
	flagQuiet          = false
	flagCpuProfile     = ""
	flagMemProfile     = ""
	flagCoarseEval     = 5.0
	flagNoCleanup      = false
	flagCompressQuery  = false
	flagBatchQueries   = false
	flagIterativeQuery = false
)

// blastArgs are all the arguments after "--blast-args".
var blastArgs []string

func init() {
	log.SetFlags(0)

	flag.StringVar(&flagMakeBlastDB, "makeblastdb",
		flagMakeBlastDB,
		"The location of the 'makeblastdb' executable.")
	flag.StringVar(&flagBlastp, "blastp",
		flagBlastp,
		"The location of the 'blastp' executable.")
	flag.StringVar(&flagBlastn, "blastn",
		flagBlastn,
		"The location of the 'blastn' executable.")
	flag.Float64Var(&flagCoarseEval, "coarse-eval", flagCoarseEval,
		"The e-value threshold for the coarse search. This will NOT\n"+
			"\tbe used on the fine search. The fine search e-value threshold\n"+
			"\tcan be set in the 'blast-args' argument.")
	flag.BoolVar(&flagNoCleanup, "no-cleanup", flagNoCleanup,
		"When set, the temporary fine BLAST database that is created\n"+
			"\twill NOT be deleted.")

	flag.IntVar(&flagGoMaxProcs, "p", flagGoMaxProcs,
		"The maximum number of CPUs that can be executing simultaneously.")
	flag.BoolVar(&flagQuiet, "quiet", flagQuiet,
		"When set, the only outputs will be errors echoed to stderr.")
	flag.StringVar(&flagCpuProfile, "cpuprofile", flagCpuProfile,
		"When set, a CPU profile will be written to the file specified.")
	flag.StringVar(&flagMemProfile, "memprofile", flagMemProfile,
		"When set, a memory profile will be written to the file specified.")
	flag.BoolVar(&flagIterativeQuery, "iterative-queries", flagIterativeQuery,
		"When set, will process queries one at a time instead of as a batch.")

	// find '--blast-args' and chop off the remainder before letting the flag
	// package have its way.
	for i, arg := range os.Args {
		if arg == "--blast-args" {
			blastArgs = os.Args[i+1:]
			os.Args = os.Args[:i]
		}
	}

	flag.Usage = usage
	flag.Parse()

	runtime.GOMAXPROCS(flagGoMaxProcs)
}

func main() {
	queryBuf := new(bytes.Buffer)  // might need more than 1 buffer
	searchBuf := new(bytes.Buffer) // might need more than 1 buffer

	if flag.NArg() != 2 {
		flag.Usage()
	}

	// If the quiet flag isn't set, enable verbose output.
	if !flagQuiet {
		cablastp.Verbose = true
	}

	inputFastaQuery, err := getInputFasta()
	if err != nil {
		fatalf("Could not read input fasta query: %s\n", err)
	}

	db, err := cablastp.NewReadDB(flag.Arg(0))
	if err != nil {
		fatalf("Could not open '%s' database: %s\n", flag.Arg(0), err)
	}
	// For query-compression mode, we first run compression on the query file
	// then coarse-coarse search, decompress both, fine-fine search.
	// otherwise, just coarse search, decompress results, fine search.
	// iterate over the query sequences in the input fasta
	// initially, only implement standard search.

	cablastp.Vprintln("\nProcessing Queries one at a time...")
	f := fasta.NewWriter(queryBuf)
	reader := fasta.NewReader(inputFastaQuery)
	for i := 0; true; i++ {
		sequence, err := reader.Read()
		if err == io.EOF {
			break
		}
		if err != nil {
			fatalf("Could not read input fasta query: %s\n", err)
		}
		origSeq := sequence.Bytes()
		n := sequence.Name
		// generate 6 ORFs
		transSeqs := cablastp.Translate(origSeq)
		for _, s := range transSeqs {
			// reduce each one
			result := seq.NewSequenceString(n, string(cablastp.Reduce(s)))
			f.Write(result)
		}

		if flagIterativeQuery {
			transQueries := bytes.NewReader(queryBuf.Bytes())
			processQueries(db, transQueries, searchBuf)
			queryBuf.Reset()
		}
	}

	if !flagIterativeQuery {
		transQueries := bytes.NewReader(queryBuf.Bytes())
		processQueries(db, transQueries, searchBuf)
	}

	cleanup(db)
}

func processQueries(
	db *cablastp.DB, transQueries *bytes.Reader, searchBuf *bytes.Buffer) error {
	// now we will read from queryBuf!
	// I think we create a NewReader from queryBuf?
	// this now needs to become the replacement for inputFastaQuery
	// so must use a different buffer for that.
	// we need a buffer for the query trans/reduce
	// and a buffer for coarse blast results

	cablastp.Vprintln("\nBlasting query on coarse database...")
	if err := blastCoarse(db, transQueries, searchBuf); err != nil {
		fatalf("Error blasting coarse database: %s\n", err)
	}

	cablastp.Vprintln("Decompressing blast hits...")
	expandedSequences, err := expandBlastHits(db, searchBuf)
	if err != nil {
		fatalf("%s\n", err)
	}

	// Write the contents of the expanded sequences to a fasta file.
	// It is then indexed using makeblastdb.
	searchBuf.Reset()
	if err := writeFasta(expandedSequences, searchBuf); err != nil {
		fatalf("Could not create FASTA input from coarse hits: %s\n", err)
	}

	// Create the fine blast db in a temporary directory
	cablastp.Vprintln("Building fine BLAST database...")
	tmpDir, err := makeFineBlastDB(db, searchBuf)
	if err != nil {
		fatalf("Could not create fine database to search on: %s\n", err)
	}

	// retrieve the cluster members for the original representative query seq

	// pass them to blastx on the expanded (fine) db

	// Finally, run the query against the fine fasta database and pass on the
	// stdout and stderr...
	cablastp.Vprintln("Blasting query on fine database...")
	if _, err := inputFastaQuery.Seek(0, os.SEEK_SET); err != nil {
		fatalf("Could not seek to start of query fasta input: %s\n", err)
	}
	if err := blastFine(db, tmpDir, inputFastaQuery); err != nil {
		fatalf("Error blasting fine database: %s\n", err)
	}

	// Delete the temporary fine database.
	if !flagNoCleanup {
		if err := os.RemoveAll(tmpDir); err != nil {
			fatalf("Could not delete fine BLAST database: %s\n", err)
		}
	}
}

func s(i int) string {
	return fmt.Sprintf("%d", i)
}

func su(i uint64) string {
	return fmt.Sprintf("%d", i)
}

func blastFine(
	db *cablastp.DB, blastFineDir string, stdin *bytes.Reader) error {

	// We pass our own "-db" flag to blastp, but the rest come from user
	// defined flags.
	flags := []string{
		"-db", path.Join(blastFineDir, cablastp.FileBlastFine),
		"-dbsize", su(db.BlastDBSize),
		"-num_threads", s(flagGoMaxProcs),
	}
	flags = append(flags, blastArgs...)

	cmd := exec.Command(flagBlastp, flags...)
	cmd.Stdin = stdin
	cmd.Stdout = os.Stdout
	cmd.Stderr = os.Stderr
	return cablastp.Exec(cmd)
}

func makeFineBlastDB(db *cablastp.DB, stdin *bytes.Buffer) (string, error) {
	tmpDir, err := ioutil.TempDir("", "cablastp-fine-search-db")
	if err != nil {
		return "", fmt.Errorf("Could not create temporary directory: %s\n", err)
	}

	cmd := exec.Command(
		flagMakeBlastDB, "-dbtype", "prot",
		"-title", cablastp.FileBlastFine,
		"-in", "-",
		"-out", path.Join(tmpDir, cablastp.FileBlastFine))
	cmd.Stdin = stdin

	cablastp.Vprintf("Created temporary fine BLAST database in %s\n", tmpDir)

	return tmpDir, cablastp.Exec(cmd)
}

func writeFasta(oseqs []cablastp.OriginalSeq, buf *bytes.Buffer) error {
	for _, oseq := range oseqs {
		_, err := fmt.Fprintf(buf, "> %s\n%s\n",
			oseq.Name, string(oseq.Residues))
		if err != nil {
			return fmt.Errorf("Could not write to buffer: %s", err)
		}
	}
	return nil
}

func expandBlastHits(
	db *cablastp.DB, blastOut *bytes.Buffer) ([]cablastp.OriginalSeq, error) {

	results := blast{}
	if err := xml.NewDecoder(blastOut).Decode(&results); err != nil {
		return nil, fmt.Errorf("Could not parse BLAST search results: %s", err)
	}

	used := make(map[int]bool, 100) // prevent original sequence duplicates
	oseqs := make([]cablastp.OriginalSeq, 0, 100)
	for _, hit := range results.Hits {
		for _, hsp := range hit.Hsps {
			someOseqs, err := db.CoarseDB.Expand(db.ComDB,
				hit.Accession, hsp.HitFrom, hsp.HitTo)
			if err != nil {
				errorf("Could not decompress coarse sequence %d (%d, %d): %s\n",
					hit.Accession, hsp.HitFrom, hsp.HitTo, err)
				continue
			}

			// Make sure this hit is below the coarse e-value threshold.
			if hsp.Evalue > flagCoarseEval {
				continue
			}

			for _, oseq := range someOseqs {
				if used[oseq.Id] {
					continue
				}
				used[oseq.Id] = true
				oseqs = append(oseqs, oseq)
			}
		}
	}
	return oseqs, nil
}

func blastCoarse(
	db *cablastp.DB, stdin *bytes.Reader, stdout *bytes.Buffer) error {

	cmd := exec.Command(
		flagBlastn,
		"-db", path.Join(db.Path, cablastp.FileBlastCoarse),
		"-num_threads", s(flagGoMaxProcs),
		"-outfmt", "5", "-dbsize", su(db.BlastDBSize))
	cmd.Stdin = stdin
	cmd.Stdout = stdout
	return cablastp.Exec(cmd)
}

func getInputFasta() (*bytes.Reader, error) {
	queryFasta, err := os.Open(flag.Arg(1))
	if err != nil {
		return nil, fmt.Errorf("Could not open '%s': %s.", flag.Arg(1), err)
	}
	bs, err := ioutil.ReadAll(queryFasta)
	if err != nil {
		return nil, fmt.Errorf("Could not read input fasta query: %s", err)
	}
	return bytes.NewReader(bs), nil
}

func cleanup(db *cablastp.DB) {
	if len(flagCpuProfile) > 0 {
		pprof.StopCPUProfile()
	}
	if len(flagMemProfile) > 0 {
		writeMemProfile(fmt.Sprintf("%s.last", flagMemProfile))
	}
	db.ReadClose()
}

func fatalf(format string, v ...interface{}) {
	fmt.Fprintf(os.Stderr, format, v...)
	os.Exit(1)
}

func errorf(format string, v ...interface{}) {
	fmt.Fprintf(os.Stderr, format, v...)
}

func writeMemProfile(name string) {
	f, err := os.Create(name)
	if err != nil {
		fatalf("%s\n", err)
	}
	pprof.WriteHeapProfile(f)
	f.Close()
}

func usage() {
	fmt.Fprintf(os.Stderr,
		"\nUsage: %s [flags] database-directory query-fasta-file "+
			"[--blast-args BLASTP_ARGUMENTS]\n",
		path.Base(os.Args[0]))
	cablastp.PrintFlagDefaults()
	os.Exit(1)
}
