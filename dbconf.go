package cablastp

import (
	"flag"
	"fmt"
	"io"

	"github.com/BurntSushi/toml"
)

type DBConf struct {
	MinMatchLen         int
	MatchKmerSize       int
	GappedWindowSize    int
	UngappedWindowSize  int
	ExtSeqIdThreshold   int
	MatchSeqIdThreshold int
	MatchExtend         int
	MapSeedSize         int
	ExtSeedSize         int
	LowComplexity       int
	SeedLowComplexity   int
	SavePlain           bool
	ReadOnly            bool
	BlastMakeBlastDB    string
	BlastDBSize         uint64
}

var DefaultDBConf = &DBConf{
	MinMatchLen:         40,
	MatchKmerSize:       4,
	GappedWindowSize:    25,
	UngappedWindowSize:  10,
	ExtSeqIdThreshold:   60,
	MatchSeqIdThreshold: 70,
	MatchExtend:         30,
	MapSeedSize:         6,
	ExtSeedSize:         0,
	LowComplexity:       10,
	SeedLowComplexity:   6,
	SavePlain:           false,
	ReadOnly:            true,
	BlastMakeBlastDB:    "makeblastdb",
	BlastDBSize:         0,
}

func LoadDBConf(r io.Reader) (conf *DBConf, err error) {
	defer func() {
		if perr := recover(); perr != nil {
			err = perr.(error)
		}
	}()
	conf = DefaultDBConf

	if _, err := toml.DecodeReader(r, &conf); err != nil {
		return nil, err
	}

	return conf, nil
}

func (flagConf *DBConf) FlagMerge(fileConf *DBConf) (*DBConf, error) {
	only := make(map[string]bool, 0)
	flag.Visit(func(f *flag.Flag) { only[f.Name] = true })

	if only["map-seed-size"] {
		return flagConf, fmt.Errorf("The map seed size cannot be changed for " +
			"an existing database.")
	}
	if only["read-only"] {
		return flagConf, fmt.Errorf("The read-only setting cannot be changed " +
			"for an existing database.")
	}

	if !only["min-match-len"] {
		flagConf.MinMatchLen = fileConf.MinMatchLen
	}
	if !only["match-kmer-size"] {
		flagConf.MatchKmerSize = fileConf.MatchKmerSize
	}
	if !only["gapped-window-size"] {
		flagConf.GappedWindowSize = fileConf.GappedWindowSize
	}
	if !only["ungapped-window-size"] {
		flagConf.UngappedWindowSize = fileConf.UngappedWindowSize
	}
	if !only["ext-seq-id-threshold"] {
		flagConf.ExtSeqIdThreshold = fileConf.ExtSeqIdThreshold
	}
	if !only["match-seq-id-threshold"] {
		flagConf.MatchSeqIdThreshold = fileConf.MatchSeqIdThreshold
	}
	if !only["match-extend"] {
		flagConf.MatchExtend = fileConf.MatchExtend
	}
	if !only["ext-seed-size"] {
		flagConf.ExtSeedSize = fileConf.ExtSeedSize
	}
	if !only["low-complexity"] {
		flagConf.LowComplexity = fileConf.LowComplexity
	}
	if !only["seed-low-complexity"] {
		flagConf.SeedLowComplexity = fileConf.SeedLowComplexity
	}
	if !only["plain"] {
		flagConf.SavePlain = fileConf.SavePlain
	}
	if !only["read-only"] {
		flagConf.ReadOnly = fileConf.ReadOnly
	}
	if !only["makeblastdb"] {
		flagConf.BlastMakeBlastDB = fileConf.BlastMakeBlastDB
	}
	if !only["dbsize"] {
		flagConf.BlastDBSize = fileConf.BlastDBSize
	}
	return flagConf, nil
}

func (dbConf DBConf) Write(w io.Writer) error {
	encoder := toml.NewEncoder(w)
	if err := encoder.Encode(dbConf); err != nil {
		return err
	}
	return nil
}
