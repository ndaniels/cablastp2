

//Compressing Sequences

func compressQueries( queryFileName string) (*DB, error) {
	queryDB, err := cablastp.NewWriteDB(dbConf, queryFileName)
	pool := startCompressWorkers( queryDB)
	orgSeqId := queryDB.ComDb.NumSequences()
	for _, arg := range flag.Args()[1:]{
		seqChan, err := cablastp.ReadOriginalSeqs(arg, ignoredResidues)
		for readSeq := range seqChan {
			orgSeqId = pool.compress(orgSeqId, readSeq.Seq)
		}
	}
	cleanup(queryDB, &pool)
	return &queryDB, nil
}



func normalSearch(){
	inputFastaQuery, err := getInputFasta()
	db, err := cablastp.NewReadDB(flag.Arg(0))
	reducedFasta, err := cablastp.ReduceQuerySeqs(inputFastaQuery)
	blastCoarse(db, reducedFasta, buf)
	expandedSequences, err := expandBlastHits(db, buf)
	writeFasta(expandedSequences, buf)
	tmpDir, err := makeFineBlastDB(db, buf)
	blastFine(db, tmpDir, inputFastaQuery)
	cleanup(db)
}

func queryCompSearch( databaseName string, queryDB *DB){
	db, err := cablastp.NewReadDB(databaseName)
	// TODO: How to get all compressed query sequences
	for _, aCompressedQuerySequence := range queryDB.allQuerySequences(){
		reader := fasta.NewReader(queryDb.coarseFasta)
		blastCoarse(db, reader, buf)
		expandedSequences, err := expandBlastHits(db, buf)
		writeFasta(expandedSequences, buf)
		tmpDir, err := makeFineBlastDB(db, buf)
		// TODO: How to get matches for a compressed seq?
		expandedQueries := queryDB.expandCompressedSeq( aCompressedQuerySeq)
		for _, anExpandedQuerySequence := range expandedQueries{
			blastFine( blastDB, anExpandedQuerySequence, buf)
		}
	}
cleanup(db)
}